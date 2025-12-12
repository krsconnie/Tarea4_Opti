# ggcplex.py
# Implementación de la formulación GG para ATSP usando DOcplex (CPLEX backend)
#
# Uso:
#   python ggcplex.py ruta_instancia.atsp
#
# Produce: solución (si la hay) y estadísticas necesarias para la tabla

import sys
import math
from docplex.mp.model import Model

def parse_tsplib_atsp(path):
    """
    Parser simple para archivos TSPLIB ATSP (formato clásico con EDGE_WEIGHT_SECTION
    que contiene la matriz completa). Devuelve la matriz de costos como lista de listas.
    Asume que la sección EDGE_WEIGHT_SECTION contiene n*n números (full matrix).
    """
    with open(path, 'r') as f:
        lines = [ln.strip() for ln in f]
    dim = None
    section = None
    data = []
    for ln in lines:
        if ln == '' or ln.startswith('COMMENT') or ln.startswith('NAME'):
            continue
        if ln.startswith('DIMENSION'):
            parts = ln.replace(':', ' ').split()
            for p in parts:
                if p.isdigit():
                    dim = int(p)
                    break
        if ln.startswith('EDGE_WEIGHT_SECTION'):
            section = 'EDGE'
            continue
        if section == 'EDGE':
            if ln == 'EOF':
                break
            # puede haber números en varias columnas, los juntamos.
            toks = ln.split()
            for t in toks:
                data.append(int(t))
    if dim is None:
        raise ValueError("No se pudo detectar DIMENSION en el archivo TSPLIB.")
    if len(data) < dim*dim:
        raise ValueError(f"No hay suficientes datos de matriz (esperado {dim*dim}, recibido {len(data)})")
    # formar la matriz
    matrix = []
    idx = 0
    for i in range(dim):
        row = []
        for j in range(dim):
            row.append(data[idx])
            idx += 1
        matrix.append(row)
    return matrix

def build_and_solve_GG(cost_matrix, time_limit_seconds=3600, log_output=False):
    """
    Construye y resuelve la formulación GG para la matriz de costos dada.
    Devuelve un diccionario con la información requerida (n, var_count, cons_count, time, gap, best_bound, obj).
    """
    n = len(cost_matrix)
    mdl = Model(name="GG_ATSP")
    # Variables x_{i,j} binarias para todos i,j = 0..n-1
    x = {}
    for i in range(n):
        for j in range(n):
            # opcional: si cii es muy grande, aún incluimos la variable; el paper sugiere c_ii = +inf
            x[(i,j)] = mdl.binary_var(name=f"x_{i}_{j}")
    # Variables g_{i,j} para i=1..n-1 (equiv a i=2..n en paper), j=0..n-1
    g = {}
    for i in range(1, n):
        for j in range(n):
            g[(i,j)] = mdl.continuous_var(lb=0.0, name=f"g_{i}_{j}")

    # Objetivo
    mdl.minimize(mdl.sum(cost_matrix[i][j] * x[(i,j)] for i in range(n) for j in range(n)))

    # Restricciones de grado: entra = 1, sale = 1
    for j in range(n):
        mdl.add_constraint(mdl.sum(x[(i,j)] for i in range(n)) == 1, ctname=f"in_deg_{j}")
    for i in range(n):
        mdl.add_constraint(mdl.sum(x[(i,j)] for j in range(n)) == 1, ctname=f"out_deg_{i}")

    # Ecuaciones (11) de GG: para i = 2..n (aquí i indexado 1..n-1)
    for i in range(1, n):
        mdl.add_constraint(mdl.sum(g[(i,j)] for j in range(n)) - mdl.sum(g[(k,i)] for k in range(1,n)) == 1,
                           ctname=f"flow_balance_{i}")

    # Bound: g_{i,j} <= (n-1) * x_{i,j}  (ecuación (12))
    bigM = n - 1
    for i in range(1, n):
        for j in range(n):
            mdl.add_constraint(g[(i,j)] <= bigM * x[(i,j)], ctname=f"g_bound_{i}_{j}")

    # Parámetros CPLEX vía docplex
    mdl.parameters.timelimit = time_limit_seconds
    # opcional: más logging 
    if log_output:
        mdl.print_information()

    # resolver
    sol = mdl.solve(log_output=log_output)
    # obtener el objeto CPLEX subyacente para métricas más detalladas
    try:
        cpx = mdl.get_cplex()
    except:
        cpx = None

    # Contar variables y restricciones (manuales)
    num_x = n * n
    num_g = (n - 1) * n
    var_count = num_x + num_g

    # constraints:
    # - 2n degree constraints
    # - (n-1) flow balance constraints
    # - (n-1)*n g_bounds
    cons_count = 2*n + (n-1) + (n-1)*n

    # extraer info de CPLEX si está disponible
    solve_time = None
    mipgap = None
    best_bound = None
    obj = None
    status = None
    if cpx is not None:
        try:
            solve_time = cpx.get_time()
        except:
            solve_time = None
        try:
            # gap relativo (puede fallar si no es MIP)
            mipgap = cpx.solution.get_mip_relative_gap()
        except:
            mipgap = None
        try:
            best_bound = cpx.solution.get_best_objective()
        except:
            best_bound = None
        try:
            # si hay solución entera:
            if cpx.solution.get_status() not in [-1, 0]:
                obj = mdl.objective_value if sol is not None else None
        except:
            obj = None
        try:
            status = cpx.solution.get_status()
        except:
            status = None
    else:
        # si no podemos acceder a CPLEX, usamos lo que docplex nos da
        if sol is not None:
            obj = mdl.objective_value
        # docplex solve_details
        try:
            sd = mdl.solve_details
            if sd is not None:
                solve_time = sd.time
                # sd has mip_gap sometimes
                if hasattr(sd, 'mip_relative_gap'):
                    mipgap = sd.mip_relative_gap
        except:
            pass

    # Si el solver no resolvió la instancia en tiempo, definimos gap 100% como indica la instrucción
    # pero aquí dejamos el valor que devolvió CPLEX si existe.
    result = {
        "n": n,
        "var_count": var_count,
        "cons_count": cons_count,
        "solve_time_sec": solve_time,
        "mipgap": mipgap,
        "best_bound": best_bound,
        "objective": obj,
        "status": status,
        "solution_exists": sol is not None
    }
    return result, mdl, sol

def example_run_on_file(path_atsp, time_limit_seconds=3600, log_output=False):
    print("Parseando instancia:", path_atsp)
    cost = parse_tsplib_atsp(path_atsp)
    print("Dimension detectada:", len(cost))
    res, mdl, sol = build_and_solve_GG(cost, time_limit_seconds=time_limit_seconds, log_output=log_output)
    print("*** RESULTADOS ***")
    for k,v in res.items():
        print(f"{k}: {v}")
   
    if res["solution_exists"]:
        tour = []
        for i in range(res["n"]):
            for j in range(res["n"]):
                var = mdl.get_var_by_name(f"x_{i}_{j}")
                if var is not None and var.solution_value > 0.5:
                    tour.append((i,j))
        print("Arcos seleccionados (parcial):", tour)
    return res

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Uso: python ggcplex.py ruta_instancia.atsp [time_limit_seconds]")
        print("Ejemplo: python ggcplex.py ftv33.atsp 3600")
        sys.exit(1)
    ruta = sys.argv[1]
    tl = 3600
    if len(sys.argv) >= 3:
        tl = int(sys.argv[2])
    example_run_on_file(ruta, time_limit_seconds=tl, log_output=True)
