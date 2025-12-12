import re
import json
import time
import math
import cplex
from pathlib import Path
from docplex.mp.model import Model

###############################################################################
# CONFIG
###############################################################################

# Definimos la ruta base como el directorio padre'
# Esto nos lleva del script al Directorio Raíz.
BASE_DIR = Path(__file__).resolve().parent.parent 

# Define la ruta del directorio de entrada: BASE_DIR / "instancias"
INPUT_DIR = BASE_DIR / "instancias"
# Define la ruta del directorio de salida: BASE_DIR / "Resultados"
OUTPUT_DIR = BASE_DIR / "Resultados"

print("Usando solver: CPLEX (docplex)")

###############################################################################
# PARSER TSPLIB
###############################################################################

def parse_matrix_file(file_path):
    """Lee matriz TSPLIB full matrix o matriz simple."""
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    # Detectar líneas de solo números (matriz simple)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    numeric = [ln for ln in lines if re.match(r"^[\d\.\-\s]+$", ln)]

    if numeric:
        rows = [re.split(r"\s+", ln) for ln in numeric]
        n = len(rows)
        if all(len(r) == n for r in rows):
            try:
                M = [[float(x) for x in row] for row in rows]
                for i in range(n):
                    # Asignar un costo muy alto a los viajes i -> i
                    M[i][i] = 1e6 
                return M
            except:
                pass

    # TSPLIB FULL MATRIX
    m = re.search(r"DIMENSION\s*:\s*(\d+)", text, re.IGNORECASE)
    if m:
        n = int(m.group(1))
        body = re.search(r"EDGE_WEIGHT_SECTION\s*(.*)", text, re.IGNORECASE | re.DOTALL)
        if body:
            nums = re.findall(r"[-+]?\d*\.\d+|\d+", body.group(1))
            if len(nums) >= n * n:
                nums = nums[:n*n]
                M = [list(map(float, nums[i*n:(i+1)*n])) for i in range(n)]
                for i in range(n):
                    # Asignar un costo muy alto a los viajes i -> i
                    M[i][i] = 1e6 
                return M

    raise ValueError(f"No se pudo interpretar {file_path}")

###############################################################################
# CONSTRUIR MODELO MTZ
###############################################################################

def build_MTZ_model(matrix, bounded=True):
    """Devuelve un modelo CPLEX MTZ (formulación de Miller-Tucker-Zemlin)."""
    n = len(matrix)
    bigM = n - 1 # El valor de Big M en la formulación MTZ
    nodes = range(n)

    mdl = Model(name=f"MTZ_{'bounded' if bounded else 'unbounded'}")

    # Variables de decisión binarias x_ij (1 si se viaja de i a j, 0 en caso contrario)
    x = {(i,j): mdl.binary_var(f"x_{i}_{j}") for i in nodes for j in nodes}

    # Variables continuas u_i (posiciones en la ruta)
    if bounded:
        # u_i acotadas entre 1 y n-1 para i > 0
        u = {i: mdl.continuous_var(lb=1, ub=n-1, name=f"u_{i}") for i in nodes if i != 0}
    else:
        # u_i no acotadas (se asume que u_0 = 0 implícitamente)
        u = {i: mdl.continuous_var(name=f"u_{i}") for i in nodes if i != 0}

    # Objetivo: Minimizar el costo total
    mdl.minimize(mdl.sum(matrix[i][j] * x[(i,j)] for i in nodes for j in nodes))

    # Restricciones de grado de entrada (entrar a cada nodo una vez)
    for j in nodes:
        mdl.add_constraint(mdl.sum(x[(i,j)] for i in nodes) == 1, ctname=f"entrada_{j}")
        
    # Restricciones de grado de salida (salir de cada nodo una vez)
    for i in nodes:
        mdl.add_constraint(mdl.sum(x[(i,j)] for j in nodes) == 1, ctname=f"salida_{i}")

    # Restricciones de eliminación de subrutas (MTZ)
    for i in nodes:
        if i == 0: continue # La variable u_0 no existe (es 0)
        for j in nodes:
            if j == 0 or j == i: continue # No para el nodo 0 o cuando i=j
            # u_i - u_j + (n-1) * x_ij <= n - 2
            mdl.add_constraint(u[i] - u[j] + bigM * x[(i,j)] <= n - 2, ctname=f"mtz_{i}_{j}")

    return mdl


###############################################################################
# EXTRACCIÓN DE MÉTRICAS DE CPLEX
###############################################################################

def get_stats_docplex(model, t0):
    """Extrae métricas desde CPLEX vía docplex."""
    # Intentar resolver el modelo
    sol = model.solve() 

    elapsed = time.time() - t0

    # Si no hay solución, devolver métricas incompletas
    if sol is None:
        return {
            "objetivo": None,
            "variables": model.number_of_variables,
            "restricciones": model.number_of_constraints,
            "tiempo": elapsed,
            "gap": 100.0,
            "best_bound": None
        }

    details = model.solve_details

    # El valor objetivo de la solución encontrada
    objective_value = sol.get_objective_value()

    # Gap correcto (en docplex se llama mip_relative_gap)
    try:
        gap = details.mip_relative_gap
        if gap is None:
            gap = 1.0 # 100% gap si no se conoce
    except:
        gap = 1.0

    gap = gap * 100  # convertir a porcentaje

    try:
        best_bound = details.best_bound
    except:
        best_bound = None

    return {
        "objetivo": objective_value,
        "variables": model.number_of_variables,
        "restricciones": model.number_of_constraints,
        "tiempo": elapsed,
        "gap": gap,
        "best_bound": best_bound
    }

###############################################################################
# SOLVER GENERAL PARA UNA INSTANCIA
###############################################################################

def solve_instance(matrix, time_limit=60):
    """Resuelve MTZ bounded y unbounded y muestra métricas en consola."""
    n = len(matrix)

    out = {
        "nodos": n,
        "modelos": {}
    }

    for bounded_flag in [True, False]:
        name = "MTZ_bounded" if bounded_flag else "MTZ_unbounded"

        print(f"\n--- {name} ---")  # encabezado en consola

        model = build_MTZ_model(matrix, bounded=bounded_flag)
        model.parameters.timelimit = time_limit

        t0 = time.time()
        stats = get_stats_docplex(model, t0)

        out["modelos"][name] = stats

        # Impresión bonita en consola
        print(f"Nodos:                 {n}")
        print(f"Variables:             {stats['variables']}")
        print(f"Restricciones:         {stats['restricciones']}")
        print(f"Objetivo:              {stats['objetivo']}")
        print(f"Tiempo (s):            {stats['tiempo']:.3f}")
        print(f"Gap (%):               {stats['gap']:.2f}")
        print(f"Best bound:            {stats['best_bound']}")

    return out


###############################################################################
# FUNCIÓN PRINCIPAL DE EJECUCIÓN
###############################################################################

def main(time_limit=60):
    """Procesa el archivo br17.atsp, resuelve el modelo MTZ y guarda los resultados."""
    
    # Lista de archivos a procesar (SOLO br17.atsp)
    target_file_name = "br17.atsp"

    results = []

    # Intentar encontrar y procesar solo br17.atsp
    file_path = INPUT_DIR / target_file_name

    if not file_path.exists():
        print(f"\nERROR: No se encontró la instancia {target_file_name} en {INPUT_DIR}")
        # Salir si el archivo no existe
        exit(1) 

    print(f"\nProcesando {target_file_name}...")
    try:
        M = parse_matrix_file(file_path)
    except Exception as e:
        print("Error al leer:", e)
        results.append({"instance": target_file_name, "error": str(e)})
        # Salir después de un error de lectura
        exit(1)

    # Resolver la instancia con el límite de tiempo especificado
    stats = solve_instance(M, time_limit=time_limit) 
    stats["instance"] = target_file_name
    results.append(stats)

    # Guardar los resultados en el archivo de resumen
    output_summary_path = OUTPUT_DIR / "summary_br17.json" # Cambiado el nombre para ser específico
    with open(output_summary_path, "w") as f:
        json.dump(results, f, indent=2)

    print("\nProcesamiento terminado.")
    print("Resumen guardado en:", output_summary_path)


###############################################################################
# PUNTO DE ENTRADA DEL SCRIPT
###############################################################################

if __name__ == "__main__":
    main(time_limit=60)