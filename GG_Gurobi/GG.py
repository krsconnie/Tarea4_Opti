import sys
import glob
import time
import os
import threading
import pandas as pd
import gurobipy as gp
from gurobipy import GRB

def iniciar_reloj(stop_flag):
    inicio = time.time()
    while not stop_flag["stop"]:
        trans = int(time.time() - inicio)
        print(f"\r⏳ Tiempo transcurrido: {trans} s", end="", flush=True)
        time.sleep(0.5)
    print()  


def leer_archivo_tsplib(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    n = 0
    raw_data = []
    reading = False

    for line in lines:
        line = line.strip()
        if line.startswith("DIMENSION"):
            n = int(line.split(":")[1].strip())
        elif line.startswith("EDGE_WEIGHT_SECTION"):
            reading = True
        elif line.startswith("EOF"):
            break
        elif reading:
            vals = [int(x) for x in line.split() if x.isdigit()]
            raw_data.extend(vals)

    matriz = []
    for i in range(n):
        matriz.append(raw_data[i * n:(i + 1) * n])
    return n, matriz


# solver GG
def solve_atsp_gavish_graves(filename, n, dist, time_limit=3600):

    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()

    model = gp.Model("ATSP_GG", env=env)
    model.setParam("TimeLimit", time_limit)

    N = range(n)
    N2 = range(1, n)

    # variables
    x = model.addVars([(i, j) for i in N for j in N if i != j],
                      vtype=GRB.BINARY, name="x")

    g = model.addVars([(i, j) for i in N2 for j in N],
                      lb=0, ub=(n - 1),
                      vtype=GRB.CONTINUOUS,
                      name="g")

    # objetivo
    model.setObjective(
        gp.quicksum(dist[i][j] * x[i, j] for i, j in x.keys()),
        GRB.MINIMIZE
    )

    # restricciones
    for i in N:
        model.addConstr(gp.quicksum(x[i, j] for j in N if j != i) == 1)

    for j in N:
        model.addConstr(gp.quicksum(x[i, j] for i in N if i != j) == 1)

    for i in N2:
        model.addConstr(
            gp.quicksum(g[i, j] for j in N)
            - gp.quicksum(g[k, i] for k in N2)
            == 1
        )

    for i in N2:
        for j in N:
            if i != j:
                model.addConstr(g[i, j] <= (n - 1) * x[i, j])
            else:
                model.addConstr(g[i, j] == 0)

    # resolver con reloj
    stop_flag = {"stop": False}
    hilo = threading.Thread(target=iniciar_reloj, args=(stop_flag,))
    hilo.start()

    inicio = time.time()
    model.optimize()
    tiempo = time.time() - inicio

    stop_flag["stop"] = True
    hilo.join()

    # datos requeridos
    vars_total = model.NumVars
    restr_total = model.NumConstrs

    if model.SolCount == 0:
        return {
            "Instancia": os.path.basename(filename),
            "Nodos": n,
            "Vars": vars_total,
            "Restr": restr_total,
            "Tiempo (s)": round(tiempo, 4),
            "Gap (%)": 100.0,
            "Best Bound": model.ObjBound,
            "Objetivo": None
        }

    return {
        "Instancia": os.path.basename(filename),
        "Nodos": n,
        "Vars": vars_total,
        "Restr": restr_total,
        "Tiempo (s)": round(tiempo, 4),
        "Gap (%)": round(model.MIPGap * 100, 4),
        "Best Bound": model.ObjBound,
        "Objetivo": model.ObjVal
    }


# manual y automatico
if __name__ == "__main__":

    BASE = r"C:\Users\ZEPHYRUS G14\Desktop\Nueva carpeta\Nueva carpeta"

    carpetas = {
        "Pequeños": os.path.join(BASE, "Pequeños"),
        "Medianos": os.path.join(BASE, "Medianos"),
        "Grandes": os.path.join(BASE, "Grandes")
    }

    MODO_MANUAL = False
    INSTANCIA_MANUAL = "rbg403.atsp"

    resultados = []

    print("\n--- Algoritmo GG (Gavish & Graves) | Gurobi ---")

    # manual(prueba de 1 instancia))
    if MODO_MANUAL:
        archivo = None
        for grupo, carpeta in carpetas.items():
            ruta = os.path.join(carpeta, INSTANCIA_MANUAL)
            if os.path.exists(ruta):
                archivo = ruta
                grupo_usado = grupo
                break

        if not archivo:
            print("ERROR: No se encontró la instancia.")
            sys.exit()

        print(f"\nProcesando instancia manual: {archivo}")

        n, dist = leer_archivo_tsplib(archivo)
        res = solve_atsp_gavish_graves(archivo, n, dist)
        res["Grupo"] = grupo_usado
        resultados.append(res)

        print("\nRESULTADO:")
        print(res)

    # automatico (todas las instancias)
    else:
        print("\n[MODO AUTOMÁTICO] Ejecutando todas las instancias...\n")

        for grupo, carpeta in carpetas.items():
            print(f"> Carpeta {grupo} ({carpeta})")

            archivos = glob.glob(os.path.join(carpeta, "*.atsp"))

            for archivo in archivos:
                print(f"   - {os.path.basename(archivo)} ... ")
                n, dist = leer_archivo_tsplib(archivo)
                res = solve_atsp_gavish_graves(archivo, n, dist)
                res["Grupo"] = grupo
                resultados.append(res)
                print("   ✓ Terminado")

    # exportar resultados   
    df = pd.DataFrame(resultados)
    df = df[[
        "Grupo", "Instancia", "Nodos",
        "Vars", "Restr",
        "Tiempo (s)", "Gap (%)",
        "Best Bound", "Objetivo"
    ]]

    df.to_csv("Resultados_GG_ATSP.csv", sep=';', index=False)

    print("\nCSV generado: Resultados_GG_ATSP.csv")
    print("¡Proceso completado!")
