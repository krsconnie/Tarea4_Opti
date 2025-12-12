import pandas as pd
import os
from gurobipy import *

# No me roben la licencia porfavor :C
options = {
    "WLSACCESSID": "291dcd15-62ac-4c10-9cf8-3195e3506067",
    "WLSSECRET": "bd688515-ce3d-4361-a937-350ca6c3990f",
    "LICENSEID": 2734085,
}
# Modo puede ser acotado o no_acotado
MODO = "acotado" # en el paper dicen que es mejor no acotarlo para el solver, pero el problema general lo formula así
 
CARPETA_INSTANCIAS = "/home/coni/Tarea4_Opti/MTZ_GUROBI/instancias" # sSte formato me funciona más que poner "instancias" sola (no lo encuentra), no sé porque
ARCHIVO_SALIDA = f'resultados_mtz_{MODO}.csv'

MIS_INSTANCIAS = [
    'br17.atsp', 'ftv33.atsp', 'ftv55.atsp', 'ftv64.atsp',  # 4 instancias pequeñas
    'ftv70.atsp', 'kro124p.atsp', 'ftv170.atsp', # 3 medianas
    'rbg323.atsp', 'rbg358.atsp', 'rbg403.atsp' # 3 grandes, aunque 403 se sale un poco del rango específicado pero no por tanto
]

def leer_instancia_atsp(filepath):
    if not os.path.exists(filepath):
        print(f"Error: archivo no encontrado {filepath}")
        return None, None

    with open(filepath, 'r') as f:
        tokens = f.read().split()

    n = 0
    start_index = -1

    for i, token in enumerate(tokens):
        if token in ["DIMENSION:"]:
            if i + 1 < len(tokens):
                if tokens[i+1] == ":":
                    n = int(tokens[i+2])
                else:
                    try:
                        n = int(tokens[i+1])
                    except ValueError:
                        pass
        
        if token == "EDGE_WEIGHT_SECTION":
            start_index = i + 1
            break
            
    if n == 0 or start_index == -1:
        print(f"Error: Parseo fallido {filepath}")
        return None, None

    valores_matriz = []
    for k in range(start_index, len(tokens)):
        token = tokens[k]
        if token == "EOF":
            break
        if token.lstrip('-').isdigit(): 
            valores_matriz.append(int(token))

    if len(valores_matriz) != n * n:
        print(f"Error: Tamaño no corresponde {filepath}")
        if len(valores_matriz) < n * n:
            return None, None

    c = []
    for i in range(n):
        fila = valores_matriz[i*n : (i+1)*n]
        c.append(fila)

    return n, c

def resolver_instancia_mtz(nombre_archivo, n, c, modo, env):
    I = [i for i in range(n)]        
    I_u = [i for i in range(1, n)]  

    try:
        mdl = Model(f'ATSP_MTZ_{nombre_archivo}', env=env)
    except GurobiError as e:
        print(f"Error creando modelo: {e}")
        return None

    mdl.setParam('TimeLimit', 3600)
    mdl.setParam('OutputFlag', 1)

    x = mdl.addVars(I, I, vtype=GRB.BINARY, name='x')
    u = mdl.addVars(I_u, vtype=GRB.CONTINUOUS, lb=0, name='u')

    mdl.setObjective(
        quicksum(c[i][j] * x[i,j] for i in I for j in I if i != j),
        GRB.MINIMIZE
    )

    for i in I:
        mdl.addConstr(quicksum(x[i,j] for j in I if i!=j) == 1)

    for j in I:
        mdl.addConstr(quicksum(x[i,j] for i in I if i!=j) == 1)

    for i in I_u:
        for j in I_u:
            if i != j:
                mdl.addConstr(u[i] - u[j] + (n - 1) * x[i,j] <= n - 2)

    if modo == "acotado":
        for i in I_u:
            mdl.addConstr(u[i] >= 1)
            mdl.addConstr(u[i] <= n - 1)

    mdl.optimize()

    if mdl.SolCount > 0:
        gap = mdl.MIPGap * 100
        obj = mdl.ObjVal 
    else:
        gap = 100.0
        obj = float('inf') 

    res = {
        "Instancia": nombre_archivo,
        "Nodos": n,
        "Variables": mdl.NumVars,
        "Restricciones": mdl.NumConstrs,
        "Tiempo_s": round(mdl.Runtime, 2),
        "Gap_Porcentaje": round(gap, 2),
        "Funcion_Objetivo": round(obj, 2)
    }

    mdl.dispose()
    return res

if __name__ == "__main__":
    resultados_lista = []
    
    if not os.path.exists(CARPETA_INSTANCIAS):
        print(f"Error: Directorio '{CARPETA_INSTANCIAS}' no encontrado.")
        exit()

    try:
        env = Env(params=options)
        env.start()
    except GurobiError as e:
        print(f"License Error: {e}")
        exit()

    for archivo in MIS_INSTANCIAS:
        ruta = os.path.join(CARPETA_INSTANCIAS, archivo)
        n, matriz_c = leer_instancia_atsp(ruta)
        
        if n and matriz_c:
            resultado = resolver_instancia_mtz(archivo, n, matriz_c, MODO, env)
            if resultado:
                resultados_lista.append(resultado)
                
                df = pd.DataFrame(resultados_lista)
                columnas_ordenadas = [
                    "Instancia", "Nodos", "Variables", "Restricciones", 
                    "Tiempo_s", "Gap_Porcentaje", "Funcion_Objetivo"
                ]
                df = df[columnas_ordenadas]
                df.to_csv(ARCHIVO_SALIDA, index=False)
        else:
            print(f"Skipping {archivo}")

    env.dispose()
    if resultados_lista:
        print(pd.DataFrame(resultados_lista)[["Instancia", "Nodos", "Gap_Porcentaje", "Funcion_Objetivo"]].to_string())