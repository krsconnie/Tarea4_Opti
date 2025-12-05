# Tarea Computacional 4 
## Objetivos:
  1. Utilizar Python junto a CPLEX y GUROBI para resolver modelos de programación lineal entera.
  2. Conocer distintas formulaciones del problema del vendedor viajero asimetrico (ATSP).
  3. Resolver instancias del ATSP de distintos tamaños.
  4. Comparar los resultados de las distintas formulaciones y solvers.

## Descripción
Sobre la implementación:
    
    Implementar las formulaciones MTZ y GG que aparecen en el artículo: Roberti, Toth. “Models and algorithms 
    for the Asymmetric Traveling Salesman Problem: an experimental comparison” (puede encontrarse en el repositorio).
    En total son 4 programas que se deben resolver las mismas instancias con un tiempo límite de 1 hora:
    
      - MTZ con CPLEX
      - MTZ con GUROBI
      - GG con CPLEX
      - GG con GUROBI
    
Sobre las instancias:

    Las instancias del ATSP son tomadas de TSPLIB2, deben seleccionarse 10 que sean distribuidas en tres tamaños 
    (pequeñas: 17–50 nodos; medianas: 60–120; grandes: 200–400), de tal manera que permitan un correcto análisis 
    de los resultados.

Sobre los resultados a extraer:

    Los resultados a extraer de la ejecución de las pruebas son: número de variables, número de restricciones, 
    tiempo de cómputo (en segundos), gap reportado por el solver (en porcentaje) y función objetivo (Best bound). 
    (Si la instancia no puede ser resuelta por el solver en el tiempo disponible, puede reportar gap 100%.)

## Informe:
Se debe incluir las 2 formulaciones (MTZ y GG), indicando sus ventajas y desventajas respecto a lo que aparece en el paper.
Los resultados deben presentarse en formato tabla, indicando por cada instancia (fila) el número de nodos; luego, por cada modelo y solver, número de variables, número de restricciones, tiempo de cómputo (en segundos), gap reportado por el      solver (MIPgap, en porcentaje) y función objetivo (Best bound).(Asumo, como en esta tabla del paper)

<div align="center">
   <img width="491" height="490" alt="image" src="https://github.com/user-attachments/assets/82c5c74d-6172-4188-96a5-81b32e9f0eae" />
   <p><em>Figura 1: Tabla de referencia del paper de Roberti y Toth</em></p>
</div>

Sobre el análisis de resultados:

    Se debe hacer mención sobre el tamaño de las formulaciones y el GAP  (contestar ¿se confirma lo indicado en el artículo?), 
    los tiempos de cómputo, comportamiento de los solvers, otros que se consideren relevantes.

Sobre las conclusiones:

    Indicar: - Qué se aprendió de las diferentes formulaciones (MTZ y GG)
             - Como se reafirma la teoría indicada en el artículo con los resultados obtenidos.
             - Qué se aprendio respecto a los solvers utilizados.
             - Qué podrían observar respecto a la escalabilidad de las instancias.
             - Qué características tenían las instancias más difíciles.
             - Que pruebas adicionales se podría realizar si se tuviera más tiempo 
               (nuevos modelos, nuevas instancias, análisis estadísticos, etc).


    
