# Simulación numérica de reacciones químicas con Julia

Este repositorio contiene el código y la documentación del proyecto
*Estudio Comparativo de Estabilidad y Rigidez de Integradores Numéricos en
Modelos de Cinética Química*, enfocado en la implementación y comparación
de métodos numéricos como Euler explícito, Euler implícito y regla trapezoidal
aplicados a sistemas de ecuaciones diferenciales ordinarias que modelan
reacciones de química cinética con significativa rigidez

## Estructura
- `code/julia/`: Código fuente en Julia con los integradores temporales y scripts principales.
- `notebooks/`: Notebooks interactivos de **Pluto.jl** con explicaciones paso a paso.
- `results/`: Figuras y tablas generadas.
- `docs/`: Informe final y presentación.

## Requisitos
- Julia 1.9 o superior.
- Paquetes: `Plots`, `Pluto`, `LinearAlgebra`, `LaTeXStrings`.

## Instalación rápida
1. Clona el repositorio:
   ```bash
   git clone https://github.com/spartan14799/simulacion-reacciones-quimicas.git
