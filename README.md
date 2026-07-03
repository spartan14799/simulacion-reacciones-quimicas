
# Comparison of Implicit and Explicit Methods for Stiff Ordinary Differential Equations

## About
This project compares implicit and explicit numerical methods for solving stiff ordinary differential equations derived from environmental chemistry models. The methods are analyzed with respect to stability, accuracy, and computational efficiency through numerical experiments. The results are presented through interactive Pluto notebooks that illustrate the implementation and behavior of each method.

## Structure

- **`src/`**: contains the implementation of the explicit and implicit numerical methods used throughout the project.
- **`notebooks/`**: contains the interactive Pluto notebooks with the numerical experiments, comparisons, and visualizations.

## Pluto Notebooks

To open and run the Pluto notebooks, first start Julia from the `notebooks/` directory.  
Install the required packages (only once):

```julia
using Pkg
Pkg.add(["Pluto", "ForwardDiff", "Plots"])
```

Load the required packages:

```julia
using Pluto, LinearAlgebra, ForwardDiff, Plots
```

Launch Pluto:

```julia
Pluto.run()
```

And once the Pluto interface opens

- Click **Open a notebook**.
- Select the desired `.jl` notebook from the `notebooks/` directory.
- Run the notebook.

## Authors
This project was developed for the course:  
**Cálculo de Ecuaciones Diferenciales Ordinarias**  
**Facultad de Ciencias**  
**Universidad Nacional de Colombia - Sede Bogotá**

Professor:

Martínez Alba, Nicolás - nmartineza@unal.edu.co

Developed by:

- Martínez Sáenz, Sebastián - sebmartinez@unal.edu.co
- Huertas Serrano, Juan Esteban - tjhuertass@unal.edu.co

## License
This project is licensed under the MIT license. See the [LICENCE](LICENSE) file for details.