### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ a2fe7350-7672-11f1-ab48-dbc0621b8d90
begin
    using Pkg, LinearAlgebra, ForwardDiff, Plots
    Pkg.activate("..") #activa el ambiente en root

    include("../src/integrators.jl")
end

# ╔═╡ 44f3e388-a052-43d8-88c6-85d8ec17c5af
md"""
Este notebook tiene como objetivo verificar la correcta implementación de los integradores numéricos explícitos e implícitos definidos en `src/integrators.jl` mediante la resolución del problema de valor inicial

```math
\begin{cases}
y'(t)=-5y(t),\\
y(0)=1.
\end{cases}
```

Que es un ejemplo canónico de una ecuación diferencial rígida y tiene solución exacta
```math
y(t)=e^{-5t},
```

Lo que permite validar comparativamente la implementación de los integradores.

"""

# ╔═╡ 4c39df53-8975-4e6e-9649-afe9946c290d
#definir la edo
f(s, t) = [-5*s[1]]

# ╔═╡ 9a957e08-e74f-49d5-94c5-beadce2b8352
#parametros
begin
    s0 = [1.0]

    tinit = 0.0
    tend = 2.0
    dt = 0.45
end

# ╔═╡ bcc935b4-d812-4a21-93ce-4cf08095acb5
md"""
## Métodos Explícitos

- Euler explícito
- Heun (Runge-Kutta 2)
- Runge-Kutta 4
"""

# ╔═╡ 5fcfe3bb-8657-424f-9269-f9630569603a
begin
	t_euler, s_euler = integrate_euler(f, s0, tinit, tend, dt)
	t_heun, s_heun = integrate_heun(f, s0, tinit, tend, dt)
	t_rk4, s_rk4 = integrate_rk4(f, s0, tinit, tend, dt)
end

# ╔═╡ db817273-47cc-4d87-9c22-c5ff8831efa8
md"""
## Métodos Implícitos
- Euler implícito
- Crank Nicolson
"""

# ╔═╡ c5759cd2-1127-41e9-bec2-2db55362f343
begin
	t_ie, s_ie = integrate_implicit_euler(f, s0, tinit, tend, dt)
	t_cn, s_cn = integrate_crank_nicolson(f, s0, tinit, tend, dt)
end

# ╔═╡ 1ff40f5b-b58f-4b9c-bbc0-dd396bc62669
md"""
- Solución exacta
"""

# ╔═╡ 7cdfbfe3-51e5-4188-92c8-cba8d2e6f070

begin
	t_exacta = range(tinit, tend, length=1000) #malla mas fina pa la exacta
	sol_exacta = exp.(-5 .* t_exacta)
end

# ╔═╡ 7d97bb6f-484a-4565-9f6d-34f3de6ec01d
md"""
## Gráficas
"""

# ╔═╡ 9a4e667d-5f7f-4cf6-86cb-f57f215481b1
begin
	plot(
	    t_euler,
	    first.(s_euler),
	    label="euler explícito",
	    lw=2,
	    xlabel="t",
	    ylabel="y(t)"
	)
	
	plot!(t_heun, first.(s_heun), label="heun", lw=2)
	plot!(t_rk4, first.(s_rk4), label="rk4", lw=2)
	plot!(t_ie, first.(s_ie), label="euler implícito", lw=2)
	plot!(t_cn, first.(s_cn), label="crank nicolson", lw=2)
	
	plot!(
	    t_exacta,
	    sol_exacta,
	    label="Exacta",
	    lw=3,
	    ls=:dash
	)
end

# ╔═╡ 1b49054a-915a-4966-b7b3-2f5e6f1bbba9
md"""
# Comparación de Integradores

A continuación se comparan las aproximaciones obtenidas con los métodos explícitos e implícitos para el PVI propuesto, junto con la solución exacta.

Al tratarse de un problema rígido, los métodos explícitos tienen problemas de estabilidad cuando el tamaño de paso es relativamente grande.

- **Euler explícito** oscila y diverge.

- **Heun** mejora el comportamiento de Euler explícito, pero no reproduce adecuadamente la solución.

- **RK4**, a pesar de ser un método de cuarto orden, resulta menos eficaz que Euler implícito, un método implícito de primer orden.

Por el contrario, los métodos implícitos son estables y aproximan cercanamente la solución exacta con el mismo tamaño de paso.

"""

# ╔═╡ a0bf1597-f1f6-4be5-892d-8b9d979f7ba1


# ╔═╡ Cell order:
# ╟─a2fe7350-7672-11f1-ab48-dbc0621b8d90
# ╟─44f3e388-a052-43d8-88c6-85d8ec17c5af
# ╠═4c39df53-8975-4e6e-9649-afe9946c290d
# ╠═9a957e08-e74f-49d5-94c5-beadce2b8352
# ╟─bcc935b4-d812-4a21-93ce-4cf08095acb5
# ╠═5fcfe3bb-8657-424f-9269-f9630569603a
# ╟─db817273-47cc-4d87-9c22-c5ff8831efa8
# ╠═c5759cd2-1127-41e9-bec2-2db55362f343
# ╟─1ff40f5b-b58f-4b9c-bbc0-dd396bc62669
# ╠═7cdfbfe3-51e5-4188-92c8-cba8d2e6f070
# ╟─7d97bb6f-484a-4565-9f6d-34f3de6ec01d
# ╠═9a4e667d-5f7f-4cf6-86cb-f57f215481b1
# ╟─1b49054a-915a-4966-b7b3-2f5e6f1bbba9
# ╠═a0bf1597-f1f6-4be5-892d-8b9d979f7ba1
