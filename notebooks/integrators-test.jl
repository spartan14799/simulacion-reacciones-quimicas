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
## Explicit Methods
"""

# ╔═╡ 5fcfe3bb-8657-424f-9269-f9630569603a
t_euler, s_euler = integrate_euler(f, s0, tinit, tend, dt)

# ╔═╡ f25f047e-fc11-4609-9097-9ecdc021bb89
t_heun, s_heun = integrate_heun(f, s0, tinit, tend, dt)

# ╔═╡ ded78a5d-ee47-4e2c-9b0a-e569c977b51c
t_rk4, s_rk4 = integrate_rk4(f, s0, tinit, tend, dt)

# ╔═╡ db817273-47cc-4d87-9c22-c5ff8831efa8
md"""
## Implicit Methods
"""

# ╔═╡ c5759cd2-1127-41e9-bec2-2db55362f343
t_ie, s_ie = integrate_implicit_euler(f, s0, tinit, tend, dt)

# ╔═╡ 4f60466c-e41c-46eb-b4be-afb14209b4b1
t_cn, s_cn = integrate_crank_nicolson(f, s0, tinit, tend, dt)

# ╔═╡ 1ff40f5b-b58f-4b9c-bbc0-dd396bc62669
md"""
## Exact Solution
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

# ╔═╡ Cell order:
# ╠═a2fe7350-7672-11f1-ab48-dbc0621b8d90
# ╠═4c39df53-8975-4e6e-9649-afe9946c290d
# ╠═9a957e08-e74f-49d5-94c5-beadce2b8352
# ╟─bcc935b4-d812-4a21-93ce-4cf08095acb5
# ╠═5fcfe3bb-8657-424f-9269-f9630569603a
# ╠═f25f047e-fc11-4609-9097-9ecdc021bb89
# ╠═ded78a5d-ee47-4e2c-9b0a-e569c977b51c
# ╟─db817273-47cc-4d87-9c22-c5ff8831efa8
# ╠═c5759cd2-1127-41e9-bec2-2db55362f343
# ╠═4f60466c-e41c-46eb-b4be-afb14209b4b1
# ╟─1ff40f5b-b58f-4b9c-bbc0-dd396bc62669
# ╠═7cdfbfe3-51e5-4188-92c8-cba8d2e6f070
# ╟─7d97bb6f-484a-4565-9f6d-34f3de6ec01d
# ╠═9a4e667d-5f7f-4cf6-86cb-f57f215481b1
