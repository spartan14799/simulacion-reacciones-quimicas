### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 255fda83-9ecd-4998-8194-0606af525e91
begin
    using Pkg, LinearAlgebra, ForwardDiff, Plots
    Pkg.activate("..") #activa el ambiente en root

    include("../src/integrators.jl")
end

# ╔═╡ 242a6a10-7704-11f1-a2dd-51d642454155
md"""
Este notebook compara los integradores numéricos implementados en `src/integrators.jl` mediante la resolución de un problema de valor inicial proveniente de un modelo de química atmosférica.
"""

# ╔═╡ 7537d8cc-90a7-4bf4-976e-501ce0b484aa
md"""
# Modelo Atmosférico

Se considera el problema de valor inicial

```math
\mathbf{u}'(t)=\mathbf{f}(t,\mathbf{u}(t)), \qquad \mathbf{u}(t_0)=\mathbf{u}_0,
```

donde

```math
\mathbf{u}(t)=
\begin{pmatrix}
u_1\\u_2\\u_3\\u_4
\end{pmatrix}, \qquad
\mathbf{u}_0=
\begin{pmatrix}
0\\
1.3\times10^8\\
5.0\times10^{11}\\
8.0\times10^{11}
\end{pmatrix},
```

y el campo vectorial está dado por

```math
\mathbf{f}(t,\mathbf{u})=
\begin{pmatrix}
k_1(t)u_3-k_2u_1\\
k_1(t)u_3-k_3u_2u_4+\sigma_2\\
k_3u_2u_4-k_1(t)u_3\\
k_2u_1-k_3u_2u_4
\end{pmatrix}.
```

Las constantes cinéticas se definen como

```math
k_1(t)=
\begin{cases}
10^{-5}e^{7\,\mathrm{solar}(t)}, & \text{día},\\
10^{-40}, & \text{noche},
\end{cases}
```

con intervalo de integración

```math
t_0=14400,\qquad t_f=504000.
```
"""

# ╔═╡ 9e06c8c0-f3a2-4514-b106-ba20821ca323

#params
begin
	k2 = 1e5
	k3 = 1e-16
	sgm2 = 1e6
end

# ╔═╡ 8e4b1d4b-6a2d-43db-b869-85c4dd35784a

#necesaria para k1
function solar(t)
    th = t/3600
    th_prom = th - 24*floor(th/24)

    if 4 <= th_prom <= 20
        return max(sin((π/16)*(th_prom - 4)), 0.0)^0.2
    else
        return 0.0
    end
end

# ╔═╡ 4a47a8b8-ff54-4f56-98a3-7354ed3bda83
function k1(t)
    s = solar(t)
    return s > 0 ? 1e-5 * exp(7*s) : 1e-40
end

# ╔═╡ 4afaa9d8-5a52-458e-a7b5-c05c625b5870
#pvi
begin
	u0 = [
	    0.0,
	    1.3e8,
	    5.0e11,
	    8.0e11]
	
	t0 = 14400.0 #4 am dia 1
	tf = 504000.0 # 8 pm día 6
	h = 600.0
end

# ╔═╡ 7ad69509-2f43-43c2-aaa8-fd926988c0d6

#funcion del sistema

f(u,t) = [
    k1(t)*u[3] - k2*u[1],
    k1(t)*u[3] - k3*u[2]*u[4] + sgm2,
    k3*u[2]*u[4] - k1(t)*u[3],
    k2*u[1] - k3*u[2]*u[4]
]

# ╔═╡ e8c537f7-970a-45f9-a4b4-e4166f35cfd2

begin
	t_ee, s_ee = integrate_euler(f, u0, t0, tf, h)
	t_h, s_h = integrate_heun(f, u0, t0, tf, h)
	t_rk4, s_rk4 = integrate_rk4(f, u0, t0, tf, h)
	
	t_ie,s_ie = integrate_implicit_euler(f,u0,t0,tf,h)
	t_cn,s_cn = integrate_crank_nicolson(f,u0,t0,tf,h)
end

# ╔═╡ 52121419-5ecf-4d6d-9e6d-2f7a987076cd
md"""
## Concentración $u_1$
"""

# ╔═╡ 764c832c-ed02-473e-a465-6dde8cff14f6

begin
	plot(t_ie, first.(s_ie), label="implicit euler", lw=2)
	plot!(t_cn, first.(s_cn), label="crank nicolson")
end

# ╔═╡ b4b6fa1a-56d8-4c9a-b738-e38aba2fb754
begin
	plot(t_ee, first.(s_ee), label="explicit euler", lw=2)
	plot!(t_h, first.(s_h), label="heun")
	plot!(t_rk4, first.(s_rk4), label = "RK4")
end

# ╔═╡ 730ab82b-62bc-4398-add0-cd5b0591c4a7
md"""
## Concentración $u_2$
"""

# ╔═╡ 68ced008-fb0e-4237-9236-c8de367b3d4b

begin
	plot(t_ie, getindex.(s_ie,2), label="implicit euler")
	plot!(t_cn, getindex.(s_cn,2), label="crank nicolson")
end

# ╔═╡ af083967-3786-438d-8ead-9c5e1b42e197
begin
	plot(t_ee, getindex.(s_ee,2), label="explicit euler", lw=2)
	plot!(t_h, getindex.(s_h,2), label="heun")
	plot!(t_rk4, getindex.(s_rk4,2), label = "RK4")
end

# ╔═╡ 95d92b1a-230e-4078-ab3a-8002a8b96c0f
md"""
## Concentración $u_3$
"""

# ╔═╡ ed3de115-48d7-46f4-8b82-4889bb222635
begin
	plot(t_ie, getindex.(s_ie,3), label="implicit euler")
	plot!(t_cn, getindex.(s_cn,3), label="crank nicolson")
end

# ╔═╡ 9efc77b3-3aa8-4bde-b00b-0b15f63108b4
begin
	plot(t_ee, getindex.(s_ee,3), label="explicit euler", lw=2)
	plot!(t_h, getindex.(s_h,3), label="heun")
	plot!(t_rk4, getindex.(s_rk4,3), label = "RK4")
end

# ╔═╡ bb50e7a4-e84a-4282-b36e-b1f13190c942

md"""
## Concentración $u_4$
"""

# ╔═╡ 7e557ce5-bad8-4c3c-8476-054e74ee96b0
begin
	plot(t_ie, getindex.(s_ie,4), label="implicit euler")
	plot!(t_cn, getindex.(s_cn,4), label="crank nicolson")
end

# ╔═╡ 796d9c8a-285a-4b41-bd93-33c336393cc4
begin
	plot(t_ee, getindex.(s_ee,3), label="explicit euler", lw=2)
	plot!(t_h, getindex.(s_h,3), label="heun")	
	plot!(t_rk4, getindex.(s_rk4,3), label = "RK4")
end

# ╔═╡ c99f18db-f9e4-410e-b620-d9967a6e5216

md"""
## Conclusiones

Los resultados obtenidos confirman el carácter rígido del sistema considerado, proveniente de un modelo de química atmosférica con escalas de tiempo ampliamente separadas.

Los métodos explícitos (Euler, Heun y RK4) presentan inestabilidad numérica para tamaños de paso moderados, manifestándose en divergencia de las soluciones y crecimiento no físico de las variables.

Por el contrario, los métodos implícitos (Euler implícito y Crank–Nicolson) muestran un comportamiento estable y reproducen adecuadamente la dinámica cualitativa del sistema incluso para pasos relativamente grandes.

Estos resultados evidencian que, en problemas rígidos, la estabilidad numérica del método es un criterio más determinante que su orden de precisión.
"""

# ╔═╡ Cell order:
# ╟─255fda83-9ecd-4998-8194-0606af525e91
# ╟─242a6a10-7704-11f1-a2dd-51d642454155
# ╟─7537d8cc-90a7-4bf4-976e-501ce0b484aa
# ╠═9e06c8c0-f3a2-4514-b106-ba20821ca323
# ╟─8e4b1d4b-6a2d-43db-b869-85c4dd35784a
# ╟─4a47a8b8-ff54-4f56-98a3-7354ed3bda83
# ╠═4afaa9d8-5a52-458e-a7b5-c05c625b5870
# ╠═7ad69509-2f43-43c2-aaa8-fd926988c0d6
# ╠═e8c537f7-970a-45f9-a4b4-e4166f35cfd2
# ╟─52121419-5ecf-4d6d-9e6d-2f7a987076cd
# ╠═764c832c-ed02-473e-a465-6dde8cff14f6
# ╠═b4b6fa1a-56d8-4c9a-b738-e38aba2fb754
# ╟─730ab82b-62bc-4398-add0-cd5b0591c4a7
# ╠═68ced008-fb0e-4237-9236-c8de367b3d4b
# ╠═af083967-3786-438d-8ead-9c5e1b42e197
# ╟─95d92b1a-230e-4078-ab3a-8002a8b96c0f
# ╠═ed3de115-48d7-46f4-8b82-4889bb222635
# ╠═9efc77b3-3aa8-4bde-b00b-0b15f63108b4
# ╟─bb50e7a4-e84a-4282-b36e-b1f13190c942
# ╠═7e557ce5-bad8-4c3c-8476-054e74ee96b0
# ╠═796d9c8a-285a-4b41-bd93-33c336393cc4
# ╠═c99f18db-f9e4-410e-b620-d9967a6e5216
