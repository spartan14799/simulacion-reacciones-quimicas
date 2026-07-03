
#explicitos
function integrate_euler(fderiv, s, tinit, tend, dt)

    s_actual = copy(s) #copia para no modificar el vector original

    historial_t = Float64[] #vector vacio
    historial_s = Vector{Float64}[]

    t = tinit
    while t <= tend

        push!(historial_t, t)
        push!(historial_s, copy(s_actual))

        dsdt = fderiv(s_actual, t)

        s_actual .+= dt .* dsdt

        t += dt
    end

    return historial_t, historial_s
end

function integrate_heun(fderiv, s, tinit, tend, dt)

    s_actual = collect(Float64, s)

    historial_t = Float64[]
    historial_s = Vector{Float64}[]

    t = tinit

    while t <= tend

        push!(historial_t, t)
        push!(historial_s, copy(s_actual))

        k1 = fderiv(s_actual, t)
        k2 = fderiv(s_actual .+ dt .* k1, t + dt)

        s_actual .+= 0.5 .* dt .* (k1 .+ k2)

        t += dt
    end

    return historial_t, historial_s

end


function integrate_rk4(fderiv, s, tinit, tend, dt)

    s_actual = copy(s)

    historial_t = Float64[]
    historial_s = Vector{Float64}[]

    t = tinit

    while t <= tend

        push!(historial_t, t)
        push!(historial_s, copy(s_actual))

        k1 = fderiv(s_actual, t)
        k2 = fderiv(s_actual .+ (dt/2) .* k1, t + dt/2)
        k3 = fderiv(s_actual .+ (dt/2) .* k2, t + dt/2)
        k4 = fderiv(s_actual .+ dt .* k3, t + dt)

        s_actual .+= (dt/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)

        t += dt
    end

    return historial_t, historial_s
end


#implicitos


function integrate_implicit_euler(fderiv, s, tinit, tend, dt)

    s_actual = copy(s)

    historial_t = Float64[]
    historial_s = Vector{Float64}[]

    max_iter = 50
    tolerancia = 1e-9

    t = tinit

    while t <= tend
        
        push!(historial_t, t)
        push!(historial_s, copy(s_actual))

        dsdt = fderiv(s_actual, t)

        #predictor de euler explicito
        s_guess = s_actual .+ dt .* dsdt
        s_next = copy(s_guess)
        #iteracion de punto fijo
        for _ in 1:max_iter

            f_next = fderiv(s_guess, t + dt)
            s_next = s_actual .+ dt .* f_next

            error = norm(s_next - s_guess)

            if error < tolerancia
                break
            end

            s_guess = s_next
        end

        s_actual = s_next
        t += dt
    end

    return historial_t, historial_s
end


function integrate_crank_nicolson(fderiv, s, tinit, tend, dt)

    s_actual = copy(s)

    historial_t = Float64[]
    historial_s = Vector{Float64}[]

    max_iter = 50
    tolerancia = 1e-9

    t = tinit

    while t <= tend

        push!(historial_t, t)
        push!(historial_s, copy(s_actual))

        dsdt = fderiv(s_actual, t)
        s_guess = s_actual .+ dt .* dsdt
        s_next = copy(s_guess)
        for _ in 1:max_iter

            f_next = fderiv(s_guess, t + dt)

            s_next = s_actual .+ 0.5 .* dt .* (dsdt .+ f_next)

            error = norm(s_next - s_guess)

            if error < tolerancia
                break
            end

            s_guess = s_next
        end

        s_actual = s_next
        t += dt
    end

    return historial_t, historial_s
end