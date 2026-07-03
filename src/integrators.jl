
#explicitos
function integrate_euler(fderiv, s, tinit, tend, dt)

    s_actual = isa(s, Vector) ? collect(Float64, s) : Float64(s)
    
    historial_t = Float64[]
    historial_s = isa(s_actual, Vector) ? Vector{Float64}[] : Float64[]
    
    t = tinit

    while t <= tend
        push!(historial_t, t)
        push!(historial_s, isa(s_actual, Vector) ? copy(s_actual) : s_actual)
        
        dsdt = fderiv(s_actual, t)

        if isa(s_actual, Vector)
            s_actual .+= dt .* dsdt
        else
            s_actual += dt * dsdt
        end
        
        t += dt
    end
    
    return historial_t, historial_s #retorna las listas llenas
end



#implicitos