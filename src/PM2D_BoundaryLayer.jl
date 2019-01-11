function CalcStagnationPoint(panels::Array{Panel2D},
                             strengths::Array{Float64},
                             oper_cond::Array{Float64})

    # Unpacking
    u_inf = unpack_oper_cond(oper_cond)
    NPAN = length(panels)

    # Determine velocity profile
    u_panels = VelocityProfile(panels,strengths,oper_cond)
    u_abs = zeros(NPAN)

    cur_min = 1000.0
    stag_panel = 0

    for n=1:NPAN
        u_abs[n] = norm(abs(norm(u_panels[n])))
        if u_abs[n] < cur_min
            cur_min = u_abs[n]
            stag_panel = n
        end
    end

    return stag_panel
end

function CalcBoundaryLayer(panels::Array{Panel2D},
                           oper_cond::Array{Float64},
                           stag_panel::Int64)

    # Unpacking
    u_inf = unpack_oper_cond(oper_cond)
    NPAN = length(panels)

    


end
