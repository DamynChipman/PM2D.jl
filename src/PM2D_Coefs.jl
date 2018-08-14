#===============================================================================
# CalcCp Function
#
#
#
#
# =============================================================================#
function CalcCp(panels::Array{Panel2D},
                 strengths::Array{Float64},
                 oper_cond::Array{Float64})

    u_inf = unpack_oper_cond(oper_cond)
    U_INF = oper_cond[1]
    NPAN = length(panels)
    Cp = Array{Float64}(NPAN)
    rC = Array{Float64}(NPAN,2)

    for n=1:NPAN
        u_ind = InducedVelocity(panels,
                                strengths,
                                oper_cond,
                                panels[n].rC)
        u_tan = dot(u_ind,panels[n].t_hat)
        Cp[n] = 1 - (u_tan/U_INF)^2
        rC[n,:] = panels[n].rC
    end

    return Cp,rC
end

#===============================================================================
# CalcCl Function
#
#
#
#
# =============================================================================#
function CalcCl(panels::Array{Panel2D},
                 strengths::Array{Float64},
                 oper_cond::Array{Float64},
                 method::String="panels")

    N_panels = length(panels)
    gamma = strengths[N_panels+1]
    U_INF = oper_cond[1]

    sum = 0.0
    for j=1:N_panels
        sum = sum + panels[j].L
    end
    Cl = ((2*gamma)/U_INF)*sum
    return Cl
end
