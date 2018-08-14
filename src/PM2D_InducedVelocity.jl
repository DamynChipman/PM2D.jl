#===============================================================================
# InducedVelocity Function
#
#
#
#
#
# =============================================================================#
function InducedVelocity(panels::Array{Panel2D},
                         strengths::Array{Float64},
                         oper_cond::Array{Float64},
                         R::Array{Float64})

    # Unpacking
    NPAN = length(panels)
    u_inf = unpack_oper_cond(oper_cond)

    u_ind = u_inf

    for j=1:NPAN
        u_ind = u_ind + CalcVelocity(panels[j],
                                     R,
                                     strengths[j],
                                     "source") +
                         CalcVelocity(panels[j],
                                      R,
                                      strengths[NPAN+1],
                                      "vortex")
    end
    return u_ind
end
