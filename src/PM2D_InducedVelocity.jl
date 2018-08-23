"""
    `InducedVelocity(panels,strengths,oper_cond,R)`

Given a body of panels, computes the induced velocity from all of the panels.
Uses a Hess-Smith style of element distribution.

# ARGUMENTS
* `panels::Array{Panel2D}`          : Array of Panel2D objects representing a body
* `strengths::Array{Float64}`       : Array of element strengths
* `oper_cond::Array{Float64}`       : Operation conditions, [U_INF,ALPHA]
* `R::Array{Float64}`               : Evaluation point

# OUTPUTS
* `u_ind::Array{Float64}`           : Vector of induced velocity
"""
function InducedVelocity(panels::Array{Panel2D},
                         strengths::Array{Float64},
                         oper_cond::Array{Float64},
                         R::Array{Float64})

    # Unpacking
    NPAN = length(panels)
    u_inf = unpack_oper_cond(oper_cond)

    u_ind = [0.0,0.0]

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
