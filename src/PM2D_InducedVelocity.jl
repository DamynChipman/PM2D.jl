"""
    `InducedVelocity(panels,strengths,R)`

Given a body of panels, computes the induced velocity from all of the panels.
Uses a Hess-Smith style of element distribution.

# ARGUMENTS
* `panels::Array{Panel2D}`          : Array of Panel2D objects representing a body
* `strengths::Array{Float64}`       : Array of element strengths
* `R::Array{Float64}`               : Evaluation point

# OUTPUTS
* `u_ind::Array{Float64}`           : Vector of induced velocity
"""
function InducedVelocity(panels::Array{Panel2D},
                         strengths::Array{Float64},
                         R::Array{Float64},
                         method::String)

    # Unpacking
    NPAN = length(panels)

    u_ind = [0.0,0.0]

    if method == "Hess-Smith"
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
    elseif method == "Vortex Sheet"
        for j=1:NPAN
            u_ind = u_ind + CalcVelocity(panels[j],
                                         R,
                                         strengths[j],
                                         "vortex")
        end

    else
        error("Invalid 'method' in InducedVelocity.")
    end

    return u_ind
end

"""
    `VelocityProfile(panels,strengths,oper_con)`

Calculates the velocity at each collocation of the body of panels

# ARGUMENTS
* `panels::Array{Panel2D}`        : Array of Panel2D objects representing a body
* `strengths::Array{Float64}`     : Array of element strengths
* `oper_cond::Array{Float64}`     : Operation conditions, [U_INF,ALPHA]

# OUTPUTS
* `u_panels::Array{Float64}`      : Array of velocity vectors at each collocation point
"""
function VelocityProfile(panels::Array{Panel2D},
                         strengths::Array{Float64},
                         oper_cond::Array{Float64},
                         method::String)

    # Unpacking
    u_inf = unpack_oper_cond(oper_cond)
    NPAN = length(panels)

    u_panels = zeros(NPAN,2)

    if method == "Hess-Smith"
        for j=1:NPAN
            u_panels[j,:] = u_inf + InducedVelocity(panels,
                                                    strengths,
                                                    panels[j].rC)
        end
    elseif method == "Vortex Sheet"
        for j=1:NPAN
            u_panels[j,:] = u_inf + InducedVelocity(panels,
                                                    strengths,
                                                    panels[j].rC)
        end
    else
        error("Invalid 'method' in InducedVelocity.")
    end

    return u_panels
end
