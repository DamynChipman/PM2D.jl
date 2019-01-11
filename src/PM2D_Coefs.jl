"""
    `CalcCp(panels,strengths,oper_cond)`

Computes the coefficient of pressure at each collocation point of panels.

# ARGUMENTS
* `panels::Array{Panel2D}`        : Array of Panel2D objects representing a body
* `strengths::Array{Float64}`     : Array of element strengths
* `oper_cond::Array{Float64}`     : Operation conditions, [U_INF,ALPHA]

# OUTPUTS
* `Cp::Array{Float64}`            : Array of pressure coefficients at each collocation point
* `rC::Array{Float64}`            : Array of points corresponding to where Cp was calculated
"""
function CalcCp(panels::Array{Panel2D},
                strengths::Array{Float64},
                oper_cond::Array{Float64})

    u_inf = unpack_oper_cond(oper_cond)
    U_INF = oper_cond[1]
    NPAN = length(panels)
    Cp = Array{Float64}(NPAN)
    rC = Array{Float64}(NPAN,2)

    for n=1:NPAN
        rC_off = panels[n].rC_offset
        factor = 0.0
        eval_pt = panels[n].rC-(factor*rC_off).*panels[n].n_hat
        u_ind = InducedVelocity(panels,
                                strengths,
                                eval_pt)
        u_tan = dot(u_ind + u_inf,panels[n].t_hat)
        Cp[n] = 1 - (u_tan/U_INF)^2
        rC[n,:] = panels[n].rC
    end

    return Cp,rC
end

"""
    `CalcCl(panels,strengths,oper_cond)`

Computes the coefficient of lift over the entire body.

# ARGUMENTS
* `panels::Array{Panel2D}`        : Array of Panel2D objects representing a body
* `strengths::Array{Float64}`     : Array of element strengths
* `oper_cond::Array{Float64}`     : Operation conditions, [U_INF,ALPHA]

# OUTPUTS
* `Cl::Float64`                   : Coefficient of lift
"""
function CalcCl(panels::Array{Panel2D},
                strengths::Array{Float64},
                oper_cond::Array{Float64})

    NPAN = length(panels)
    gamma = strengths[NPAN+1]
    U_INF = oper_cond[1]

    sum = 0.0
    for j=1:NPAN
        sum = sum + panels[j].L
    end
    Cl = ((2*gamma)/U_INF)*sum
    return Cl
end
