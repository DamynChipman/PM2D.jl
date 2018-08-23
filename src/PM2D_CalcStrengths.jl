"""
    `CalcStrengths(panels,oper_cond)`

Given a list of Panel2D objects and the operation conditions, calculates the
strengths of the elements. Implements a Hans-Smith panel method. The Hans-Smith
method uses a distribution of constant strength source panels with varying
strengths for each panel and a constant strength vortex distribution across
the entire body. This allows for easy implementation of the Kutta Condition.

# ARGUMENTS
* `panels::Array{Panel2D}`       : Array of Panel2D objects representing a body
* `oper_cond::Array{Float64}`    : Operation conditions, [U_INF,ALPHA]

# OUTPUTS
* `A::Array{Float64}`            : Coefficient matrix, size: (N+1)x(N+1)
* `b::Array{Float64}`            : RHS vector
* `strengths::Array{Float64}`    : Array of element strengths
"""
function CalcStrengths(panels::Array{Panel2D},
                       oper_cond::Array{Float64})

    NPAN = length(panels)
    unit_str = 1.0
    u_inf = unpack_oper_cond(oper_cond)

    A = zeros(NPAN+1,NPAN+1)
    b = zeros(NPAN+1)

    # Tangent Flow Conditions
    for i=1:NPAN
        n_hat_pan = panels[i].n_hat
        for j=1:NPAN
            u_pan = CalcVelocity(panels[j],
                                 panels[i].rC,
                                 unit_str,
                                 "source")
            A[i,j] = dot(u_pan,n_hat_pan)
        end

        u_vor = [0.0, 0.0]
        for j=1:NPAN
            u_vor = u_vor + CalcVelocity(panels[j],
                                         panels[i].rC,
                                         unit_str,
                                         "vortex")
        end
        A[i,NPAN+1] = dot(u_vor,n_hat_pan)
    end

    # Kutta Condition
    t_hat_1 = panels[1].t_hat
    t_hat_N = panels[NPAN].t_hat
    for j=1:NPAN
        u_pan_1 = CalcVelocity(panels[j],
                               panels[1].rC,
                               unit_str,
                               "source")
        u_pan_N = CalcVelocity(panels[j],
                               panels[NPAN].rC,
                               unit_str,
                               "source")
        A[NPAN+1,j] = dot(u_pan_1,t_hat_1) + dot(u_pan_N,t_hat_N)
    end

    u_vor_1 = [0.0, 0.0]
    u_vor_N = [0.0, 0.0]
    for j=1:NPAN
        u_vor_1 = u_vor_1 + CalcVelocity(panels[j],
                                         panels[1].rC,
                                         unit_str,
                                         "vortex")
        u_vor_N = u_vor_N + CalcVelocity(panels[j],
                                         panels[NPAN].rC,
                                         unit_str,
                                         "vortex")
    end
    A[NPAN+1,NPAN+1] = dot(u_vor_1,t_hat_1) + dot(u_vor_N,t_hat_N)

    # RHS Vector
    for i=1:NPAN
        b[i] = -dot(u_inf,panels[i].n_hat)
    end

    b[NPAN+1] = -(dot(u_inf,t_hat_1) + dot(u_inf,t_hat_N))


    # Strengths Vector
    strengths = A\b

    return A,b,strengths
end
