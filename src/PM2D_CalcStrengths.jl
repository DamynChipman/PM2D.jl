#===============================================================================
# CalcStrengths Function
#
#   Given a list of Panel2D objects and the operation conditions (i.e. u_inf),
#   calculates the strengths of the elements. Implements the Kutta Condition
#   using a distribution of source and vortex elements. Solves using A\b
#   notation.
#
#   CalcStrengths(f_vel,panels,u_inf,doKutta) -> A,b,strengths
#       Inputs:
#       panels::Array{Panel2D} - list of Panel2D objects corresponding to the
#                                locations of panels
#       oper_cond::Array{Float64} - Operation conditions = {U_INF, AoA}
#       check_BC::Bool - Optional argument to print the values of the boundary
#                        conditions (tanget flow and Kutta condition)
#       method::String - Method for element distribution. Options are: "panels",
#                        "HansSmith". Defaults to "panels"
#
#       Outputs:
#       A - coefficient matrix
#       b - RHS vector
#       strengths - list of strengths corresponding to the strengths of the
#                   elements
# =============================================================================#
function CalcStrengths(panels::Array{Panel2D},
                       oper_cond::Array{Float64})

    NPAN = length(panels)
    unit_str = 1.0
    u_inf = unpack_oper_cond(oper_cond)

    A = zeros(NPAN+1,NPAN+1)
    b = zeros(NPAN+1)

    unit_str = 1.0

    # Tangent Flow Conditions
    for i=1:NPAN
        n_hat_pan = panels[i].n_hat
        for j=1:NPAN
            u_pan = CalcVelocity(panels[i],
                                 panels[j].rC,
                                 unit_str,
                                 "source")
            A[j,i] = dot(u_pan,n_hat_pan)
        end

        u_vor = [0.0, 0.0]
        for j=1:NPAN
            u_vor = u_vor + CalcVelocity(panels[i],
                                         panels[j].rC,
                                         unit_str,
                                         "vortex")
        end
        A[NPAN+1,i] = dot(u_vor,n_hat_pan)
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
        A[j,NPAN+1] = dot(u_pan_1,t_hat_1) + dot(u_pan_N,t_hat_N)
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
        b[i] = dot(u_inf,panels[i].n_hat)
    end

    b[NPAN+1] = dot(u_inf,t_hat_1) + dot(u_inf,t_hat_N)

    # --- Strengths Vector ---
    strengths = transpose(A)\b

    return A,b,strengths
end
