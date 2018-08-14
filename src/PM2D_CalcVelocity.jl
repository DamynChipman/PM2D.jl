#===============================================================================
# CalcVelocity Function
#
#
#
#
#
# =============================================================================#
function CalcVelocity(elem_panel::Panel2D,
                      eval_pt::Array{Float64},
                      str::Float64,
                      panel_type::String)

    # Unpacking variables
    pan = elem_panel
    x,y = eval_pt[1],eval_pt[2]
    A = str/(2*pi)

    velocity = Array{Float64}(2)

    # Unpacking Variables
    x1 = pan.r1[1]
    y1 = pan.r1[2]
    x2 = pan.r2[1]
    y2 = pan.r2[2]

    R1 = [x1,y1]
    R2 = [x2,y2]
    R = [x,y]

    TRANS(A) = [[cos(A) sin(A)];
                [-sin(A) cos(A)]]

    R_prime = R - R1
    R_panel = TRANS(elem_panel.theta)*R_prime

    x = R_panel[1]
    y = R_panel[2]
    L = pan.L

    if panel_type == "source"
        velocity[1] = 0.5*A*log(((x - L)^2 + y^2)/(x^2 + y^2))
        velocity[2] = -A*(atan2(y,x - L) - atan2(y,x))
        velocity = -TRANS(-pan.theta)*velocity
    elseif panel_type == "vortex"
        velocity[2] = -0.5*A*log(((x - L)^2 + y^2)/(x^2 + y^2))
        velocity[1] = -A*(atan2(y,x - L) - atan2(y,x))
        velocity = -TRANS(-pan.theta)*velocity
    elseif panel_type == "doublet"
        x_int_D(x0) = ((x - x0)*y)/(((x - x0)^2 + y^2)^2)
        y_int_D(x0) = ((x - x0)^2 - y^2)/(((x - x0)^2 + y^2)^2)

        velocity[1] = A*quadgk(x_int_D,x1,x2)[1]
        velocity[2] = -A*quadgk(x_int_D,x1,x2)[1]
    end

    return velocity
end
