#===============================================================================
# Panel2D Type
#
#   Represents a 2-D Panel given the starting and end (geometry) point of the
#   panel
#
#   Panel2D(star_pt,r2)
#       Inputs:
#       r1 - (x,y) location of starting point
#       r2 - (x,y) location of ending point
#
#       Properties:
#       r1::Array{Float64} - (x,y) location of starting point
#       r2::Array{Float64} - (x,y) location of ending point
#       r0::Array{Float64} - (x,y) location of element point
#       rC::Array{Float64} - (x,y) location of collocation point (BC)
#       L::Float64 - length of panel
#       theta::Float64 - angle of panel relative to x-axis
#       n_hat::Array{Float64} - normal vector
#       t_hat::Array{Float64} - tanget vector
# =============================================================================#
type Panel2D

    # -- Properties --
    r1::Array{Float64}
    r2::Array{Float64}
    r0::Array{Float64}
    rC::Array{Float64}
    L::Float64
    theta::Float64
    n_hat::Array{Float64}
    t_hat::Array{Float64}

    # -- Constructor --
    function Panel2D(r1::Array{Float64},
                     r2::Array{Float64},
                     r0_location::Float64=0.25,
                     rC_location::Float64=0.5,
                     rC_offset::Float64=0.01,
                     orientation::String="CCW")

        x1, y1 = r1[1], r1[2]
        x2, y2 = r2[1], r2[2]

        L = sqrt((y2 - y1)^2 + (x2 - x1)^2)
        theta = atan2((y2 - y1),(x2 - x1))

        t_hat = (r2 - r1)/norm(r2 - r1)
        if orientation == "CW"
            n_hat = [t_hat[2],-t_hat[1]]
        elseif orientation == "CCW"
            n_hat = -[t_hat[2],-t_hat[1]]
        else
            options = ["CCW, ","CW"]
            error("Invalid orientation for Panel2D. Options are: "*options)
        end

        x0 = (r0_location*L*cos(theta)) + x1
        y0 = (r0_location*L*sin(theta)) + y1
        xC = (rC_location*L*cos(theta)) + x1 + rC_offset*n_hat[1]
        yC = (rC_location*L*sin(theta)) + y1 + rC_offset*n_hat[2]

        new(r1,r2,[x0,y0],[xC,yC],L,theta,n_hat,t_hat)
    end
end
