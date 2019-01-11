"""
    `Panel2D(r1,r1,rC_location,rC_offset,orientation)`

Object representing the geometry of a 2D panel. Points are in a global reference frame.

# ARGUMENTS
* `r1::Array{Float64}`        : Starting point of panel
* `r2::Array{Float64}`        : Ending point of panel
* `rC_location::Float64=0.5`  : Location along panel (between 0 and 1) of collocation point
* `rC_offset::Float64=0.01`   : Offset of collocation point from panel in direction of normal vector
* `orientation::String`       : Panel (and collection of panels) orientation; either "CCW" or "CW"

# PROPERTIES
* `r1::Array{Float64}`        : Starting point of panel
* `r2::Array{Float64}`        : Ending point of panel
* `rC::Array{Float64}`        : Collocation point along panel
* `L::Float64`                : Length of panel
* `theta::Float64`            : Angle of panel relative to global X-axis
* `n_hat::Array{Float64}`     : Normal vector
* `t_hat::Array{Float64}`     : Tangent vector
"""
type Panel2D

    # -- Properties --
    r1::Array{Float64}
    r2::Array{Float64}
    rC::Array{Float64}
    L::Float64
    theta::Float64
    n_hat::Array{Float64}
    t_hat::Array{Float64}
    rC_offset::Float64

    # -- Constructor --
    function Panel2D(r1::Array{Float64},
                     r2::Array{Float64};
                     rC_location::Float64=0.5,
                     rC_offset::Float64=1e-4,
                     orientation::String="CCW")

        x1, y1 = r1[1], r1[2]
        x2, y2 = r2[1], r2[2]

        L = sqrt((y2 - y1)^2 + (x2 - x1)^2)
        theta = atan2((y2 - y1),(x2 - x1))

        t_hat = (r2 - r1)/norm(r2 - r1)
        n_hat = [t_hat[2],-t_hat[1]]

        xC = (rC_location*L*cos(theta)) + x1 + rC_offset*n_hat[1]
        yC = (rC_location*L*sin(theta)) + y1 + rC_offset*n_hat[2]

        new(r1,r2,[xC,yC],L,theta,n_hat,t_hat,rC_offset)
    end
end
