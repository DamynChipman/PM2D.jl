"""
    `NACA_airfoil(numb, N, c)`

NACA Four-Digit Series Airfoil. Returns two lists of ordered pairs representing
a NACA Four-Digit Airfoil.

# ARGUMENTS
* `numb::String`     : Four digit series. Symmetric airfoil given by "00xx"
* `N::Int64`         : Number of points
* `c::Float64=1.0`   : Length of airfoil. Defaults to 1.0 for appropiate scaling
"""
function NACA_airfoil(numb::String,
                      N::Int64,
                      c::Float64=1.0)

    # Unpackage the digits
    m = parse(Int,numb[1])*.01
    p = parse(Int,numb[2])*.1
    tau = (parse(Int,numb[3])*10 + parse(Int,numb[4]))*.01

    println(" NACA digits: m = ",m,"  p = ",p,"  tau = ",tau)

    # Define camber geometery
    if p != 0
        del_x1 = (p*c - 0)/(N/4)
        del_x2 = (c - p*c)/(N/4)

        x1 = 0:del_x1:(p*c)
        Y1 = ((m)./(p^2)).*(2*p*(x1./c) - (x1./c).^2)
        dY_dx1 = ((2*m)/(p^2)).*(p - x1./c)

        x2 = (p*c):del_x2:(c)
        Y2 = ((m)./((1 - p)^2)).*((1 - 2*p) + (2*p).*(x2./c) - (x2./c).^2)
        dY_dx2 = ((m)/((1 - p)^2)).*((1 - 2*p) + 2*p.*(x2./c) - (x2./c).^2)

        camber = [x1 Y1; x2 Y2]

        x = zeros(length(x1)+length(x2)-1)
        Y = zeros(length(Y1)+length(Y2)-1)
        dY_dx = zeros(length(dY_dx1)+length(dY_dx2)-1)
        for i=1:length(x1)
            x[i] = x1[i]
            Y[i] = Y1[i]
            dY_dx[i] = dY_dx1[i]
        end
        for i=2:length(x2)
            x[i+length(x1)-1] = x2[i]
            Y[i+length(x1)-1] = Y2[i]
            dY_dx[i+length(x1)-1] = dY_dx2[i]
        end
    else
        del_x = c/(N/2)
        x = 0:del_x:c
        Y = zeros(length(x),1)
        camber = [x Y]
        dY_dx = zeros(length(x),1)
    end

    # Define airfoil thickness
    T(x) = 5*tau*(0.2969*x^.5 - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)

    # Determine x and y coordinates of airfoil surface
    x_upper, x_lower, y_upper, y_lower = 0. * x, 0. * x, 0. * x, 0. * x
    for i=1:length(x)
        theta = atan(dY_dx[i])

        x_upper[i] = x[i] - T(x[i])*sin(theta)
        x_lower[i] = x[i] + T(x[i])*sin(theta)
        y_upper[i] = Y[i] + T(x[i])*cos(theta)
        y_lower[i] = Y[i] - T(x[i])*cos(theta)
    end

    upper = zeros(length(x_upper),2)
    lower = zeros(length(x_lower),2)
    for i=1:length(upper[:,1])
        upper[i,:] = [x_upper[i],y_upper[i]]
        lower[i,:] = [x_lower[i],y_lower[i]]
    end

    if m == 0 && p == 0
        upper[1,:] = [0.0 0.0]
        lower[1,:] = [0.0 0.0]
    end

    # Redefine the orientation of the airfoil
    airfoil = zeros(2*length(x_lower)-2,2)
    for n=1:length(upper[:,1])
        airfoil[n,1] = upper[(length(upper[:,1])+1)-n,1]
        airfoil[n,2] = upper[(length(upper[:,1])+1)-n,2]
    end
    for n=2:(length(upper[:,1])-1)
        airfoil[n+length(upper[:,1])-1,1] = lower[n,1]
        airfoil[n+length(upper[:,1])-1,2] = lower[n,2]
    end

    return upper,lower,airfoil,camber
end

#===============================================================================
# NACA_panels Function
#
#   Creates a list of Panel2D objects corresponding to a NACA airfoil with
#   specified upper and lower coordinates from NACA_airfoil function
#
#   Note: upper and lower should "start left and go right", i.e. form
#         a line beginning with small x values towards larger x values
#
#   NACA_panels(upper,lower) -> panels
#       Inputs:
#       upper - list of (x,y) corresponding to the positions on the upper side
#               of an airfoil
#       lower - list of (x,y) corresponding to the positions on the lower side
#               of an airfoil
#       Outputs:
#       panels - list of Panel2D objects for the given arifoil. Starts at the
#                upper trailing edge and ends with the lower trailing edge
# =============================================================================#
function NACA_panels(upper::Array{Float64},
                     lower::Array{Float64},
                     panel_type::String)

    N_panels = (length(upper[:,1]) - 1) + (length(lower[:,1]) - 1)
    N_side = Int64(N_panels/2)
    N_pts = length(upper[:,1])
    panels = Array{Panel2D}(N_panels)
    for n=0:(N_side-1)
        panels[n+1] = Panel2D(upper[N_pts-n,:],upper[N_pts-n-1,:])
        panels[n+1+N_side] = Panel2D(lower[n+1,:],lower[n+2,:])
    end
    return panels
end

function NACA_camber(camber::Array{Float64},
                     panel_type::String)

    N_panels = length(camber[:,1]) - 1
    panels = Array{Panel2D}(N_panels)
    for n=1:N_panels
        panels[n] = Panel2D(camber[n,:],camber[n+1,:])
    end
    return panels
end

"""
    `NACA_body(X,Y,orientation,r0_location,rC_location,rC_offset,refine_TE)`

Creates an array of Panel2D objects based on a closed body given by a set of
X and Y points. Can be in either clockwise of counterclockwise orientations but
must be a closed, continuous body.

# ARGUMENTS
* `X::Array{Float64}`           : List of X body points
* `Y::Array{Float64}`           : List of Y body points
* `orientation::String="CCW"`   : Panel orientation; either "CCW" or "CW"
* `rC_location::Float64=0.5`    : Location along panel (between 0 and 1) of collocation point
* `rC_offset::Float64=0.01`     : Offset of collocation point from panel in direction of normal vector
* `refine_TE::Bool=false`       : Option to place collocation points of trailing edges closer to trailing edge

# OUTPUTS
* `panels::Array{Panel2D}`      : Array of Panel2D objects representing a body
"""
function NACA_body(X::Array{Float64},
                   Y::Array{Float64};
                   orientation::String="CCW",
                   rC_location::Float64=0.5,
                   rC_offset::Float64=1e-4,
                   refine_TE=false)

    # Iterate around surface to create panels
    NPTS = length(X)
    NPAN = NPTS - 1
    panels = Array{Panel2D}(NPAN)
    for n=1:NPTS-1
        if refine_TE
            if n == 1
                rC_location = 0.1
            elseif n == NPAN
                rC_location = 0.9
            end
            r1 = [X[n];Y[n]]
            r2 = [X[n+1];Y[n+1]]
            panels[n] = Panel2D(r1,r2,rC_location=rC_location,
                                      rC_offset=rC_offset,
                                      orientation=orientation)
        else
            r1 = [X[n];Y[n]]
            r2 = [X[n+1];Y[n+1]]
            panels[n] = Panel2D(r1,r2,rC_location=rC_location,
                                      rC_offset=rC_offset,
                                      orientation=orientation)
        end
    end
    return panels
end
