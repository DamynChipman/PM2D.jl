module PM2D
#= ------------------------- 2D Panel Method Module ----------------------------
#   @author: Damyn Chipman
#
#   This package employs potential flow theory to solve Laplace's equation
#   around an object such as an airfoil. Solutions to Laplace's equation
#   include sources, doublets and vortices. These solution elements are
#   gathered along panels, or discrete lines along the airfoil. By including
#   the interactions between all other panels, the velocity profile can be
#   found around the airfoil.
#
#
=#



# Begin Class Definitions ------------------------------------------------------

#===============================================================================
# Panel2D Type
#
#   Represents a 2-D Panel given the starting and end (geometry) point of the
#   panel
#
#   Panel2D(star_pt,end_pt)
#       Inputs:
#       start_pt - (x,y) location of starting point
#       end_pt - (x,y) location of ending point
#
#       Properties:
#       start_pt::Array{Float64} - (x,y) location of starting point
#       end_pt::Array{Float64} - (x,y) location of ending point
#       r0::Array{Float64} - (x,y) location of element point
#       rC::Array{Float64} - (x,y) location of collocation point (BC)
#       L::Float64 - length of panel
#       theta::Float64 - angle of panel relative to x-axis
#       n_hat::Array{Float64} - normal vector
#       t_hat::Array{Float64} - tanget vector
# =============================================================================#
type Panel2D

    # -- Properties --
    start_pt::Array{Float64}
    end_pt::Array{Float64}
    r0::Array{Float64}
    rC::Array{Float64}
    L::Float64
    theta::Float64
    n_hat::Array{Float64}
    t_hat::Array{Float64}

    # -- Constructor --
    function Panel2D(start_pt,end_pt)

        x1, y1 = start_pt[1], start_pt[2]
        x2, y2 = end_pt[1], end_pt[2]

        L = sqrt((y2 - y1)^2 + (x2 - x1)^2)
        theta = atan2((y2 - y1),(x2 - x1))

        x0 = (.25*L*cos(theta)) + x1
        y0 = (.25*L*sin(theta)) + y1
        xC = (.75*L*cos(theta)) + x1
        yC = (.75*L*sin(theta)) + y1

        n_hat = [sin(theta), cos(theta)]
        t_hat = (end_pt - start_pt)/norm(end_pt - start_pt)

        new(start_pt,end_pt,[x0,y0],[xC,yC],L,theta,n_hat,t_hat)

    end
end

# Begin Function Definitions ---------------------------------------------------

#===============================================================================
# NACA_airfoil Function
#
#   Creates two lists of ordered pairs representing a NACA Four
#   digit airfoil
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
function NACA_airfoil(numb::Int64,c=1)

    # Unpackage the digits
    m = digits(numb)[4]*.01
    p = digits(numb)[3]*.1
    tau = (digits(numb)[2]*10 + digits(numb)[1])*.01

    # Define camber geometery
    x1 = 0:0.01:(p*c)
    Y1 = ((m)./(p^2)).*(2*p*(x1./c) - (x1./c).^2)
    dY_dx1 = ((2*m)/(p^2)).*(p - x1./c)

    x2 = (p*c):0.01:(c)
    Y2 = ((m)./((1 - p)^2)).*((1 - 2*p) + (2*p).*(x2./c) - (x2./c).^2)
    dY_dx2 = ((m)/((1 - p)^2)).*((1 - 2*p) + 2*p.*(x2./c) - (x2./c).^2)

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

    # Define airfoil thickness
    T(x) = 5*tau*(0.2969*x^.5 - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)

    # Determine x and y coordinates of airfoil surface
    x_upper, x_lower, y_upper, y_lower = 0.*x, 0.*x, 0.*x, 0.*x
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

    return upper,lower
end

#===============================================================================
# NACA_panels Function
#
#   Creates a list of Panel2D objects corresponding to a NACA airfoil with
#   specified upper and lower coordinates from NACA_airfoil function
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
function NACA_panels(upper,lower)
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

#===============================================================================
# calc_potential Function
#
#
#
#
#
# =============================================================================#
function calc_potential(panel_type::String,str,elem_pt,coll_pt)
    if panel_type == "source"
        potential = (str/(2*pi))*log(sqrt(norm(elem_pt - coll_pt)))
    elseif panel_type == "vortex"
        potential = -(str/(2*pi))*(atan((elem_pt[2] - coll_pt[2])/
                                        (elem_pt[1] - coll_pt[1])))
    elseif panel_type == "doublet"
        potential = -(1/(2*pi))*((dot(strength,elem_pt))/
                                 (norm(elem_pt - coll_pt)))
    else
        error("Invalid panel_type. Options are: source, vortex, or doublet.")
    end
    return potential
end

#===============================================================================
# calc_velocity Function
#
#
#
#
#
# =============================================================================#
function calc_velocity(panel_type::String,str,elem_pt,coll_pt)
    velocity = Array{Float64}(2)
    if panel_type == "source"
        velocity[1] = (str/(2*pi))*((elem_pt[1] - coll_pt[1])/
                    ((elem_pt[1] - coll_pt[1])^2 + (elem_pt[2] - coll_pt[2])^2))
        velocity[2] = (str/(2*pi))*((elem_pt[2] - coll_pt[2])/
                    ((elem_pt[1] - coll_pt[1])^2 + (elem_pt[2] - coll_pt[2])^2))
    elseif panel_type == "vortex"
        velocity[1] = (str/(2*pi))*((elem_pt[2] - coll_pt[2])/
                                    (norm(elem_pt - coll_pt)^2))
        velocity[2] = (-str/(2*pi))*((elem_pt[1] - coll_pt[1])/
                                    (norm(elem_pt - coll_pt)^2))
    elseif panel_type == "doublet"
        error("Still in development")
    else
        error("Invalid panel_type. Options are: source, vortex, or doublet.")
    end
    return velocity
end

#===============================================================================
# calc_strengths Function
#
#   Given a callable velocity vector function and a list of Panel2D objects,
#   calculates the strengths of the elements. If doKutta is true, inserts the
#   Kutta condition row on the second row of the coefficient matrix. Solves
#   using b\A notation.
#
#   calc_strengths(f_vel,panels,u_inf,doKutta) -> strengths
#       Inputs:
#       f_vel - callable velocity vector function. Should return a vector.
#       panels - list of Panel2D objects corresponding to the locations of
#                panels
#       u_inf - Operation condition for the free stream velocity
#       doKutta - boolean to determine if Kutta condition should be implemented
#       Outputs:
#       strengths - list of strengths corresponding to the strengths of the
#                   elements
# =============================================================================#
function calc_strengths(f_vel,panels,u_inf,doKutta)
    N_panels = length(panels)

    A = Array{Float64}(N_panels,N_panels)
    for i=1:N_panels
        for j=1:N_panels
            A[i,j] = dot(f_vel(panels[i].r0,panels[j].rC),panels[i].n_hat)
        end
    end

    b = Array{Float64}(N_panels)
    for i=1:N_panels
        b[i] = -dot(u_inf,panels[i].n_hat)
    end

    if doKutta == true
        # Implement the Kutta Condition

        kutta = Array{Float64}(N_panels)
        for j=1:N_panels
            kutta[j] = dot(f_vel(panels[1].r0,panels[j].rC),panels[1].t_hat) +
                       dot(f_vel(panels[N_panels].r0,panels[j].rC),panels[N_panels].t_hat)
        end
        A[2,:] = kutta
        b[2] = 0.0
        strengths = b\A
    else
        strengths = b\A
    end
    return strengths
end

end
