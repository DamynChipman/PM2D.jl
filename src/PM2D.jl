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
# Source2D Type
#
#   Represents a 2-D Point Source Singularity element
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
type Source2D

    # -- Properties --
    strength::Float64
    elem_loc::Array{Float64}
    coll_loc::Array{Float64}
    potential::Float64
    velocity::Array{Float64}

    # -- Constructor --
    function Source2D(strength,elem_loc)
        # Error Handling:

        new(strength,elem_loc)
    end
end

# --- calc_potential Method --
#
#
#
#
function _calc_potential(obj::Source2D,coll_point)
    # Calculate potential
    obj.potential = (obj.strength/(2*pi))*
                    log(sqrt(norm(obj.elem_loc - coll_point)))
end

# -- calc_velocity Method --
#
#
#
#
function _calc_velocity(obj::Source2D,coll_point)
    # Initialize velocity vector
    obj.velocity = [0;0]

    # Compute x and z velocity components
    obj.velocity[1] = (obj.strength/(2*pi))*
                      ((obj.elem_loc[1] - coll_point[1])/
                       ((obj.elem_loc[1] - coll_point[1])^2 +
                        (obj.elem_loc[2] - coll_point[2])^2))
    obj.velocity[2] = (obj.strength/(2*pi))*
                      ((obj.elem_loc[2] - coll_point[2])/
                       ((obj.elem_loc[1] - coll_point[1])^2 +
                        (obj.elem_loc[2] - coll_point[2])^2))
end

#===============================================================================
# Vortex2D Type
#
#   Represents a 2-D Point Vortex Singularity element
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
type Vortex2D

    # -- Properties --
    strength::Array{Float64}
    elem_loc::Array{Float64}
    coll_loc::Array{Float64}
    potential::Float64
    velocity::Array{Float64}

    # -- Constructor --
    function Vortex2D(strength,elem_loc)
        # Error Handling:

        new(strength,elem_loc)
    end
end

# -- calc_potential Method --
#
#
#
#
function _calc_potential(obj::Vortex2D,coll_point)
    # Calculate potential
    obj.potential = -(obj.strength/(2*pi))*
                    (atan((obj.elem_loc[2] - coll_point[2])/
                         (obj.elem_loc[1] - coll_point[1])))
end

# -- calc_velocity Method --
#
#
#
#
function _calc_velocity(obj::Vortex2D,coll_point)
    # Initialize velocity vector
    obj.velocity = [0;0]

    # Compute x and z velocity components
    obj.velocity[1] = (obj.strength/(2*pi))*
                      ((obj.elem_loc[2] - coll_point[2])/
                       (norm(obj.elem_loc - coll_point)^2))
   obj.velocity[2] = (-obj.strength/(2*pi))*
                     ((obj.elem_loc[1] - coll_point[1])/
                      (norm(obj.elem_loc - coll_point)^2))
end

#===============================================================================
# Doublet2D Type
#
#   Represents a 2-D Point Doublet Singularity element
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
type Doublet2D

    # -- Properties --
    strength::Array{Float64}
    elem_loc::Array{Float64}
    coll_loc::Array{Float64}
    potential::Float64
    velocity::Array{Float64}

    # -- Constructor --
    function Doublet2D(strength,elem_loc)
        # Error Handling:

        new(strength,elem_loc)
    end
end

# -- calc_potential Method --
#
#
#
#
function _calc_potential(obj::Doublet2D,coll_point)
    # Calculate potential
    obj.potential = -(1/(2*pi))*
                    ((dot(obj.strength,obj.elem_loc))/
                     (norm(obj.elem_loc - coll_point)))
end

# -- calc_velocity Method --
#
#
#
#
function _calc_velocity(obj::Doublet2D,coll_point)
    # Initialize velocity vector
    obj.velocity = [0;0]

    # Compute x and z velocity components

end

#===============================================================================
# Panel2D Type
#
#   Represents a 2-D Panel
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
type Panel2D

    # -- Properties --
    panel_type::String
    start_pt::Array{Float64}
    end_pt::Array{Float64}
    r0::Array{Float64}
    rC::Array{Float64}
    L::Float64
    theta::Float64

    # -- Constructor --
    function Panel2D(panel_type,start_pt,end_pt)

        x1, y1 = start_pt[1], start_pt[2]
        x2, y2 = end_pt[1], end_pt[2]

        L = sqrt((y2 - y1)^2 + (x2 - x1)^2)
        theta = atan((y2 - y1)/(x2 - x1))

        x0 = (.25*L*cos(theta)) + x1
        y0 = (.25*L*sin(theta)) + y1
        xC = (.75*L*cos(theta)) + x1
        yC = (.75*L*sin(theta)) + y2

        new(panel_type,start_pt,end_pt,[x0,y0],[xC,yC],L,theta)

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
# calc_coeff Function
#
#   Given x and y coordinates of an airfoil, calculate the influence
#   coefficients and reutrn a matrix containing them
#
#   [insert some kind of explanation here]
#
#
# =============================================================================#
function calc_coeff(upper_profile,lower_profile)
    N_panels = length(upper_profile[1,:]) - 1
    upper_airfoil = zeros(N_panels)
    lower_airfoil = zeros(N_panels)
    println(length(upper_profile))
    println(length(N_panels))
    for i=1:(N_panels)
        upper_airfoil[i] = Panel2D("vortex",upper_profile[i],upper_profile[i+1])
        lower_airfoil[i] = Panel2D("vortex",lower_profile[i],lower_profile[i+1])
    end



end


end
