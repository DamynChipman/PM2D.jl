#=========================== 2D Panel Method Module ============================
#   HEADER FILE
#
#   @author: Damyn Chipman
#
#   Employs potential flow theory to solve Laplace's equation
#   around an object such as an airfoil. Solutions to Laplace's equation
#   include sources, doublets and vortices. These solution elements are
#   gathered along panels, or discrete lines along the airfoil. By including
#   the interactions between all other panels, the velocity profile can be
#   found around the airfoil.
#
=#
module PM2D

# Imports ======================================================================

# Global Variables =============================================================
const DEG2RAD = pi/180

# Headers ======================================================================
for header_name in ["Panel2D",
                    "CalcVelocity",
                    "CalcStrengths",
                    "InducedVelocity",
                    "NACA",
                    "Coefs"]

    include("PM2D_"*header_name*".jl")
end

# Helper Functions =============================================================

#===============================================================================
# unpack_oper_cond Function
#
#   Function to unpack operation conditions and return a vector for the free-
#   stream velocity in the global reference frame.
#
#   unpack_oper_cond(oper_cond) -> u_inf
#       Inputs:
#       oper_cond - operation conditions = [U_INF,AoA]
#
#       Outputs:
#       u_inf - vector of free-stream velocity
# =============================================================================#
function unpack_oper_cond(oper_cond)
    U_INF = oper_cond[1]
    AoA = oper_cond[2]
    alpha = AoA*DEG2RAD
    u_inf = U_INF.*[cos(alpha);sin(alpha)]
    return u_inf
end

#===============================================================================
# array_to_matrix Function
#
#   Converts a m-by-n array into an m element array of arrays (matrix)
#
#   array_to_matrix(array) -> matrix
#       Inputs:
#       array - m-by-n array
#
#       Outputs:
#       matrix - m element array of arrays (matrix)
# =============================================================================#
function array_to_matrix(array)
    matrix = [array[1,:]]
    for n=2:length(array[:,1])
        push!(matrix,array[n,:])
    end
    return matrix
end

#===============================================================================
# read_dat_file Function
#
#   Reads a .dat or .txt file of airfoil data from airfoiltools.com in "selig"
#   data format (CCW orientation starting at UTE)
#
#   read_dat_file() ->
#       Inputs:
#
#       Outputs:
#
# =============================================================================#
function readcontour(file_name; header_len=1, delim=" ", path="",
                      output="arrays")
  x, y = Float64[], Float64[]

  open(joinpath(path,file_name)) do f
    for (i,line) in enumerate(eachline(f))

      # Ignores header
      if i<=header_len
        nothing
      # Parses each line
      else
        this_x, this_y = split(line, delim; keep=false)
        push!(x, parse(Float64, this_x))
        push!(y, parse(Float64, this_y))
      end

    end
  end

  if output=="arrays"
    return x,y
  elseif output=="matrix"
    xy = [x,y]
    return [xy[j][i] for i in 1:size(x,1), j in 1:2]
  else
    error("Invalid `output` argument $(output).")
  end
end

#===============================================================================
# PanelMethod Function
#
#   Wrapper function for the default options for the implemented panel method.
#   Takes a matrix of points outlining the body in a counter-clockwise
#   orientation of the form:
#   XY = [[X1, Y1],
#         [X2, Y2],
#         [...],
#         [XN, YN]]
#   as well as the operation conditions.
#   Returns a vector of strengths corresponding to the strengths of the elements
#   for each panel.
#
#   PanelMethod(XY,oper_cond) -> strengths
#       Inputs:
#       XY - matrix of points outlining body in CCW orientation
#       oper_cond - operation conditions = [U_INF,AoA]
#
#       Outputs:
#       strengths - list of strengths corresponding to the strengths of the
#                   elements
# =============================================================================#
function PanelMethod(XY,oper_cond)
    X = zeros(length(XY))
    Y = zeros(length(XY))
    for n=1:length(XY)
        X[n] = XY[n][1]
        Y[n] = XY[n][2]
    end
    panels = NACA_body(X,Y)
    A,b,strengths = CalcStrengths(panels,oper_cond)
    return strengths,panels
end

#=============================== END OF MODULE ================================#
end
