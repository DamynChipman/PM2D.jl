"""
=========================== 2D PANEL METHOD ====================================
HEADER FILE

@author: Damyn Chipman

PM2D is a 2D Panel Method package. The panel method employs potential flow
theory to solve for element strengths of singularity elements distributed along
a body (such as an airfoil). The tangent flow boundary conditions are imposed
on collocation points distributed along the body. A physical Kutta condition is
also imposed to ensure that the airflow at the trailing edges leaves smoothly.

The Hans-Smith method for element distribution is used. The body has a constant
source distribution with varying strengths at each panel as well as a constant
vortex distribution with a single vortex strength for each panel. This allows
for easy implementation of the Kutta condition by introducing one more unknown
strength to solve for then panels.
"""
module PM2D

# Imports ======================================================================
using QuadGK

# Global Variables =============================================================
const DEG2RAD = pi/180

# Headers ======================================================================
files = ["Panel2D","CalcVelocity","CalcStrengths",
         "InducedVelocity","NACA","Coefs","BoundaryLayer","CalcVortexSheet"]
for header_name in files
    include("PM2D_"*header_name*".jl")
end

# Helper Functions =============================================================

"""
    `unpack_oper_cond(oper_cond)`

Function to unpack operation conditions and return a vector for the free-
stream velocity in the global reference frame.

# ARGUMENTS
* `oper_cond::Array{Float64}`     : Operation conditions, [U_INF,ALPHA]

# OUTPUTS
* `u_inf::Array{Float64}`         : Vector of free stream velocity in
                                    global frame
"""
function unpack_oper_cond(oper_cond)
    U_INF = oper_cond[1]
    AoA = oper_cond[2]
    alpha = AoA*DEG2RAD
    u_inf = U_INF.*[cos(alpha);sin(alpha)]
    return u_inf
end

"""
    `array_to_matrix(array)`

Converts a m-by-n array into an m element array of arrays (matrix)

# ARGUMENTS
* `array::Array`         : MxN array

# OUTPUTS
* `matrix::Array{Array}` : M element array of arrays (matrix)
"""
function array_to_matrix(array)
    matrix = [array[1,:]]
    for n=2:length(array[:,1])
        push!(matrix,array[n,:])
    end
    return matrix
end

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

"""
    `PanelMethod(X_body,Y_body,oper_cond,outputs)`

Wrapper function for the panel method implemented by PM2D. Performs the panel
method with the default options. These options are:

rC_location = 0.5 : Collocation points location on panels
rC_offset = 1e-4  : Offset from panel in normal direction
refine_TE = false : Option to refine collocation points on trailing edges

# ARGUMENTS
* `X_body::Array{Float64}`         : X coordinates of body in CCW direction
* `Y_body::Array{Float64}`         : Y coordinates of body in CCW direction
* `oper_cond::Array{Float64}`      : Operation conditions, [U_INF,ALPHA]
* `outputs::Array{String}`         : Options for outputs of panel method.
                                     Options include:
                                     "Panels", "Coef Matrix", "RHS Vector",
                                     "Strengths", "Pres Coef", "Lift Coef"
"""
function PanelMethod(X_body::Array{Float64},Y_body::Array{Float64},
                     oper_cond::Array{Float64},
                     outputs::Array{String},
                     method::String="Hess-Smith")

    panels = NACA_body(X_body,Y_body)
    NPAN = length(panels)

    to_return = []
    options = ["Panels",
               "Coef Matrix",
               "RHS Vector",
               "Strengths",
               "Pres Coef",
               "Lift Coef",
               "Vortex Sheet"]
    for out in outputs
        if out == options[1]

            push!(to_return,panels)

        elseif out == options[2]

            A,b,strengths = CalcStrengths(panels,oper_cond)
            push!(to_return,A)

        elseif out == options[3]

            A,b,strengths = CalcStrengths(panels,oper_cond)
            push!(to_return,b)

        elseif out == options[4]

            A,b,strengths = CalcStrengths(panels,oper_cond)
            push!(to_return,strengths)

        elseif out == options[5]

            A,b,strengths = CalcStrengths(panels,oper_cond)
            Cp,rC = CalcCp(panels,strengths,oper_cond)
            push!(to_return,[Cp,rC])

        elseif out == options[6]

            A,b,strengths = CalcStrengths(panels,oper_cond)
            Cl = CalcCl(panels,strengths,oper_cond)
            push!(to_return,Cl)

        elseif out == options[7]
            #println("In here")

            RHS = zeros(NPAN-1)
            for i=1:NPAN-1
                RHS[i] = dot(2,panels[i].t_hat,1,unpack_oper_cond(oper_cond),1)
            end
            push!(to_return,CalcVortexSheetCoef(panels,RHS))

        else

            error("Invalid output option for PanelMethod. See ? for options.")

        end
    end

    return to_return
end

#=============================== END OF MODULE ================================#
end
