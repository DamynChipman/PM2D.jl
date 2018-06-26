# NACA2412 Session
using PM2D
using VTKtools
vtk = VTKtools

# Create a NACA2412 Airfoil and create a list of Panel2D type panels
println(" -- Creating Airfoil --")
upper,lower = PM2D.NACA_airfoil("2412",1.0,400)
panels = PM2D.NACA_panels(upper,lower)
N_panels = length(panels)

# Create wrapper velocity function for calc_strength function
function u_vortex(r0,rC)
    u = PM2D.calc_velocity("source",1.0,r0,rC)
    return u
end

# Operation conditions
u_inf = [-10.0;0.0;0.0]

# Calculate the strengths (Panel Method)
println(" -- Calculating Strengths --")
strengths,A,b = PM2D.calc_strengths(u_vortex,panels,u_inf[1:2],true)

# Visualize results
println(" -- Visualizing Results --")
org_points = vcat([[ [side[i,1], side[i,2], 0] for i in 1:size(side,1)] for side in [upper, lower]]...)
line = org_points
lines = [line]
points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
vtk.generateVTK("airfoil", points; lines=vtk_lines)

function wrap_u(X)
    u = u_inf
    for j=1:N_panels
        u = u + vcat(PM2D.calc_velocity("source",strengths[j],panels[j].r0,X),0)
    end
    return u
end

P_min = [-.5,-.4]
P_max = [1.5,.4]
NDIVS = [200,100]
fgrid = vtk.Grid(P_min, P_max, NDIVS)

vtk.calculate_field(fgrid, wrap_u, "U", "vector")
vtk.save(fgrid, "mygrid")

println(" Opening ParaView... ")
run(`paraview --data="mygrid.vtk;airfoil.vtk"`)
