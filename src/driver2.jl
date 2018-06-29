using PM2D

visualize = true
if visualize
    using VTKtools
    vtk = VTKtools
end

# Create one source panel for analysis
println(" -- Creating panel -- ")
source_pan = PM2D.Panel2D("source",[0.0; 0.0],[1.0; 0.0])

# Function for velocity
pan_strgth = 5.0
function s_vel(r,r0)
    velocity = PM2D.calc_velocity(source_pan.panel_type,pan_strgth,r0,r)
end

u_inf = [-10.0; 0.0; 0.0]

# Visualize results
# println(" -- Generating Grid -- ")
# org_points = [0.0 0.0; 1.0 0.0]
# line = org_points
# lines = [line]
# points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
# vtk.generateVTK("airfoil", points; lines=vtk_lines)

function wrap_u(X)
    u = u_inf
    u = u + vcat(s_vel(X,source_pan.r0),0)
    return u
end

println(" -- Calculating Field -- ")
P_min = [-1.0,-1.0]
P_max = -P_min
NDIVS = [100,100]
fgrid = vtk.Grid(P_min, P_max, NDIVS)

vtk.calculate_field(fgrid, wrap_u, "U", "vector")
vtk.save(fgrid, "mygrid")

println(" Opening ParaView... ")
run(`paraview --data="mygrid.vtk"`)
