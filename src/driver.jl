import PM2D
import VTKtools
vtk = VTKtools

pan1 = PM2D.Panel2D([0.0 0.0],[1.0 0.0])
pan2 = PM2D.Panel2D([1.0 0.0],[2.0 0.0])
pan3 = PM2D.Panel2D([2.0 0.0],[3.0 0.0])
pan4 = PM2D.Panel2D([3.0 0.0],[4.0 0.0])
panels = [pan1, pan2, pan3, pan4]

# upper,lower = PM2D.NACA_airfoil("2412",1.0,1000)
# panels = PM2D.NACA_panels(upper,lower)
N_panels = length(panels)

function u_vortex(r0,rC)
    u = PM2D.calc_velocity("vortex",1.0,r0,rC)
    return u
end

u_inf = [-30.0;5.0;0.0]
strengths,A,b = PM2D.calc_strengths(u_vortex,panels,u_inf[1:2],false)

# Visualize results
# org_points = vcat([[ [side[i,1], side[i,2], 0] for i in 1:size(side,1)] for side in [upper, lower]]...)
# line = org_points
# lines = [line]
# points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
# vtk.generateVTK("airfoil", points; lines=vtk_lines)

function wrap_u(X)
    u = u_inf
    for j=1:N_panels
        u = u + vcat(PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,X),0)
    end
    return u
end

P_min = [-.5,-.4]
P_max = [1.5,.4]
NDIVS = [200,100]
fgrid = vtk.Grid(P_min, P_max, NDIVS)

vtk.calculate_field(fgrid, wrap_u, "U", "vector")
vtk.save(fgrid, "mygrid")

run(`paraview --data="mygrid.vtk"`)

# upper,lower = PM2D.NACA_airfoil(2412)
# panels = PM2D.NACA_panels(upper,lower)
#
# function u_vortex(r0,rC)
#     u = PM2D.calc_velocity("vortex",1.0,r0,rC)
#     return u
# end
#
# strengths = PM2D.calc_strengths(u_vortex,panels,[-5.0;0.0],true)
# println(strengths)


#
# function f_vel(r0,r)
#     vel = [0.0; 0.0]
#     vel[1] = (1/(2*pi))*((r[2] - r0[2])/(norm(r - r0)^2))
#     vel[2] = -(1/(2*pi))*((r[1] - r0[1])/(norm(r - r0)^2))
#     return vel
# end
#
# str = PM2D.calc_strengths(f_vel,panels,5.*[sin(pi/12);cos(pi/12)],true)
# println(str)
