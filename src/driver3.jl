using PM2D

visualize = false
if visualize
    using VTKtools
    vtk = VTKtools
end

println(" -- Creating Airfoil --")
upper,lower,camber = PM2D.NACA_airfoil("2412",1.0,11)
surf_panels = PM2D.NACA_panels(upper,lower,"source")
camb_panels = PM2D.NACA_camber(camber,"vortex")

panels = [surf_panels; camb_panels]
N_panels = length(panels)
u_inf = [-10.0; 0.0; 0.0]

# Calculate the strengths (Panel Method)
println(" -- Calculating Strengths --")
strengths,A,b = PM2D.calc_strengths(panels,u_inf[1:2],false)

for n=1:N_panels
    println("Panel ",n," --------------")
    println(" Type: ",panels[n].panel_type)
    println(" Location: ",round(panels[n].start_pt,2),"   ",round(panels[n].end_pt,2))
    println(" n_hat: ",round(panels[n].n_hat[1],2),"   ",round(panels[n].n_hat[2],2))
    println(" t_hat: ",round(panels[n].t_hat[1],2),"   ",round(panels[n].t_hat[2],2))
    println(" Strength = ",strengths[n])
end

# Visualize results
if visualize
    println(" -- Visualizing Results --")
    org_points = vcat([[ [side[i,1], side[i,2], 0] for i in 1:size(side,1)] for side in [upper, lower]]...)
    line = org_points
    lines = [line]
    points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
    vtk.generateVTK("airfoil", points; lines=vtk_lines)

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

    println(" Opening ParaView... ")
    run(`paraview --data="mygrid.vtk;airfoil.vtk"`)
end
