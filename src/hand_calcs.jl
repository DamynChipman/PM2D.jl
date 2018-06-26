# Hand Calcs for a 4-panel surface
# --- Hand Calcs ---
println(" ===================== Hand Calcs ===================== ")

# Define velocity function
function f_vel(x,y,x0,y0)
    velocity = [0.0;0.0]
    velocity[1] = (1/(2*pi))*((y - y0)/((x-x0)^2 + (y-y0)^2))
    velocity[2] = -(1/(2*pi))*((x - x0)/((x-x0)^2 + (y-y0)^2))
    return velocity
end

# Define geometry
x0_panels = [1.75,0.75,0.25,1.25]
y0_panels = [0.25,0.75,-0.25,-0.75]
xC_panels = [1.25,0.25,0.75,1.75]
yC_panels = [0.75,0.25,-0.75,-0.25]

n1_hat = (sqrt(2)/2).*[1;1]
n2_hat = (sqrt(2)/2).*[-1;1]
n3_hat = (sqrt(2)/2).*[-1;-1]
n4_hat = (sqrt(2)/2).*[1;-1]
n_hats = [n1_hat n2_hat n3_hat n4_hat]

t1_hat = (sqrt(2)/2).*[1;-1]
t2_hat = (sqrt(2)/2).*[1;1]
t3_hat = (sqrt(2)/2).*[1;-1]
t4_hat = (sqrt(2)/2).*[1;1]
t_hats = [t1_hat t2_hat t3_hat t4_hat]

# Create coefficient matrix
A = [dot((f_vel(xC_panels[1],yC_panels[1],x0_panels[1],y0_panels[1])),n1_hat) dot((f_vel(xC_panels[1],yC_panels[1],x0_panels[2],y0_panels[2])),n1_hat) dot((f_vel(xC_panels[1],yC_panels[1],x0_panels[3],y0_panels[3])),n1_hat) dot((f_vel(xC_panels[1],yC_panels[1],x0_panels[4],y0_panels[4])),n1_hat);
     dot((f_vel(xC_panels[2],yC_panels[2],x0_panels[1],y0_panels[1])),n2_hat) dot((f_vel(xC_panels[2],yC_panels[2],x0_panels[2],y0_panels[2])),n2_hat) dot((f_vel(xC_panels[2],yC_panels[2],x0_panels[3],y0_panels[3])),n2_hat) dot((f_vel(xC_panels[2],yC_panels[2],x0_panels[4],y0_panels[4])),n2_hat);
     dot((f_vel(xC_panels[3],yC_panels[3],x0_panels[1],y0_panels[1])),n3_hat) dot((f_vel(xC_panels[3],yC_panels[3],x0_panels[2],y0_panels[2])),n3_hat) dot((f_vel(xC_panels[3],yC_panels[3],x0_panels[3],y0_panels[3])),n3_hat) dot((f_vel(xC_panels[3],yC_panels[3],x0_panels[4],y0_panels[4])),n3_hat);
     dot((f_vel(xC_panels[4],yC_panels[4],x0_panels[1],y0_panels[1])),n4_hat) dot((f_vel(xC_panels[4],yC_panels[4],x0_panels[2],y0_panels[2])),n4_hat) dot((f_vel(xC_panels[4],yC_panels[4],x0_panels[3],y0_panels[3])),n4_hat) dot((f_vel(xC_panels[4],yC_panels[4],x0_panels[4],y0_panels[4])),n4_hat)]

println(" -- Coefficient Matrix -- ")
println(A[1,:])
println(A[2,:])
println(A[3,:])
println(A[4,:])

# Create RHS vector
u_inf = [-10.0;0.0]
RHS = -[dot(u_inf,n1_hat);
        dot(u_inf,n2_hat);
        dot(u_inf,n3_hat);
        dot(u_inf,n4_hat)]

println(" -- RHS Vector --")
println(RHS)

# Kutta Condition
doKutta = false
if doKutta
     kutta = Array{Float64}(4)
     for j=1:4
          kutta[j] = dot(f_vel(xC_panels[1],yC_panels[1],x0_panels[j],y0_panels[j]),t1_hat) +
                     dot(f_vel(xC_panels[4],yC_panels[4],x0_panels[j],y0_panels[j]),t4_hat)
     end
     A[2,:] = kutta
     RHS[2] = 0.0

     println(" -- Kutta Condition Imposed --")
     println(" -- Coefficient Matrix -- ")
     println(A[1,:])
     println(A[2,:])
     println(A[3,:])
     println(A[4,:])
     println(" -- RHS Vector --")
     println(RHS)
end

# Solve for strengths
strengths = A\RHS

println(" -- Strengths Vector -- ")
println(strengths)

# Check BC 1
u1_ind,u2_ind,u3_ind,u4_ind = [0.0;0.0],[0.0;0.0],[0.0;0.0],[0.0;0.0]

for j=1:4
     u1_ind = u1_ind + f_vel(xC_panels[1],yC_panels[1],x0_panels[j],y0_panels[j]).*strengths[j]
     u2_ind = u2_ind + f_vel(xC_panels[2],yC_panels[2],x0_panels[j],y0_panels[j]).*strengths[j]
     u3_ind = u3_ind + f_vel(xC_panels[3],yC_panels[3],x0_panels[j],y0_panels[j]).*strengths[j]
     u4_ind = u4_ind + f_vel(xC_panels[4],yC_panels[4],x0_panels[j],y0_panels[j]).*strengths[j]
end

induced = [dot(u1_ind,n1_hat);dot(u2_ind,n2_hat);dot(u3_ind,n3_hat);dot(u4_ind,n4_hat)]
fs = [dot(u_inf,n1_hat);dot(u_inf,n2_hat);dot(u_inf,n3_hat);dot(u_inf,n4_hat)]
diff = [abs(induced[1]) - abs(fs[1]);
        abs(induced[2]) - abs(fs[2]);
        abs(induced[3]) - abs(fs[3]);
        abs(induced[4]) - abs(fs[4])]

println("Difference between Induced and Free-Stream (HC): ",diff)

# Check BC 2

visualize = false
if visualize
     using VTKtools
     vtk = VTKtools

     # Visualize results
     upper = [0.0 0.0; 1.0 1.0; 2.0 0.0]
     lower = [0.0 0.0; 1.0 -1.0; 2.0 0.0]
     org_points = vcat([[ [side[i,1], side[i,2], 0] for i in 1:size(side,1)] for side in [upper, lower]]...)
     line = org_points
     lines = [line]
     points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
     vtk.generateVTK("airfoil", points; lines=vtk_lines)

     function wrap_u(X)
         u = [u_inf[1];u_inf[2];0.0]
         #u = zeros(3,1)
         for j=1:4
             u = u + vcat(f_vel(X[1],X[2],x0_panels[j],y0_panels[j]).*strengths[j],0)
         end
         return u
     end

     P_min = [-2.5,-2.5]
     P_max = [3.5,2.5]
     NDIVS = [100,50]
     fgrid = vtk.Grid(P_min, P_max, NDIVS)

     vtk.calculate_field(fgrid, wrap_u, "U", "vector")
     vtk.save(fgrid, "mygrid")

     run(`paraview --data="mygrid.vtk;airfoil.vtk"`)
end


# PM2D
println(" ======================== PM2D ======================== ")
using PM2D

# Create panels
upper = [0.0 0.0; 1.0 1.0; 2.0 0.0]
lower = [0.0 0.0; 1.0 -1.0; 2.0 0.0]
panels = PM2D.NACA_panels(upper,lower)
N_panels = length(panels)

# Create wrapper velocity function for calc_strength function
function u_vortex(r0,rC)
    u = PM2D.calc_velocity("vortex",1.0,r0,rC)
    return u
end

# Calculate the strengths (Panel Method)
strengths,A,b = PM2D.calc_strengths(u_vortex,panels,u_inf[1:2],false)

println(" -- Coefficient Matrix -- ")
println(A[1,:])
println(A[2,:])
println(A[3,:])
println(A[4,:])

println(" -- RHS Vector --")
println(b)

println(" -- Strengths Vector -- ")
println(strengths)

# Check BC
u1_ind,u2_ind,u3_ind,u4_ind = [0.0;0.0],[0.0;0.0],[0.0;0.0],[0.0;0.0]

for j=1:4
     u1_ind = u1_ind + PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,panels[1].rC)
     u2_ind = u2_ind + PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,panels[2].rC)
     u3_ind = u3_ind + PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,panels[3].rC)
     u4_ind = u4_ind + PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,panels[4].rC)
end

induced = [dot(u1_ind,panels[1].n_hat);
           dot(u2_ind,panels[2].n_hat);
           dot(u3_ind,panels[3].n_hat);
           dot(u4_ind,panels[4].n_hat)]
fs = [dot(u_inf,panels[1].n_hat);
           dot(u_inf,panels[2].n_hat);
           dot(u_inf,panels[3].n_hat);
           dot(u_inf,panels[4].n_hat)]
diff = [abs(induced[1]) - abs(fs[1]);
        abs(induced[2]) - abs(fs[2]);
        abs(induced[3]) - abs(fs[3]);
        abs(induced[4]) - abs(fs[4])]

println("Difference between Induced and Free-Stream (PM2D): ",diff)



# Compare Geometries
println(" ===================== Geometries ===================== ")
for i=1:4
     println("Panel ",i)
     println(" r0:     HC = ",[x0_panels[i] y0_panels[i]],"  PM2D = ",panels[i].r0)
     println(" rC:     HC = ",[xC_panels[i] yC_panels[i]],"  PM2D = ",panels[i].rC)
     println(" n_hat:  HC = ",[n_hats[1,i],n_hats[2,i]],"  PM2D = ",panels[i].n_hat)
     println(" t_hat:  HC = ",[t_hats[1,i],t_hats[2,i]],"  PM2D = ",panels[i].t_hat)
end

if visualize
     # Visualize results
     org_points = vcat([[ [side[i,1], side[i,2], 0] for i in 1:size(side,1)] for side in [upper, lower]]...)
     line = org_points
     lines = [line]
     points, vtk_lines, vtk_values = vtk.lines2vtk(lines)
     vtk.generateVTK("airfoil", points; lines=vtk_lines)

     function wrap_u(X)
         u = [u_inf[1]; u_inf[2]; 0.0]
         for j=1:N_panels
             u = u + vcat(PM2D.calc_velocity("vortex",strengths[j],panels[j].r0,X),0)
         end
         return u
     end

     P_min = [-2.5,-2.5]
     P_max = [3.5,2.5]
     NDIVS = [100,50]
     fgrid = vtk.Grid(P_min, P_max, NDIVS)

     vtk.calculate_field(fgrid, wrap_u, "U", "vector")
     vtk.save(fgrid, "mygrid")

     run(`paraview --data="mygrid.vtk;airfoil.vtk"`)
end
