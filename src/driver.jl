import PM2D

upper,lower = PM2D.NACA_airfoil(2412)
panels = PM2D.NACA_panels(upper,lower)

function u_vortex(r0,rC)
    u = PM2D.calc_velocity("vortex",1.0,r0,rC)
    return u
end

strengths = PM2D.calc_strengths(u_vortex,panels,[-5.0;0.0],true)
println(strengths)

# pan1 = PM2D.Panel2D("vortex",[0.0 0.0],[1.0 0.0])
# pan2 = PM2D.Panel2D("vortex",[1.0 0.0],[2.0 0.0])
# pan3 = PM2D.Panel2D("vortex",[2.0 0.0],[3.0 0.0])
# pan4 = PM2D.Panel2D("vortex",[3.0 0.0],[4.0 0.0])
# panels = [pan1 pan2 pan3 pan4]
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
