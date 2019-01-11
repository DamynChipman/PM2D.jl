# ===== PM2D Verifications: Pressure Distribution =====

# == Imports ==
using PM2D

# == Set up geometry and grids ==
N = 200                             # Number of points
NPAN = N - 1                        # Number of panels
del_theta = (2*pi)/(N-1)            # Spacing in theta
R = 1.0                             # Radius
X = zeros(N)                        # X-coordinates
Y = zeros(N)                        # Y-coordinates
theta = zeros(NPAN-1)               # Angle from positive X-axis
CP_analytical = zeros(NPAN-1)       # Analytical CP
CP_PM2D = zeros(NPAN-1)             # PM2D (Semi-Analytical) CP
CP_SVPM = zeros(NPAN-1)             # SVPM CP
for i=1:N
    X[i] = R*cos((i-1)*del_theta)
    Y[i] = R*sin((i-1)*del_theta)
    if i != N-1
        theta[i] = (i-1)*del_theta
    end
end

# == CP: Analytical ==
for i=1:N-1
    CP_analytical[i] = 1 - 4*sin(theta[i])^2
end

# == CP: PM2D w/ Summation of Linear Vortex Lines (Semi-Analytical) ==
oper_cond = [10.0,0.0]
outputs = PM2D.PanelMethod(X,Y,oper_cond,["Vortex Sheet","Panels"])
panels = outputs[2]
alpha_coefs = outputs[1]

# Wrapper function for
function u_gamma(X)
    U = PM2D.CalcVortexSheetVelocity(panels,alpha_coefs,X,sigma=0.2,a=18.5*pi)
    return U
end

for i=1:NPAN-1
    CP_PM2D[i] = 1 - (norm(u_gamma([X[i],Y[i],0.0]))/(oper_cond[1]))^2
end
