import PM2D

upper,lower = PM2D.NACA_airfoil(2412)
PM2D.calc_coeff(upper,lower)
