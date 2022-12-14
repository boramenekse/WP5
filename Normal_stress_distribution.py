#yeah
import sympy as smp
import matplotlib.pyplot as plt
import numpy as np
import Copy_Variables_for_crossection as var
from Moment_of_inertia_z_span import return_parameters
from Moment_of_inertia_x_span import print_fit
import Centroid_x_direction as cgx
import Centroid_z_direction as cgz

#Moment Distribution
#loading case 16 and loading case 12
#jump in the moment diagram 11.69m

a = smp.symbols('a', real=True)
h = smp.symbols('h', real=True)
d = smp.symbols('d', real=True, positive=True)
t = smp.symbols('t', real=True)
x = smp.symbols('x', real=True)
y = smp.symbols('y', real=True)
z = smp.symbols('z', real=True)
r = smp.symbols('r', real=True, positive=True)
alpha = smp.symbols('\u03B1', real=True)
beta = smp.symbols('\u03B2', real=True)
rho = smp.symbols('\u03C1', real=True)
sigma = smp.symbols('\u03C3', real=True)
theta = smp.symbols('\u03B8', real=True)
eta = smp.symbols('\u03B7', real=True)
xi = smp.symbols('\u03BE', real=True)
ixx = smp.symbols('Ixx', real=True)
iyy = smp.symbols('Iyy', real=True)
ixy = smp.symbols('Ixy', real=True)
mx = smp.symbols('Mx', real=True)
my = smp.symbols('My', real=True)

M1_12 = [-0.19569866609295541, 6.6801804985865445, 347.650403262873, -3851.0392032213094, -529912.4293065271, 13329458.866177564] 
M2_12 = [-0.19569866608662626, 6.680180498386035, 347.65040326524246, -3851.039203234463, -684065.9293064886, 14011505.048648307]
M1_16 = [0.40016072099856226, -9.397942238957981, -1047.384217723905, 6708.870613191617, 1582191.8673072485, -33488163.089606084] 
M2_16 = [0.40016072098555516, -9.397942238637903, -1047.3842177269048, 6708.870613212022, 1967575.6173070923, -39113307.35963459] 

# Izz
par_list_z = return_parameters()
izz_fun = smp.nsimplify(round(par_list_z[0], 6))*y**3+smp.nsimplify(round(par_list_z[1], 6))*y**2+smp.nsimplify(round(par_list_z[2], 6))*y+smp.nsimplify(round(par_list_z[3], 6))

# Ixx
par_list_x = print_fit()
ixx_fun = smp.nsimplify(round(par_list_x[0], 6))*y**3+smp.nsimplify(round(par_list_x[1], 6))*y**2+smp.nsimplify(round(par_list_x[2], 6))*y+smp.nsimplify(round(par_list_x[3], 6))

#Normal force diagram (until 11.69)
#force is -500000 N (compression)
md_fun1_12 = smp.nsimplify(round(M1_12[0], 6))*y**5+smp.nsimplify(round(M1_12[1], 6))*y**4+smp.nsimplify(round(M1_12[2], 6))*y**3+smp.nsimplify(round(M1_12[3], 6))*y**2+smp.nsimplify(round(M1_12[4], 6))*y+smp.nsimplify(round(M1_12[5], 6))
md_fun2_12 = smp.nsimplify(round(M2_12[0], 6))*y**5+smp.nsimplify(round(M2_12[1], 6))*y**4+smp.nsimplify(round(M2_12[2], 6))*y**3+smp.nsimplify(round(M2_12[3], 6))*y**2+smp.nsimplify(round(M2_12[4], 6))*y+smp.nsimplify(round(M2_12[5], 6))

md_fun1_16 = smp.nsimplify(round(M1_16[0], 6))*y**5+smp.nsimplify(round(M1_16[1], 6))*y**4+smp.nsimplify(round(M1_16[2], 6))*y**3+smp.nsimplify(round(M1_16[3], 6))*y**2+smp.nsimplify(round(M1_16[4], 6))*y+smp.nsimplify(round(M1_16[5], 6))
md_fun2_16 = smp.nsimplify(round(M2_16[0], 6))*y**5+smp.nsimplify(round(M2_16[1], 6))*y**4+smp.nsimplify(round(M2_16[2], 6))*y**3+smp.nsimplify(round(M2_16[3], 6))*y**2+smp.nsimplify(round(M2_16[4], 6))*y+smp.nsimplify(round(M2_16[5], 6))

heaviside = smp.Heaviside(y-0.35*0.5*var.Span, 0)
mx_fun_12 = md_fun1_12*(1-heaviside)+md_fun2_12*heaviside
mx_fun_16 = md_fun1_16*(1-heaviside)+md_fun2_16*heaviside
Ny_fun = -500000*(1-heaviside)
mz_fun = (-6008660+514000*y)*(1-heaviside)

span = np.linspace(0, 0.5*var.Span, 1000, endpoint=True)
#plt.figure()
#plt.plot(span, -smp.lambdify([x], mx_fun_12)(span[0:]))
#plt.plot(span, -smp.lambdify([x], mx_fun_16)(span[0:]))
#plt.plot(span, smp.lambdify([x], Ny_fun)(span[0:]))
#plt.plot(span, smp.lambdify([x], mz_fun)(span[0:]))
#plt.show()



def sigma_z(mx, mz, ixx, izz):
  exp = ((mx*z)/ixx)+((mz*x)/izz)
  return exp

c_fun = var.Chord_root*(1 + (var.Taper_ratio-1)*(y/(0.5*var.Span)))
spar_fr_len_fun = var.Spar_fr_len_root*(1+(var.Taper_ratio-1)*(y/(0.5*var.Span)))
# alpha and beta are proportionality constants. Since we don't have a function for the centroid location as a funtion of y, we may move on by assuming the ratios will stay the same throughout the entire half-span
alpha = cgx.Centroid_x/(0.55*var.Chord_root)
beta = cgz.Centroid_z/var.Spar_fr_len_root
cgx_fun = alpha*0.55*c_fun
cgz_fun = beta*spar_fr_len_fun

stress = sigma_z(mx_fun_12, mz_fun, ixx_fun, izz_fun)
# WE are calculating the top left corner location, which coordinates are (-cgx, -cgz)
stress_top_corner_left = stress.subs([(z, -cgz.Centroid_z), (x, -cgx.Centroid_x), (y, 0)])
# This time for the top right corner, whose coordinates are (sheet_top_length*cos(Sheet_top_Angle), cgz - sheet_top_lenght*sin(sheet_top_angle))
stress_top_corner_right = stress.subs([(z, -cgz.Centroid_z+var.Sheet_top_len*smp.sin(var.Sheet_top_angle)), (x, var.Sheet_top_len*smp.cos(var.Sheet_top_angle)-cgx.Centroid_x), (y, 0)])

print(var.Spar_fr_len)
print(var.Spar_fr_len_root)
print(var.Spar_fr_len_tip)
print(var.Chord_root*var.Taper_ratio)
print(cgx.Centroid_x, cgz.Centroid_z) # These values should be for the root
print(stress.subs(y, 0))
print(stress.subs(y, 0.5*var.Span))
print(stress.subs([(z, -cgz_fun), (x, -cgx_fun)]).simplify())
stress_fun = smp.lambdify([y], stress.subs([(z, -cgz_fun), (x, -cgx_fun)]).simplify())
plt.figure()
plt.plot(span, stress_fun(span[0:]))
plt.show()
print(stress_top_corner_left/1e6)
print(stress_top_corner_right/1e6)
print(abs(stress_top_corner_right) > abs(stress_top_corner_left))