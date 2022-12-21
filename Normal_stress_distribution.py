#yeah
import sympy as smp
import matplotlib.pyplot as plt
import numpy as np
import Variables_for_crossection_geometry as var
from Moment_of_inertia_z_span import return_parameters
from Moment_of_inertia_x_span import print_fit
import Centroid_x_direction as cgx
import Centroid_z_direction as cgz
from Normal_Stress_2 import print_results

ph_list = print_results()

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
izz_fun1 = par_list_z[0][0]*y**3+par_list_z[0][1]*y**2+par_list_z[0][2]*y+par_list_z[0][3]
izz_fun2 = par_list_z[1][0]*y**3+par_list_z[1][1]*y**2+par_list_z[1][2]*y+par_list_z[1][3]
izz_fun3 = par_list_z[2][0]*y**3+par_list_z[2][1]*y**2+par_list_z[2][2]*y+par_list_z[2][3]

# Ixx
par_list_x = print_fit()
ixx_fun1 = par_list_x[0][0]*y**3+par_list_x[0][1]*y**2+par_list_x[0][2]*y+par_list_x[0][3]
ixx_fun2 = par_list_x[1][0]*y**3+par_list_x[1][1]*y**2+par_list_x[1][2]*y+par_list_x[1][3]
ixx_fun3 = par_list_x[2][0]*y**3+par_list_x[2][1]*y**2+par_list_x[2][2]*y+par_list_x[2][3]

#Normal force diagram (until 11.69)
#force is -500000 N (compression)
md_fun1_12 = M1_12[0]*y**5+M1_12[1]*y**4+M1_12[2]*y**3+M1_12[3]*y**2+M1_12[4]*y+M1_12[5]
md_fun2_12 = M2_12[0]*y**5+M2_12[1]*y**4+M2_12[2]*y**3+M2_12[3]*y**2+M2_12[4]*y+M2_12[5]

md_fun1_16 = M1_16[0]*y**5+M1_16[1]*y**4+M1_16[2]*y**3+M1_16[3]*y**2+M1_16[4]*y+M1_16[5]
md_fun2_16 = M2_16[0]*y**5+M2_16[1]*y**4+M2_16[2]*y**3+M2_16[3]*y**2+M2_16[4]*y+M2_16[5]

heaviside = smp.Heaviside(y-11.69, 1)
mx_fun_12 = md_fun1_12-md_fun1_12*heaviside+md_fun2_12*heaviside
mx_fun_16 = md_fun1_16-md_fun1_16*heaviside+md_fun2_16*heaviside
Ny_fun = -500000 + 500000*heaviside
mz_fun = (-6008660+514000*y-(-6008660+514000*y)*heaviside)

span = np.linspace(0, 0.5*var.Span, 1000, endpoint=True)
# plt.figure()
# plt.plot(span, -smp.lambdify([y], mx_fun_12)(span[0:]))
# plt.plot(span, -smp.lambdify([y], mx_fun_16)(span[0:]))
# # plt.plot(span, smp.lambdify([y], Ny_fun)(span[0:]))
# plt.plot(span, smp.lambdify([y], mz_fun)(span[0:]))
# plt.show()



def sigma_z(mx, mz, ixx, izz):
  exp = ((mx*z)/ixx+(mz*x)/izz)
  return exp

# # alpha and beta are proportionality constants. Since we don't have a function for the centroid location as a funtion of y, we may move on by assuming the ratios will stay the same throughout the entire half-span
# alpha = cgx.Centroid_x/(0.55*var.Chord_root)
# beta = cgz.Centroid_z/var.Spar_fr_len_root
# cgx_fun = alpha*0.55*c_fun
# cgz_fun = beta*spar_fr_len_fun

stress_12_1 = sigma_z(mx_fun_12, mz_fun, ixx_fun1, izz_fun1)
stress_12_2 = sigma_z(mx_fun_12, mz_fun, ixx_fun2, izz_fun2)
stress_12_3 = sigma_z(mx_fun_12, mz_fun, ixx_fun3, izz_fun3)
stress_16_1 = sigma_z(mx_fun_16, mz_fun, ixx_fun1, izz_fun1)
stress_16_2 = sigma_z(mx_fun_16, mz_fun, ixx_fun2, izz_fun2)
stress_16_3 = sigma_z(mx_fun_16, mz_fun, ixx_fun3, izz_fun3)
# WE are calculating the top left corner location, which coordinates are (-cgx, -cgz)
# stress_top_corner_left = stress_12.subs([(z, -cgz.Centroid_z), (x, -cgx.Centroid_x), (y, 0)])
# # This time for the top right corner, whose coordinates are (sheet_top_length*cos(Sheet_top_Angle), cgz - sheet_top_lenght*sin(sheet_top_angle))
# stress_top_corner_right = stress_12.subs([(z, -cgz.Centroid_z+var.Sheet_top_len*smp.sin(var.Sheet_top_angle)), (x, var.Sheet_top_len*smp.cos(var.Sheet_top_angle)-cgx.Centroid_x), (y, 0)])

def curve_fit(a, b, c, d, x):
  exp = a*x**3+b*x**2+c*x+d
  return exp

fit1x = curve_fit(float(ph_list[0][0]), float(ph_list[0][1]), float(ph_list[0][2]), float(ph_list[0][3]), y)
fit1z = curve_fit(float(ph_list[1][0]), float(ph_list[1][1]), float(ph_list[1][2]), float(ph_list[1][3]), y)
fit2x = curve_fit(float(ph_list[2][0]), float(ph_list[2][1]), float(ph_list[2][2]), float(ph_list[2][3]), y)
fit2z = curve_fit(float(ph_list[3][0]), float(ph_list[3][1]), float(ph_list[3][2]), float(ph_list[3][3]), y)
fit3x = curve_fit(float(ph_list[4][0]), float(ph_list[4][1]), float(ph_list[4][2]), float(ph_list[4][3]), y)
fit3z = curve_fit(float(ph_list[5][0]), float(ph_list[5][1]), float(ph_list[5][2]), float(ph_list[5][3]), y)

c_fun = var.Chord_root*(1 + (var.Taper_ratio-1)*(y/(0.5*var.Span)))
spar_re_len_fun = var.Spar_re_len_root*(1+(var.Taper_ratio-1)*(y/(0.5*var.Span)))
br_z1 = (0.55*c_fun*smp.tan(var.Sheet_top_angle))+spar_re_len_fun-fit1z
br_x1 = 0.55*c_fun-fit1x
br_z2 = (0.55*c_fun*smp.tan(var.Sheet_top_angle))+spar_re_len_fun-fit2z
br_x2 = 0.55*c_fun-fit2x
br_z3 = (0.55*c_fun*smp.tan(var.Sheet_top_angle))+spar_re_len_fun-fit3z
br_x3 = 0.55*c_fun-fit3x

# print(var.Spar_fr_len)
# print(var.Spar_fr_len_root)
# print(var.Spar_fr_len_tip)
# print(var.Chord_root*var.Taper_ratio)
# print(cgx.Centroid_x, cgz.Centroid_z) # These values should be for the root
# print(stress_12.subs(y, 0))
# print(stress_12.subs(y, 0.5*var.Span))
stress_fun12_1 = smp.lambdify([y], stress_12_1.subs([(z, -fit1z), (x, -fit1z)]).simplify())
stress_fun12_2 = smp.lambdify([y], stress_12_2.subs([(z, -fit2z), (x, -fit2z)]).simplify())
stress_fun12_3 = smp.lambdify([y], stress_12_3.subs([(z, -fit3z), (x, -fit3x)]).simplify())
stress_fun16_1 = smp.lambdify([y], stress_16_1.subs([(z, br_z1), (x, br_x1)]).simplify())
stress_fun16_2 = smp.lambdify([y], stress_16_2.subs([(z, br_z2), (x, br_x2)]).simplify())
stress_fun16_3 = smp.lambdify([y], stress_16_3.subs([(z, br_z3), (x, br_x3)]).simplify())
plt.figure()
plt.title('Stress distribution at loading case 12 for philosophy 1')
plt.plot(span, stress_fun12_1(span[0:]))
plt.figure()
plt.title('Stress distribution at loading case 12 for philosophy 2')
plt.plot(span, stress_fun12_2(span[0:]))
plt.figure()
plt.title('Stress distribution at loading case 12 for philosophy 3')
plt.plot(span, stress_fun12_3(span[0:]))
plt.figure()
plt.title('Stress distribution at loading case 16 for philosophy 1')
plt.plot(span, stress_fun16_1(span[0:]))
plt.figure()
plt.title('Stress distribution at loading case 16 for philosophy 2')
plt.plot(span, stress_fun16_2(span[0:]))
plt.figure()
plt.title('Stress distribution at loading case 16 for philosophy 3')
plt.plot(span, stress_fun16_3(span[0:]))
plt.show()
# plt.figure()
# plt.plot(span, smp.lambdify([y], fit1x)(span[0:]))
# plt.show()
# plt.figure()
# plt.plot(span, smp.lambdify([y], fit1z)(span[0:]))
# plt.show()
# plt.figure()
# plt.plot(span, smp.lambdify([y], fit2x)(span[0:]))
# plt.show()
# plt.figure()
# plt.plot(span, smp.lambdify([y], fit2z)(span[0:]))
# plt.show()
# print(stress_top_corner_left/1e6)
# print(stress_top_corner_right/1e6)
# print(abs(stress_top_corner_right) > abs(stress_top_corner_left))