#yeah
import sympy as smp
import matplotlib.pyplot as plt
import numpy as np
import Copy_Variables_for_crossection as var

#Moment Distribution
#loading case 16 and loading case 12
#jump in the moment diagram 11.69m

a = smp.symbols('a', real=True)
h = smp.symbols('h', real=True)
d = smp.symbols('d', real=True, positive=True)
t = smp.symbols('t', real=True)
x = smp.symbols('x', real=True)
y = smp.symbols('y', real=True)
r = smp.symbols('r', real=True, positive=True)
alpha = smp.symbols('\u03B1', real=True)
beta = smp.symbols('\u03B2', real=True)
rho = smp.symbols('\u03C1', real=True)
sigma = smp.symbols('\u03C3', real=True)
theta = smp.symbols('\u03B8', real=True)
eta = smp.symbols('\u03B7', real=True)
xi = smp.symbols('\u03BE', real=True)
ixx = smp.symbols('Ixx', real=True)
iyy = smp.symbols('Ixy', real=True)
ixy = smp.symbols('Iyy', real=True)
mx = smp.symbols('Mx', real=True)
my = smp.symbols('My', real=True)

M1_12 = [-0.19569866609295541, 6.6801804985865445, 347.650403262873, -3851.0392032213094, -529912.4293065271, 13329458.866177564] 
M2_12 = [-0.19569866608662626, 6.680180498386035, 347.65040326524246, -3851.039203234463, -684065.9293064886, 14011505.048648307]
M1_16 = [0.40016072099856226, -9.397942238957981, -1047.384217723905, 6708.870613191617, 1582191.8673072485, -33488163.089606084] 
M2_16 = [0.40016072098555516, -9.397942238637903, -1047.3842177269048, 6708.870613212022, 1967575.6173070923, -39113307.35963459] 



#Normal force diagram (until 11.69)
#force is -500000 N (compression)
md_fun1_12 = smp.nsimplify(round(M1_12[0], 6))*x**5+smp.nsimplify(round(M1_12[1], 6))*x**4+smp.nsimplify(round(M1_12[2], 6))*x**3+smp.nsimplify(round(M1_12[3], 6))*x**2+smp.nsimplify(round(M1_12[4], 6))*x+smp.nsimplify(round(M1_12[5], 6))
md_fun2_12 = smp.nsimplify(round(M2_12[0], 6))*x**5+smp.nsimplify(round(M2_12[1], 6))*x**4+smp.nsimplify(round(M2_12[2], 6))*x**3+smp.nsimplify(round(M2_12[3], 6))*x**2+smp.nsimplify(round(M2_12[4], 6))*x+smp.nsimplify(round(M2_12[5], 6))

md_fun1_16 = smp.nsimplify(round(M1_16[0], 6))*x**5+smp.nsimplify(round(M1_16[1], 6))*x**4+smp.nsimplify(round(M1_16[2], 6))*x**3+smp.nsimplify(round(M1_16[3], 6))*x**2+smp.nsimplify(round(M1_16[4], 6))*x+smp.nsimplify(round(M1_16[5], 6))
md_fun2_16 = smp.nsimplify(round(M2_16[0], 6))*x**5+smp.nsimplify(round(M2_16[1], 6))*x**4+smp.nsimplify(round(M2_16[2], 6))*x**3+smp.nsimplify(round(M2_16[3], 6))*x**2+smp.nsimplify(round(M2_16[4], 6))*x+smp.nsimplify(round(M2_16[5], 6))

heaviside = smp.Heaviside(x-11.69, 0)
mx_fun_12 = md_fun1_12-md_fun1_12*heaviside+md_fun2_12*heaviside
mx_fun_16 = md_fun1_16-md_fun1_16*heaviside+md_fun2_16*heaviside
Ny_fun = -500000 + 500000*heaviside
mz_fun = -6008660+514000*x-(-6008660+514000*x)*heaviside

span = np.linspace(0, 0.5*var.Span, 1000, endpoint=True)
plt.figure()
plt.plot(span, -smp.lambdify([x], mx_fun_12)(span[0:]))
plt.plot(span, -smp.lambdify([x], mx_fun_16)(span[0:]))
plt.plot(span, smp.lambdify([x], Ny_fun)(span[0:]))
plt.plot(span, smp.lambdify([x], mz_fun)(span[0:]))
plt.show()

def sigma_z(mx, my, ixx, iyy, ixy, x, y):
  exp = ((mx*iyy-my*ixy)*y+(my*ixx-mx*ixy)*x)/(ixx*iyy-ixy**2)
  return exp
