#also check compressive failure for this component
import sympy as smp
import Copy_Variables_for_crossection as var
import numpy as np
import matplotlib.pyplot as plt

sigma_cr = smp.symbols('\u03C3_{cr}', real=True)
theta = smp.symbols('\u03B8', real=True)
kc = smp.symbols('k_c', real=True)
emod = smp.symbols('E', real=True, positive=True)
gmod = smp.symbols('G', real=True, positive=True)
v = smp.symbols('v', real=True, positive=True)
t = smp.symbols('t', real=True, positive=True)
ttopr = smp.symbols('t_{top_r}', real=True, positive=True)
tbotr = smp.symbols('t_{bottom_r}', real=True, positive=True)
b = smp.symbols('b', real=True, positive=True)
b_2 = smp.symbols('b/2', real=True, positive=True)
x = smp.symbols('x')
y = smp.symbols('y')
z = smp.symbols('z')
a = smp.symbols('a', real=True) 
h = smp.symbols('h', real=True)
d = smp.symbols('d', real=True, positive=True) 
cr = smp.symbols('c_r', real=True, positive=True) 
taper = smp.symbols('\u03BB', real=True, positive=True) 

emod = var.e_mod
gmod = var.g_mod
b_2 = smp.Rational(1, 2)*smp.nsimplify(round(var.Span, 7))
cr = smp.nsimplify(round(var.Chord_root, 7))
g_e_relation = smp.Eq(gmod, smp.Rational(1, 2)*(emod/(1+v))) 
v = smp.solve(g_e_relation, v)[0]
taper = smp.nsimplify(round(var.Taper_ratio, 7))
ttopr = smp.nsimplify(round(var.Sheet_top_th_root, 7))
tbotr = smp.nsimplify(round(var.Sheet_bottom_th_root, 7))
ttopfun = ttopr*(1 + (taper-1)*(y/b_2))
ttopfun = tbotr*(1 + (taper-1)*(y/b_2))
theta_top = smp.nsimplify(round(var.Sheet_top_angle, 7))
theta_bot = smp.nsimplify(round(var.Sheet_bottom_angle, 7))

def sigma(kc, t, b):
  expr = smp.Rational(1, 12) * ((smp.pi**2 * kc * emod)/(1 - v**2)) * (t/b)**2
  return expr.simplify()

def b(y, theta):
  c_local = cr - cr*(1-taper)*(y/(b_2))
  '''Wingbox starts at 20 percent of the local chord and end at 75 percent of it
     So b is 55 percent of the local chord length -> 11/20
     Sheets are at an angle
  '''
  b_y = c_local*smp.Rational(11, 20)/smp.cos(theta)
  return b_y.simplify()

a_b = np.linspace(0, 5, 1000, endpoint=True)
ab = smp.symbols('a_b', real=True, positive=True)
m = smp.symbols('m', real=True, positive=True)
def k(m, a_b):
  expr = (1/a_b)**2*(m+(1/m)*a_b**2)**2
  return expr

m1 = k(1, ab)
m2 = k(2, ab)
m3 = k(3, ab)
m4 = k(4, ab)
m5 = k(5, ab)
m6 = k(6, ab)
kfun = smp.Piecewise((m1, (ab>=0) & (ab <= float(smp.solve(m1-m2, ab)[0]))), (m2, (ab > float(smp.solve(m1-m2, ab)[0])) & (ab <= float(smp.solve(m2-m3, ab)[0]))), (m3, (ab > float(smp.solve(m2-m3, ab)[0])) & (ab <= float(smp.solve(m3-m4, ab)[0]))), (m4, (ab > float(smp.solve(m3-m4, ab)[0])) & (ab <= float(smp.solve(m4-m5, ab)[0]))), (m5, (ab > float(smp.solve(m4-m5, ab)[0])) & (ab <= 5)), (0, True))
plt.figure()
plt.ylim((0, 16))
plt.xlim((0, 5))
plt.plot(a_b, smp.lambdify([ab], kfun)(a_b[0:]))
plt.show()