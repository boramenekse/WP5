#also check compressive failure for this component
import sympy as smp
import Copy_Variables_for_crossection as var

sigma_cr = smp.symbols('\u03C3_{cr}', real=True)
kc = smp.symbols('k_c', real=True)
emod = smp.symbols('E', real=True, positive=True)
gmod = smp.symbols('G', real=True, positive=True)
v = smp.symbols('v', real=True, positive=True)
t = smp.symbols('t', real=True, positive=True)
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

# emod = var.e_mod
# gmod = var.g_mod
# b_2 = 0.5*var.Span
# cr = var.Chord_root
# g_e_relation = smp.Eq(gmod, smp.Rational(1, 2)*(emod/(1+v))) 
# v = smp.solve(g_e_relation, v)[0]
# taper = var.Taper_ratio

def sigma(kc, t, b):
  expr = smp.Rational(1, 12) * ((smp.pi**2 * kc * emod)/(1 - v**2)) * (t/b)**2
  return expr

def b(y):
  b_y = cr - cr*(1-taper)*(y/(b_2))
  return b_y