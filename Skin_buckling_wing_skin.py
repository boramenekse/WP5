#also check compressive failure for this component
import sympy as smp
import Copy_Variables_for_crossection as var

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
b_2 = 0.5*var.Span
cr = var.Chord_root
g_e_relation = smp.Eq(gmod, smp.Rational(1, 2)*(emod/(1+v))) 
v = smp.solve(g_e_relation, v)[0]
taper = var.Taper_ratio
ttopr = var.Sheet_top_th_root
tbotr = var.Sheet_bottom_th_root
ttopfun = ttopr*(1 + (var.Taper_ratio-1)*(y/b_2))
ttopfun = tbotr*(1 + (var.Taper_ratio-1)*(y/b_2))
theta_top = var.Sheet_top_angle
theta_bot = var.Sheet_bottom_angle

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
