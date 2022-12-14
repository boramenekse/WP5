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

heaviside = smp.Heaviside(y-11.69, 1)
nd = (514000*y - (514000*y - 6008660)*heaviside - 6008660)*(2.68003769312149e-19*y**3 - 1.6822756916907e-17*y**2 - 0.0704116492628627*y + 3.26464295018786)/(-8.95576123019816e-5*y**3 + 0.00863008340506688*y**2 - 0.291525157620533*y + 3.45400374877192) + (7.87266384897094e-20*y**3 - 3.47925199872395e-18*y**2 - 0.00511997686115935*y + 0.23738822396431)*(0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y + (0.400160720985555*y**5 - 9.3979422386379*y**4 - 1047.3842177269*y**3 + 6708.87061321202*y**2 + 1967575.61730709*y - 39113307.3596346)*heaviside - (0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y - 33488163.0896061)*heaviside - 33488163.0896061)/(-3.69724112051317e-6*y**3 + 0.000356279030822228*y**2 - 0.0120351444942484*y + 0.142592956886914)

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
tbotfun = tbotr*(1 + (taper-1)*(y/b_2))
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

# KS for C with SS loaded edges 
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
# plt.figure()
# plt.ylim((0, 16))
# plt.xlim((0, 5))
# plt.plot(a_b, smp.lambdify([ab], kfun)(a_b[0:]))
# plt.show()

sfc = 1.5

import math
import numpy as np


def zbar(p, no_str, fr_t_root, re_t_root, top_sheet_t_root, bottom_sheet_t_root):
    #input parameters

    #Taper ratio and span
    Taper_ratio = 0.279881175
    Span = 66.77675839 #[m]
    Chord_root = 11.59426354 #[m]

    #material properties for AL6061-T6
    e_mod = 69*10**9 #[pa]
    g_mod = 26*10**9 #[pa]

    #stringers
    Str_A = 0.001 #[m^2]
    Str_N = no_str #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...

    #spanwise location
    y = p

    #GEOMETRIES AT THE ROOT

    #spars root (and tip for front length)
    Spar_fr_len_root = (0.045100 + 0.045200)*Chord_root
    Spar_fr_len_tip = Spar_fr_len_root * Taper_ratio
    Spar_fr_th_root= fr_t_root

    Spar_re_len_root = (0.033900 + 0.016200)*Chord_root
    Spar_re_th_root = re_t_root

    #sheets root 
    Sheet_top_len_root = math.sqrt((0.045100 - 0.033900)**2 + 0.55**2)*Chord_root
    Sheet_top_th_root = top_sheet_t_root

    Sheet_bottom_len_root = math.sqrt((0.045200 - 0.016200)**2 + 0.55**2)*Chord_root
    Sheet_bottom_th_root = bottom_sheet_t_root

    Sheet_top_angle = math.atan((0.045100-0.033900)/(0.75-0.2))
    Sheet_bottom_angle = math.atan((0.045200-0.016200)/(0.75-0.2))

    #Ratio's
    Constant = (Spar_fr_len_tip - Spar_fr_len_root)/(Span/2)

    Ratio_Sheet_top_len = Spar_fr_len_root / Sheet_top_len_root
    Ratio_Spar_re_len = Spar_fr_len_root / Spar_re_len_root
    Ratio_Sheet_bottom_len = Spar_fr_len_root / Sheet_bottom_len_root
    Ratio_Sheet_top_th = Spar_fr_len_root / Sheet_top_th_root
    Ratio_Sheet_bottom_th = Spar_fr_len_root / Sheet_bottom_th_root
    Ratio_Spar_re_th = Spar_fr_len_root / Spar_re_th_root
    Ratio_Spar_fr_th = Spar_fr_len_root / Spar_fr_th_root



    #VALUES FOR GEOMETRIES THROUGHOUT THE SPAN

    #spars along span
    Spar_fr_len = Spar_fr_len_root + Constant * y 
    Spar_fr_th = Spar_fr_len / Ratio_Spar_fr_th

    Spar_re_th = Spar_fr_len / Ratio_Spar_re_th
    Spar_re_len = Spar_fr_len / Ratio_Spar_re_len

    #sheets along span
    Sheet_top_th = Spar_fr_len / Ratio_Sheet_top_th
    Sheet_top_len = Spar_fr_len / Ratio_Sheet_top_len

    Sheet_bottom_th = Spar_fr_len / Ratio_Sheet_bottom_th
    Sheet_bottom_len = Spar_fr_len / Ratio_Sheet_bottom_len

    #Check:
    #print(Spar_fr_len)
    #print(Spar_fr_th)
    #print(Spar_re_len)
    #print(Spar_re_th)
    #print(Sheet_top_len)
    #print(Sheet_top_th)
    #print(Sheet_bottom_len)
    #print(Sheet_bottom_th)

    #Total area and distance to stringer from start of sheet
    total_area = Spar_fr_len * Spar_fr_th + Spar_re_len * Spar_re_th + Sheet_bottom_len * Sheet_bottom_th + Sheet_top_len * Sheet_top_th + Str_A * Str_N
    Str_dis_top = Sheet_top_len / ((Str_N / 2) + 1)
    Str_dis_bottom = Sheet_bottom_len / ((Str_N / 2) + 1)

    #Z tilda calculations for every part
    z_tilda_Spar_fr = Spar_fr_len / 2
    z_tilda_Spar_re = (Sheet_top_len) * math.sin(Sheet_top_angle) + Spar_re_len / 2
    z_tilda_Sheet_top = (Sheet_top_len / 2) * math.sin(Sheet_top_angle)
    z_tilda_Sheet_bottom = (Spar_fr_len) - ((Sheet_bottom_len / 2) * math.sin(Sheet_bottom_angle))

    z_tilda_Str_bottom = []
    for i in range(int(Str_N / 2)):
        z_tilda_Str_N = Spar_fr_len - (Str_dis_bottom * (i+1) * math.sin(Sheet_bottom_angle))
        z_tilda_Str_bottom.append(z_tilda_Str_N)

    z_tilda_Str_top = []
    for i in range(int(Str_N / 2)):
        z_tilda_Str_P = Str_dis_top * (i+1) * math.sin(Sheet_top_angle)
        z_tilda_Str_top.append(z_tilda_Str_P)

    #Area's of seperate parts
    Spar_fr_A = Spar_fr_len*Spar_fr_th
    Spar_re_A = Spar_re_len*Spar_re_th
    Sheet_top_A = Sheet_top_len*Sheet_top_th
    Sheet_bottom_A = Sheet_bottom_len*Sheet_bottom_th

    #Area multiplied by z tilda for every part
    sum_of_products_Str_top = 0
    for i in range(int(Str_N / 2)):
        sum_of_products_Str_top += z_tilda_Str_top[i] * Str_A

    sum_of_products_Str_bottom = 0
    for i in range(int(Str_N / 2)):
        sum_of_products_Str_bottom += z_tilda_Str_bottom[i] * Str_A

    sum_of_products_Spar = z_tilda_Spar_fr * Spar_fr_A + z_tilda_Spar_re * Spar_re_A

    sum_of_products_Sheet = z_tilda_Sheet_top * Sheet_top_A + z_tilda_Sheet_bottom * Sheet_bottom_A

    sum_of_products = sum_of_products_Sheet + sum_of_products_Spar + sum_of_products_Str_bottom + sum_of_products_Str_top

    #Calculation of centroid
    Centroid_z = sum_of_products / total_area
    return Centroid_z

def Moi_z_wingbox(p, no_str, fr_t_root, re_t_root, top_sheet_t_root, bottom_sheet_t_root):
    #input parameters

    #Taper ratio and span
    Taper_ratio = 0.279881175
    Span = 66.77675839 #[m]
    Chord_root = 11.59426354 #[m]

    #material properties for AL6061-T6
    e_mod = 69*10**9 #[pa]
    g_mod = 26*10**9 #[pa]

    #stringers
    Str_A = 0.001 #[m^2]
    Str_N = no_str #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...

    #spanwise location
    y = p

    #GEOMETRIES AT THE ROOT

    #spars root (and tip for front length)
    Spar_fr_len_root = (0.045100 + 0.045200)*Chord_root
    Spar_fr_len_tip = Spar_fr_len_root * Taper_ratio
    Spar_fr_th_root= fr_t_root

    Spar_re_len_root = (0.033900 + 0.016200)*Chord_root
    Spar_re_th_root = re_t_root

    #sheets root 
    Sheet_top_len_root = math.sqrt((0.045100 - 0.033900)**2 + 0.55**2)*Chord_root
    Sheet_top_th_root = top_sheet_t_root

    Sheet_bottom_len_root = math.sqrt((0.045200 - 0.016200)**2 + 0.55**2)*Chord_root
    Sheet_bottom_th_root = bottom_sheet_t_root

    Sheet_top_angle = math.atan((0.045100-0.033900)/(0.75-0.2))
    Sheet_bottom_angle = math.atan((0.045200-0.016200)/(0.75-0.2))


    #Ratio's
    Constant = (Spar_fr_len_tip - Spar_fr_len_root)/(Span/2)

    Ratio_Sheet_top_len = Spar_fr_len_root / Sheet_top_len_root
    Ratio_Spar_re_len = Spar_fr_len_root / Spar_re_len_root
    Ratio_Sheet_bottom_len = Spar_fr_len_root / Sheet_bottom_len_root
    Ratio_Sheet_top_th = Spar_fr_len_root / Sheet_top_th_root
    Ratio_Sheet_bottom_th = Spar_fr_len_root / Sheet_bottom_th_root
    Ratio_Spar_re_th = Spar_fr_len_root / Spar_re_th_root
    Ratio_Spar_fr_th = Spar_fr_len_root / Spar_fr_th_root



    #VALUES FOR GEOMETRIES THROUGHOUT THE SPAN

    #spars along span
    Spar_fr_len = Spar_fr_len_root + Constant * y 
    Spar_fr_th = Spar_fr_len / Ratio_Spar_fr_th

    Spar_re_th = Spar_fr_len / Ratio_Spar_re_th
    Spar_re_len = Spar_fr_len / Ratio_Spar_re_len

    #sheets along span
    Sheet_top_th = Spar_fr_len / Ratio_Sheet_top_th
    Sheet_top_len = Spar_fr_len / Ratio_Sheet_top_len

    Sheet_bottom_th = Spar_fr_len / Ratio_Sheet_bottom_th
    Sheet_bottom_len = Spar_fr_len / Ratio_Sheet_bottom_len

    #Check:
    #print(Spar_fr_len)
    #print(Spar_fr_th)
    #print(Spar_re_len)
    #print(Spar_re_th)
    #print(Sheet_top_len)
    #print(Sheet_top_th)
    #print(Sheet_bottom_len)
    #print(Sheet_bottom_th)

    #Total area and distance to stringer from start of sheet
    total_area = Spar_fr_len*Spar_fr_th + Spar_re_len*Spar_re_th + Sheet_bottom_len*Sheet_bottom_th + Sheet_top_len*Sheet_top_th + Str_A*Str_N
    Str_dis_top = Sheet_top_len / (Str_N/2 + 1)
    Str_dis_bottom = Sheet_bottom_len / (Str_N/2 + 1)

    #Z tilda calculations for every part
    x_tilda_Spar_fr = 0
    x_tilda_Spar_re = Sheet_top_len * math.cos(Sheet_top_angle)
    x_tilda_Sheet_top = (Sheet_top_len/2) * math.cos(Sheet_top_angle)
    x_tilda_Sheet_bottom = (Sheet_bottom_len/2) * math.cos(Sheet_bottom_angle)

    x_tilda_Str_bottom = []
    for i in range(int(Str_N/2)):
        x_tilda_Str_bottom.append( (Str_dis_bottom * (i+1)) * math.cos(Sheet_bottom_angle))

    x_tilda_Str_top = []
    for i in range(int(Str_N/2)):
        x_tilda_Str_top.append( (Str_dis_top * (i+1)* math.cos(Sheet_top_angle)) )

    #Area's of seperate parts
    Spar_fr_A = Spar_fr_len*Spar_fr_th
    Spar_re_A = Spar_re_len*Spar_re_th
    Sheet_top_A = Sheet_top_len*Sheet_top_th
    Sheet_bottom_A = Sheet_bottom_len*Sheet_bottom_th

    #Area multiplied by z tilda for every part
    sum_of_products_Str = 0
    for i in range(int(Str_N/2)):
        sum_of_products_Str = sum_of_products_Str + Str_A * x_tilda_Str_bottom[i]
        sum_of_products_Str = sum_of_products_Str + Str_A * x_tilda_Str_top[i]

    sum_of_products = sum_of_products_Str + Spar_fr_A*x_tilda_Spar_fr + Spar_re_A*x_tilda_Spar_re + Sheet_top_A*x_tilda_Sheet_top + Sheet_bottom_A*x_tilda_Sheet_bottom

    #Calculation of centroid
    Centroid_x = sum_of_products/total_area
    return Centroid_x

no_list = [0, 18, 34]
fr_t_list = [0.045, 0.121, 0.100]
re_t_list = [0.045, 0.121, 0.100]
top_t_list = [0.063, 0.056, 0.057]
bottom_t_list = [0.063, 0.056, 0.057]

Span_y_z = np.arange(0.0, (var.Span / 2), 1.0)
index = 0

ph1x = []
ph2x = []
ph3x = []
ph1z = []
ph2z = []
ph3z = []
for i in Span_y_z:
  x = Moi_z_wingbox(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  z = zbar(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  ph1x.append(x)
  ph1z.append(z)
index+=1
for i in Span_y_z:
  x = Moi_z_wingbox(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  z = zbar(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  ph2x.append(x)
  ph2z.append(z)
index+=1
for i in Span_y_z:
  x = Moi_z_wingbox(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  z = zbar(i, no_list[index], fr_t_list[index], re_t_list[index], top_t_list[index], bottom_t_list[index])
  ph3x.append(x)
  ph3z.append(z)

fit1x = np.polyfit(Span_y_z, ph1x, 3)
fit1z = np.polyfit(Span_y_z, ph1z, 3)
fit2x = np.polyfit(Span_y_z, ph2x, 3)
fit2z = np.polyfit(Span_y_z, ph2z, 3)
fit3x = np.polyfit(Span_y_z, ph3x, 3)
fit3z = np.polyfit(Span_y_z, ph3z, 3)

def print_results():
  return fit1x, fit1z, fit2x, fit2z, fit3x, fit3z