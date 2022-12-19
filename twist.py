import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt
import Variables_for_crossection_geometry as Var
from scipy.optimize import curve_fit
from scipy import integrate
import sympy as smp
import Moment_of_inertia_x_span as MOIX
import Moment_of_inertia_z_span as MOIZ


torsion_stiffness_z_span = []

#Torsian distribution constants for LC-8: 
#Part 1: [-2.3564272917644833, -110.84913656162982, -63250.84560607412, 3855696.551646198, -58140844.54168405] 
#Part 2:[27.128436138924133, -1751.353291452992, -26414.77577940539, 3572871.2685643807, -58379080.83130088


def t_stiffness_z_span(p):
    #input parameters

    #Taper ratio and span
    Taper_ratio = Var.Taper_ratio
    Span = Var.Span
    Chord_root = 11.59426354 #[m]

    #material properties for AL6061-T6
    e_mod = Var.e_mod #[gpa]
    g_mod = Var.g_mod #[gpa]

    #stringers
    Str_A = Var.Str_A
    Str_N = Var.Str_N #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...

    #spanwise location
    y = p

    #GEOMETRIES AT THE ROOT

    #spars root (and tip for front length)
    Spar_fr_len_root = Var.Spar_fr_len_root
    Spar_fr_len_tip = Var.Spar_fr_len_tip
    Spar_fr_th_root= Var.Spar_fr_th_root

    Spar_re_len_root = Var.Spar_re_len_root
    Spar_re_th_root = Var.Spar_re_th_root

    #sheets root 
    Sheet_top_len_root = Var.Sheet_top_len_root
    Sheet_top_th_root = Var.Sheet_top_th_root

    Sheet_bottom_len_root = Var.Sheet_bottom_len_root
    Sheet_bottom_th_root = Var.Sheet_bottom_th_root

    Sheet_top_angle = Var.Sheet_top_angle
    Sheet_bottom_angle = Var.Sheet_bottom_angle



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

    #Variables
    a= Spar_fr_len 
    b= Spar_re_len 
    h= Sheet_top_len
    l= Sheet_bottom_len
    t_1= Sheet_bottom_th 
    t_2= Spar_re_th 
    t_3= Sheet_top_th 
    t_4= Spar_fr_th 
    

    #A_m= a * 0.55 * Chord_root - (a - math.sin(Sheet_top_angle) * h - b) * 0.55 * Chord_root * 0.5 - (a - math.sin(Sheet_bottom_angle) * l - b) * 0.55 * Chord_root * 0.5
    A_m = ((a + b) / 2) * h * math.cos(Sheet_top_angle)
    dst_1=l/t_1
    dst_2=b/t_2
    dst_3=h/t_3
    dst_4=a/t_4

    #stiffness formula:
    GJ= g_mod*(4*A_m**2)/(dst_1+dst_2+dst_3+dst_4)
    #GJ = g_mod*(MOIX.Moi_x_wingbox(y) + MOIZ.Moi_z_wingbox(y))

    return GJ


#Calculating all the different values along the span for both x and y values
Span_y_z = np.linspace(0, 0.5*Var.Span, 1000, endpoint=True) #Span array for torsional stiffness (GJ) diagram



for i in Span_y_z:
    torsion_stiffness_z_span.append(t_stiffness_z_span(i))



#Set up for test function:
def test_function_GJ(x, A, B, C, D):
    y = A*(x**3) + B*(x**2) + C*x + D
    return y 

Parameters_1, coVariance = curve_fit(test_function_GJ, Span_y_z, torsion_stiffness_z_span)

Fit_A = Parameters_1[0]
Fit_B = Parameters_1[1]
Fit_C = Parameters_1[2]
Fit_D = Parameters_1[3]

Fit_y = test_function_GJ(Span_y_z, Fit_A, Fit_B, Fit_C, Fit_D)
#print('The values for A, B, C and D in the function Ax^3 + Bx^2 + Cx + D are:', Fit_A, Fit_B, Fit_C, Fit_D)

#Plotting the results (WP4.3a)
#plt.plot(Span_y_z, torsion_stiffness_z_span, 'o', label='Data')
#plt.plot(Span_y_z, Fit_y, '-', label='Fit')
#plt.xlabel("Spanwise location")
#plt.ylabel("Torsional Stiffness")
#plt.show()

Span_y_1 = np.linspace(0, 0.35*0.5*Var.Span, 1000, endpoint=True) # Span arrange for twist differential before discontinuity
Span_y_2 = np.linspace(0.35*0.5*Var.Span, 0.5*Var.Span, 1000, endpoint=True) # Span arrange for twist differential after discontinuity (0.35 * 0.5 * Var.Span)
#print(Span_y_1[-1], Span_y_2[0])
#Checking the torque distribution
def torque_b(p):
    #Torsian distribution constants for LC-8,12,16: Part 1: 

    #LC_1 = [-2.35642729183374, -110.84913656004754, -63250.84560608597, 3855696.5516462256, -57608532.93666329] #8 
    LC_1 = [0.05255743793122159, -13.820176577363013, -2166.4464360103175, 191852.8205224776, -297957.8252965813] #12
    #LC_1 = [-0.7447294322108541, 76.08911284422432, 4845.807207531824, -497278.2940762989, 6429799.316955603] #16

    #Toque digram before discontinuity 
    T1_1 = LC_1[0]
    T1_2 = LC_1[1] 
    T1_3 = LC_1[2] 
    T1_4 = LC_1[3]
    T1_5 = LC_1[4] 

    torque_b_z = (T1_1*p**4+T1_2*p**3+T1_3*p**2+T1_4*p+T1_5)
    return torque_b_z

def torque_a(p):
    #Torsian distribution constants for LC-8,12,16: Part 2:
    #LC_2 = [27.128436138924133, -1751.353291452992, -26414.77577940539, 3572871.2685643807, -58379080.83130088] #8
    LC_2 = [0.05255743793178574, -13.82017657737858, -2166.4464360101647, 191852.82052247645, -3541466.425956474] #12
    #LC_2 = [-0.7447294322224564, 76.08911284444716, 4845.8072075307255, -497278.2940763008, 9294745.002470266] #16

    #Toque digram after discontinuity
    T2_1 = LC_2[0]
    T2_2 = LC_2[1] 
    T2_3 = LC_2[2] 
    T2_4 = LC_2[3] 
    T2_5 = LC_2[4]

    torque_a_z = (T2_1*p**4+T2_2*p**3+T2_3*p**2+T2_4*p+T2_5)
    return torque_a_z

torque_b_d = []
torque_a_d = []

for i in Span_y_z:
    torque_b_d.append(torque_b(i))

for i in Span_y_z:
    torque_a_d.append(torque_a(i))

#plt.plot(Span_y_z, torque_b_d, 'o', label='Data Torque before')
#plt.plot(Span_y_z, torque_a_d, 'o', label='Data Torque after')
#plt.legend()
#plt.show()


#dthet/dz functions:
def dtheta_dz_z(p):
    #Torsian distribution constants for LC-8,12,16: Part 2:
    #LC_2 = [27.128436138924133, -1751.353291452992, -26414.77577940539, 3572871.2685643807, -58379080.83130088] #8
    LC_2 = [0.05255743793178574, -13.82017657737858, -2166.4464360101647, 191852.82052247645, -3541466.425956474] #12
    #LC_2 = [-0.7447294322224564, 76.08911284444716, 4845.8072075307255, -497278.2940763008, 9294745.002470266] #16

    #Toque digram after discontinuity
    T2_1 = LC_2[0]
    T2_2 = LC_2[1] 
    T2_3 = LC_2[2] 
    T2_4 = LC_2[3] 
    T2_5 = LC_2[4]

    dtheta_dz = (T2_1*p**4+T2_2*p**3+T2_3*p**2+T2_4*p+T2_5)/(test_function_GJ(p, Fit_A, Fit_B, Fit_C, Fit_D))

    return dtheta_dz


def dtheta_dz_1(p):

#Torsian distribution constants for LC-8,12,16: Part 1: 
    #LC_1 = [-2.35642729183374, -110.84913656004754, -63250.84560608597, 3855696.5516462256, -57608532.93666329] #8 
    LC_1 = [0.05255743793122159, -13.820176577363013, -2166.4464360103175, 191852.8205224776, -297957.8252965813] #12
    #LC_1 = [-0.7447294322108541, 76.08911284422432, 4845.807207531824, -497278.2940762989, 6429799.316955603] #16

    #Toque digram before discontinuity 
    T1_1 = LC_1[0]
    T1_2 = LC_1[1] 
    T1_3 = LC_1[2] 
    T1_4 = LC_1[3]
    T1_5 = LC_1[4]

    dtheta_dz = (T1_1*p**4+T1_2*p**3+T1_3*p**2+T1_4*p+T1_5)/(test_function_GJ(p, Fit_A, Fit_B, Fit_C, Fit_D))

    return dtheta_dz


def dtheta_dz_2(p):
    #Torsian distribution constants for LC-8,12,16: Part 2:
    #LC_2 = [27.128436138924133, -1751.353291452992, -26414.77577940539, 3572871.2685643807, -58379080.83130088] #8
    LC_2 = [0.05255743793178574, -13.82017657737858, -2166.4464360101647, 191852.82052247645, -3541466.425956474] #12
    #LC_2 = [-0.7447294322224564, 76.08911284444716, 4845.8072075307255, -497278.2940763008, 9294745.002470266] #16

    #Toque digram after discontinuity
    T2_1 = LC_2[0]
    T2_2 = LC_2[1] 
    T2_3 = LC_2[2] 
    T2_4 = LC_2[3] 
    T2_5 = LC_2[4]

    dtheta_dz = (T2_1*p**4+T2_2*p**3+T2_3*p**2+T2_4*p+T2_5)/(test_function_GJ(p, Fit_A, Fit_B, Fit_C, Fit_D))

    return dtheta_dz

int_dtheta_dz_z = []
int_dtheta_dz_1 = []
int_dtheta_dz_2 = []

for i in Span_y_z:
    int_dtheta_dz_z.append(dtheta_dz_z(i))    

for i in Span_y_1:
   int_dtheta_dz_1.append(dtheta_dz_1(i)) 

for i in Span_y_2:
   int_dtheta_dz_2.append(dtheta_dz_2(i)) 

#Set up test  function for differential equation z:
def test_function_d_dz_z(x, A, B, C, D, E, F, G, H):
    y = A*(x**7) + B*(x**6) + C*(x**5) + D*(x**4) + E*(x**3) + F*(x**2) + G*(x) + H
    return y
Parameters_z, covariance = curve_fit(test_function_d_dz_z, Span_y_z, int_dtheta_dz_z)

d_dz_z_A = Parameters_z[0]
d_dz_z_B = Parameters_z[1]
d_dz_z_C = Parameters_z[2]
d_dz_z_D = Parameters_z[3]
d_dz_z_E = Parameters_z[4]
d_dz_z_F = Parameters_z[5]
d_dz_z_G = Parameters_z[6]
d_dz_z_H = Parameters_z[7]


#Set up for test function for diferential equation 1:
def test_function_d_dz_1(x, A, B, C, D):
    y = A*(x**3) + B*(x**2) + C*x + D
    return y 

Parameters_2, coVariance = curve_fit(test_function_d_dz_1, Span_y_1, int_dtheta_dz_1)

d_dz_1_A = Parameters_2[0]
d_dz_1_B = Parameters_2[1]
d_dz_1_C = Parameters_2[2]
d_dz_1_D = Parameters_2[3]

#Set up for test function for diferential equation 2:
def test_function_d_dz_2(x, A, B, C, D, E, F):
    y = A*(x**5) + B*(x**4) + C*(x**3) + D*(x**2) + E*x + F
    return y 

Parameters_3, coVariance = curve_fit(test_function_d_dz_2, Span_y_2, int_dtheta_dz_2)

d_dz_2_A = Parameters_3[0]
d_dz_2_B = Parameters_3[1]
d_dz_2_C = Parameters_3[2]
d_dz_2_D = Parameters_3[3]
d_dz_2_E = Parameters_3[4]
d_dz_2_F = Parameters_3[5]

#print(d_dz_1_A, d_dz_1_B, d_dz_1_C, d_dz_1_D) #Check
#print(d_dz_2_A, d_dz_2_B, d_dz_2_C, d_dz_2_D, d_dz_2_E, d_dz_2_F) #Check

Fit_d_dz_z = test_function_d_dz_z(Span_y_z, d_dz_z_A, d_dz_z_B, d_dz_z_C, d_dz_z_D, d_dz_z_E, d_dz_z_F, d_dz_z_G, d_dz_z_H)
Fit_d_dz_1 = test_function_d_dz_1(Span_y_1, d_dz_1_A, d_dz_1_B, d_dz_1_C, d_dz_1_D)
Fit_d_dz_2 = test_function_d_dz_2(Span_y_2, d_dz_2_A, d_dz_2_B, d_dz_2_C, d_dz_2_D, d_dz_2_E, d_dz_2_F)

#plt.plot(Span_y_z, int_dtheta_dz_z, 'o', label='Data Z')
#plt.plot(Span_y_1, int_dtheta_dz_1, 'o', label='Data 1')
#plt.plot(Span_y_2, int_dtheta_dz_2, 'o', label='Data 2')
#plt.plot(Span_y_z, Fit_d_dz_z, '-', label='Fit Z')
#plt.plot(Span_y_1, Fit_d_dz_1, '-', label='Fit 1')
#plt.plot(Span_y_2, Fit_d_dz_2, '-', label='Fit 2')
#plt.legend()
#plt.xlabel("Spanwise location")
#plt.ylabel("dtheta/dz")
#plt.show()

#Integrating dtheta/dz (Bora's method)
diff_function_1 = np.polyfit(Span_y_1, int_dtheta_dz_1, 3) # sympy fit for dtheta_dz_1
coef_diff_func_1_x3= float(diff_function_1[0])
coef_diff_func_1_x2= float(diff_function_1[1])
coef_diff_func_1_x1= float(diff_function_1[2])
coef_diff_func_1_x0= float(diff_function_1[3])
diff_function_2 = np.polyfit(Span_y_2, int_dtheta_dz_2, 5) # sympy fit for dtheta_dz_2
coef_diff_func_2_x5= float(diff_function_2[0])
coef_diff_func_2_x4= float(diff_function_2[1])
coef_diff_func_2_x3= float(diff_function_2[2])
coef_diff_func_2_x2= float(diff_function_2[3])
coef_diff_func_2_x1= float(diff_function_2[4])
coef_diff_func_2_x0= float(diff_function_2[5])

def function_1(x): #approximation of dtheta_dz_1 (used for integrating in sympy)
    return coef_diff_func_1_x3*x**3+coef_diff_func_1_x2*x**2+coef_diff_func_1_x1*x+coef_diff_func_1_x0

def function_2(x): #approximation of dtheta_dz_2 (used for integrating in sympy)
    return coef_diff_func_2_x5*x**5+coef_diff_func_2_x4*x**4+coef_diff_func_2_x3*x**3+coef_diff_func_2_x2*x**2+coef_diff_func_2_x1*x+coef_diff_func_1_x0

x = smp.symbols('x', real = True)

result_1 = smp.integrate(function_1(x), x) # = theta_1
result_2 = smp.integrate(function_2(x), x) # = theta_2 - constant of integration
#print(result)

bora_func_1 = smp.lambdify([x], result_1) # converting sympy object to python function
bora_func_2 = smp.lambdify([x], result_2) # converting sympy object to python function

c_of_int_b = bora_func_1(0.35*0.5*Var.Span) - bora_func_2(0.35*0.5*Var.Span)

#Stan's method of integration
def j(y):
    int_fun_1 = d_dz_1_A*y**3 + d_dz_1_B*y**2 + d_dz_1_C*y + d_dz_1_D
    return int_fun_1

def k(y):
    int_fun_2 = d_dz_2_A*y**5 + d_dz_2_B*y**4 + d_dz_2_C*y**3 + d_dz_2_D*y**2 + d_dz_2_E*y + d_dz_2_F
    return int_fun_2

def l(y):
    int_fun_z = d_dz_z_A*y**7 + d_dz_z_B*y**6 + d_dz_z_C*y**5 + d_dz_z_D*y**4 + d_dz_z_E*y**3 + d_dz_z_F*y**2 + d_dz_z_G*y + d_dz_z_H
    return int_fun_z

int_z_z = []
for i in Span_y_z:
    estimatef1_1_z, errorf1_1_z = sp.integrate.quad(l, 0, i)
    int_z_z.append(estimatef1_1_z)

int_z_1 = []
for i in Span_y_1:
    estimatef1_1_1, errorf1_1_1 = sp.integrate.quad(j, 0, i)
    int_z_1.append(estimatef1_1_1)

int_z_2_no_c = []
for i in Span_y_2:
    estimatef1_1_1, errorf1_1_1 = sp.integrate.quad(k, Span_y_2[0], i)
    int_z_2_no_c.append(estimatef1_1_1)

c_of_int_s = int_z_1[-1]  - int_z_2_no_c[0] 

int_z_2 = []
for i in int_z_2_no_c:
    i += c_of_int_s
    int_z_2.append(i)

print(int_z_z[-1])
#print(c_of_int_s, c_of_int_b)

# plt.plot(Span_y_1, bora_func_1(Span_y_1[0:]))
# plt.plot(Span_y_2, bora_func_2(Span_y_2[0:]) + c_of_int_b)
plt.plot(Span_y_z, int_z_z)
#plt.plot(Span_y_1, int_z_1)
#plt.plot(Span_y_2, int_z_2)
plt.xlabel("Spanwise location")
plt.ylabel("theta(z)")
plt.show()