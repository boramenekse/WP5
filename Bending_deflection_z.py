import Moment_of_inertia_x_span as momx
from Moment_of_inertia_x_span import print_fit
from math import *
import Variables_for_crossection_geometry as Var
from Variables_for_crossection_geometry import print_load_cases
import numpy as np
import scipy as sp
import sympy as smp
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

result_t = print_load_cases()
result_fit = print_fit()
#Span in a list
Span_y = np.arange(0.0, (Var.Span / 2) + 1, 1.0)

#CHECK:
#print(Span_y)

#METHOD 1 OF GETTING THE NON INTEGRATED FUNCTION
#function f(y) = -(M/E*MOI) split up in separate fractions, and put in a list 
def f_1(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_1 = const*y**5
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_1 = -1*((Mx_1)/(Var.e_mod * Moixx))
    return int_fun_1

def f_2(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_2 = const*y**4
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_2 = -1*((Mx_2)/(Var.e_mod * Moixx))
    return int_fun_2

def f_3(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_3 = const*y**3
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_3 = -1*((Mx_3)/(Var.e_mod * Moixx))
    return int_fun_3

def f_4(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_4 = const*y**2
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_4 = -1*((Mx_4)/(Var.e_mod * Moixx))
    return int_fun_4

def f_5(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_5 = const*y**1
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_5 = -1*((Mx_5)/(Var.e_mod * Moixx))
    return int_fun_5

def f_6(y, const, fit_a, fit_b, fit_c, fit_d):
    Mx_6 = const
    Moixx = (fit_a*y**3 + fit_b*y**2 + fit_c*y + fit_d)
    int_fun_6 = -1*((Mx_6)/(Var.e_mod * Moixx))
    return int_fun_6

int_0 = []
for t in result_t:
    for fit in result_fit:
        int = []
        for i in range(0, len(Span_y)):
            int_0_n = f_1(i, t[0], fit[0], fit[1], fit[2], fit[3]) + f_2(i, t[1], fit[0], fit[1], fit[2], fit[3]) + f_3(i, t[2], fit[0], fit[1], fit[2], fit[3]) + f_4(i, t[3], fit[0], fit[1], fit[2], fit[3]) + f_5(i, t[4], fit[0], fit[1], fit[2], fit[3]) + f_6(i, t[5], fit[0], fit[1], fit[2], fit[3])
            int.append(int_0_n)
        int_0.append(int) 

#setting up test function and fitting
def test_function_0(x, A, B, C, D, E, F):
    y = A*x**5 + B*x**4 + C*x**3 + D*x**2 + E*x + F
    return y

parameters0, covariance0 = curve_fit(test_function_0, Span_y, int_0[0])
parameters1, covariance1 = curve_fit(test_function_0, Span_y, int_0[1])
parameters2, covariance2 = curve_fit(test_function_0, Span_y, int_0[2])
parameters3, covariance3 = curve_fit(test_function_0, Span_y, int_0[3])
parameters4, covariance4 = curve_fit(test_function_0, Span_y, int_0[4])
parameters5, covariance5 = curve_fit(test_function_0, Span_y, int_0[5])
parameters6, covariance6 = curve_fit(test_function_0, Span_y, int_0[6])
parameters7, covariance7 = curve_fit(test_function_0, Span_y, int_0[7])
parameters8, covariance8 = curve_fit(test_function_0, Span_y, int_0[8])
param_list = [parameters0, parameters1, parameters2, parameters3, parameters4, parameters5, parameters6, parameters7, parameters8]

A00 = parameters0[0]
B00 = parameters0[1]
C00 = parameters0[2]
D00  =parameters0[3]
E00 = parameters0[4]
F00 = parameters0[5]

A01 = parameters1[0]
B01 = parameters1[1]
C01 = parameters1[2]
D01  =parameters1[3]
E01 = parameters1[4]
F01 = parameters1[5]

A02 = parameters2[0]
B02 = parameters2[1]
C02 = parameters2[2]
D02  =parameters2[3]
E02 = parameters2[4]
F02 = parameters2[5]

A03 = parameters3[0]
B03 = parameters3[1]
C03 = parameters3[2]
D03  =parameters3[3]
E03 = parameters3[4]
F03 = parameters3[5]

A04 = parameters4[0]
B04 = parameters4[1]
C04 = parameters4[2]
D04  =parameters4[3]
E04 = parameters4[4]
F04 = parameters4[5]

A05 = parameters5[0]
B05 = parameters5[1]
C05 = parameters5[2]
D05  =parameters5[3]
E05 = parameters5[4]
F05 = parameters5[5]

A06 = parameters6[0]
B06 = parameters6[1]
C06 = parameters6[2]
D06  =parameters6[3]
E06 = parameters6[4]
F06 = parameters6[5]

A07 = parameters7[0]
B07 = parameters7[1]
C07 = parameters7[2]
D07  =parameters7[3]
E07 = parameters7[4]
F07 = parameters7[5]

A08 = parameters8[0]
B08 = parameters8[1]
C08 = parameters8[2]
D08  =parameters8[3]
E08 = parameters8[4]
F08 = parameters8[5]

Fit_00 = test_function_0(Span_y, A00, B00, C00, D00, E00, F00)
Fit_01 = test_function_0(Span_y, A01, B01, C01, D01, E01, F01)
Fit_02 = test_function_0(Span_y, A02, B02, C02, D02, E02, F02)
Fit_03 = test_function_0(Span_y, A03, B03, C03, D03, E03, F03)
Fit_04 = test_function_0(Span_y, A04, B04, C04, D04, E04, F04)
Fit_05 = test_function_0(Span_y, A05, B05, C05, D05, E05, F05)
Fit_06 = test_function_0(Span_y, A06, B06, C06, D06, E06, F06)
Fit_07 = test_function_0(Span_y, A07, B07, C07, D07, E07, F07)
Fit_08 = test_function_0(Span_y, A08, B08, C08, D08, E08, F08)

#CHECK:

#first integration
#METHOD 1 INTEGRATING WITH ORIGINAL FUNCTION
#estimatef1 = 0
#int_1_1 = []
#for i in range(0, len(Span_y)):
    #estimatef1_1, errorf1_1 = sp.integrate.quad(f_1,1,i)
    #estimatef1_2, errorf1_2 = sp.integrate.quad(f_2,1,i)
    #estimatef1_3, errorf1_3 = sp.integrate.quad(f_3,1,i)
    #estimatef1_4, errorf1_4 = sp.integrate.quad(f_4,1,i)
    #estimatef1_5, errorf1_5 = sp.integrate.quad(f_5,1,i)
    #estimatef1_6, errorf1_6 = sp.integrate.quad(f_6,1,i)
    #estimatef1 = estimatef1_1 + estimatef1_2 + estimatef1_3 + estimatef1_4 + estimatef1_5 + estimatef1
    #int_1_1.append(estimatef1)

#METHOD 2.1 INTEGRATING WITH FITTED FUNCTION Fit_0 IN PYTHON
def j(y, a0, b0, c0, d0, e0, f0):
    int_fun_0 = a0*y**5 + b0*y**4 + c0*y**3 + d0*y**2 + e0*y + f0
    return int_fun_0

int_1_2 = []
for par in param_list:
    list_t = [] 
    for i in range(0, len(Span_y / 2)):
        estimatef1_1_1, errorf1_1_1 = sp.integrate.quad(j, 0, i, args=(par[0], par[1], par[2], par[3], par[4], par[5]))
        list_t.append(estimatef1_1_1)
    int_1_2.append(list_t)

#METHOD 2.2 INTEGRATING WITH FITTED FUNCTION Fit_0 MANUALLY, REMEMBER CONSTANT IS 0
def k(y, a0, b0, c0, d0, e0, f0):
    int_1_manual_n = (1/6)*a0*y**6 + (1/5)*b0*y**5 + (1/4)*c0*y**4 + (1/3)*d0*y**3 + (1/2)*e0*y**2 + f0*y
    return int_1_manual_n

# int_1_manual = k(Span_y)

#making the test function and fitting it for first integration:
def test_function_1(x, A, B, C, D, E, F, G):
    y = A*x**6 + B*x**5 + C*x**4 + D*x**3 + E*x**2 + F*x + G
    return y

parameters1, covariance1= curve_fit(test_function_1, Span_y, int_1_2[0]) #change here int_1_(1/)2/manual to choose the used method
parameters2, covariance1= curve_fit(test_function_1, Span_y, int_1_2[1]) #change here int_1_(1/)2/manual to choose the used method
parameters3, covariance1= curve_fit(test_function_1, Span_y, int_1_2[2]) #change here int_1_(1/)2/manual to choose the used method
parameters4, covariance1= curve_fit(test_function_1, Span_y, int_1_2[3]) #change here int_1_(1/)2/manual to choose the used method
parameters5, covariance1= curve_fit(test_function_1, Span_y, int_1_2[4]) #change here int_1_(1/)2/manual to choose the used method
parameters6, covariance1= curve_fit(test_function_1, Span_y, int_1_2[5]) #change here int_1_(1/)2/manual to choose the used method
parameters7, covariance1= curve_fit(test_function_1, Span_y, int_1_2[6]) #change here int_1_(1/)2/manual to choose the used method
parameters8, covariance1= curve_fit(test_function_1, Span_y, int_1_2[7]) #change here int_1_(1/)2/manual to choose the used method
parameters9, covariance1= curve_fit(test_function_1, Span_y, int_1_2[8]) #change here int_1_(1/)2/manual to choose the used method
param_list_2 = [parameters1, parameters2, parameters3, parameters4, parameters5, parameters6, parameters7, parameters8, parameters9]

A11 = parameters1[0]
B11 = parameters1[1]
C11 = parameters1[2]
D11 = parameters1[3]
E11 = parameters1[4]
F11 = parameters1[5]
G11 = parameters1[6]

A12 = parameters2[0]
B12 = parameters2[1]
C12 = parameters2[2]
D12 = parameters2[3]
E12 = parameters2[4]
F12 = parameters2[5]
G12 = parameters2[6]

A13 = parameters3[0]
B13 = parameters3[1]
C13 = parameters3[2]
D13 = parameters3[3]
E13 = parameters3[4]
F13 = parameters3[5]
G13 = parameters3[6]

A14 = parameters4[0]
B14 = parameters4[1]
C14 = parameters4[2]
D14 = parameters4[3]
E14 = parameters4[4]
F14 = parameters4[5]
G14 = parameters4[6]

A15 = parameters5[0]
B15 = parameters5[1]
C15 = parameters5[2]
D15 = parameters5[3]
E15 = parameters5[4]
F15 = parameters5[5]
G15 = parameters5[6]

A16 = parameters6[0]
B16 = parameters6[1]
C16 = parameters6[2]
D16 = parameters6[3]
E16 = parameters6[4]
F16 = parameters6[5]
G16 = parameters6[6]

A17 = parameters7[0]
B17 = parameters7[1]
C17 = parameters7[2]
D17 = parameters7[3]
E17 = parameters7[4]
F17 = parameters7[5]
G17 = parameters7[6]

A18 = parameters8[0]
B18 = parameters8[1]
C18 = parameters8[2]
D18 = parameters8[3]
E18 = parameters8[4]
F18 = parameters8[5]
G18 = parameters8[6]

A19 = parameters9[0]
B19 = parameters9[1]
C19 = parameters9[2]
D19 = parameters9[3]
E19 = parameters9[4]
F19 = parameters9[5]
G19 = parameters9[6]

Fit_1 = test_function_1(Span_y, A11, B11, C11, D11, E11, F11, G11)
Fit_2 = test_function_1(Span_y, A12, B12, C12, D12, E12, F12, G12)
Fit_3 = test_function_1(Span_y, A13, B13, C13, D13, E13, F13, G13)
Fit_4 = test_function_1(Span_y, A14, B14, C14, D14, E14, F14, G14)
Fit_5 = test_function_1(Span_y, A15, B15, C15, D15, E15, F15, G15)
Fit_6 = test_function_1(Span_y, A16, B16, C16, D16, E16, F16, G16)
Fit_7 = test_function_1(Span_y, A17, B17, C17, D17, E17, F17, G17)
Fit_8 = test_function_1(Span_y, A18, B18, C18, D18, E18, F18, G18)
Fit_9 = test_function_1(Span_y, A19, B19, C19, D19, E19, F19, G19)

#CHECK:

#second integration
#METHOD 1 INTEGRATING WITH FITTED FUNCTION Fit_1 IN PYTHON
def g(y, a1, b1, c1, d1, e1, f1, g1):
    int_fun_1 = a1*y**6 + b1*y**5 + c1*y**4 + d1*y**3 + e1*y**2 + f1*y + g1
    return int_fun_1

int_2_2 = []
for par in param_list_2:
    int_list = []
    for i in range(0, len(Span_y)):
        estimatef2, errorf2 = sp.integrate.quad(g, 0 , i, args=(par[0], par[1], par[2], par[3], par[4], par[5], par[6]))
        int_list.append(estimatef2)
    int_2_2.append(int_list)

#METHOD 2 INTEGRATING WITH FITTED FUNCTION Fit_1 MANUALLY REMEMBER CONSTANT IS 0
def l(y, a1, b1, c1, d1, e1, f1, g1):
    int_2_manual_n = a1*(1/7)*y**7 + (1/6)*b1*y**6 + (1/5)*c1*y**5 + (1/4)*d1*y**4 + (1/3)*e1*y**3 + (1/2)*f1*y**2 + g1*y
    return int_2_manual_n

# int_2_manual = l(Span_y)

#making the test function and fitting it for second integration:
def test_function_2(x, A, B, C, D, E, F, G, H):
    y = A*x**7 + B*x**6 + C*x**5 + D*x**4 + E*x**3 + F*x**2 + G*x + H
    return y

parameters21, covariance21 = curve_fit(test_function_2, Span_y, int_2_2[0]) #change here int_2_2/manual to choose the used method
parameters22, covariance22 = curve_fit(test_function_2, Span_y, int_2_2[1]) #change here int_2_2/manual to choose the used method
parameters23, covariance23 = curve_fit(test_function_2, Span_y, int_2_2[2]) #change here int_2_2/manual to choose the used method
parameters24, covariance24 = curve_fit(test_function_2, Span_y, int_2_2[3]) #change here int_2_2/manual to choose the used method
parameters25, covariance25 = curve_fit(test_function_2, Span_y, int_2_2[4]) #change here int_2_2/manual to choose the used method
parameters26, covariance26 = curve_fit(test_function_2, Span_y, int_2_2[5]) #change here int_2_2/manual to choose the used method
parameters27, covariance27 = curve_fit(test_function_2, Span_y, int_2_2[6]) #change here int_2_2/manual to choose the used method
parameters28, covariance28 = curve_fit(test_function_2, Span_y, int_2_2[7]) #change here int_2_2/manual to choose the used method
parameters29, covariance29 = curve_fit(test_function_2, Span_y, int_2_2[8]) #change here int_2_2/manual to choose the used method

A21 = parameters21[0]
B21 = parameters21[1]
C21 = parameters21[2]
D21 = parameters21[3]
E21 = parameters21[4]
F21 = parameters21[5]
G21 = parameters21[6]
H21 = parameters21[7]

A22 = parameters22[0]
B22 = parameters22[1]
C22 = parameters22[2]
D22 = parameters22[3]
E22 = parameters22[4]
F22 = parameters22[5]
G22 = parameters22[6]
H22 = parameters22[7]

A23 = parameters23[0]
B23 = parameters23[1]
C23 = parameters23[2]
D23 = parameters23[3]
E23 = parameters23[4]
F23 = parameters23[5]
G23 = parameters23[6]
H23 = parameters23[7]

A24 = parameters24[0]
B24 = parameters24[1]
C24 = parameters24[2]
D24 = parameters24[3]
E24 = parameters24[4]
F24 = parameters24[5]
G24 = parameters24[6]
H24 = parameters24[7]

A25 = parameters25[0]
B25 = parameters25[1]
C25 = parameters25[2]
D25 = parameters25[3]
E25 = parameters25[4]
F25 = parameters25[5]
G25 = parameters25[6]
H25 = parameters25[7]

A26 = parameters26[0]
B26 = parameters26[1]
C26 = parameters26[2]
D26 = parameters26[3]
E26 = parameters26[4]
F26 = parameters26[5]
G26 = parameters26[6]
H26 = parameters26[7]

A27 = parameters27[0]
B27 = parameters27[1]
C27 = parameters27[2]
D27 = parameters27[3]
E27 = parameters27[4]
F27 = parameters27[5]
G27 = parameters27[6]
H27 = parameters27[7]

A28 = parameters28[0]
B28 = parameters28[1]
C28 = parameters28[2]
D28 = parameters28[3]
E28 = parameters28[4]
F28 = parameters28[5]
G28 = parameters28[6]
H28 = parameters28[7]

A29 = parameters29[0]
B29 = parameters29[1]
C29 = parameters29[2]
D29 = parameters29[3]
E29 = parameters29[4]
F29 = parameters29[5]
G29 = parameters29[6]
H29 = parameters29[7]

Fit_21 = test_function_2(Span_y, A21, B21, C21, D21, E21, F21, G21, H21)
Fit_22 = test_function_2(Span_y, A22, B22, C22, D22, E22, F22, G22, H22)
Fit_23 = test_function_2(Span_y, A23, B23, C23, D23, E23, F23, G23, H23)
Fit_24 = test_function_2(Span_y, A24, B24, C24, D24, E24, F24, G24, H24)
Fit_25 = test_function_2(Span_y, A25, B25, C25, D25, E25, F25, G25, H25)
Fit_26 = test_function_2(Span_y, A26, B26, C26, D26, E26, F26, G26, H26)
Fit_27 = test_function_2(Span_y, A27, B27, C27, D27, E27, F27, G27, H27)
Fit_28 = test_function_2(Span_y, A28, B28, C28, D28, E28, F28, G28, H28)
Fit_29 = test_function_2(Span_y, A29, B29, C29, D29, E29, F29, G29, H29)

#CHECK:

#total CHECK of derivatives
def h(y, a2, b2, c2, d2, e2, f2, g2):
    orig_first_int_for_def = 7*a2*y**6 + 6*b2*y**5 + 5*c2*y**4 + 4*d2*y**3 + e2*3*y**2 + f2*y + g2
    return orig_first_int_for_def
# orig_first_int_for = h(Span_y)

def i(y, a2, b2, c2, d2, e2, f2):
    orig_for_def = 6*7*a2*y**5 + 5*6*b2*y**4 + 4*5*c2*y**3 + 3*4*d2*y**2 + 2*e2*3*y + f2
    return orig_for_def
# orig_for = i(Span_y)

#plt.plot(Span_y, orig_first_int_for, '-', label='First Deriv of defl')
#plt.plot(Span_y, orig_for, '-', label='Second Deriv of defl')

#printing deflection at the tip:
print('The deflection at the tip is for 1:', test_function_2(Var.Span / 2, A21, B21, C21, D21, E21, F21, G21, H21), '[m]')
print('The deflection at the tip is for 2:', test_function_2(Var.Span / 2, A22, B22, C22, D22, E22, F22, G22, H22), '[m]')
print('The deflection at the tip is for 3:', test_function_2(Var.Span / 2, A23, B23, C23, D23, E23, F23, G23, H23), '[m]')
print('The deflection at the tip is for 4:', test_function_2(Var.Span / 2, A24, B24, C24, D24, E24, F24, G24, H24), '[m]')
print('The deflection at the tip is for 5:', test_function_2(Var.Span / 2, A25, B25, C25, D25, E25, F25, G25, H25), '[m]')
print('The deflection at the tip is for 6:', test_function_2(Var.Span / 2, A26, B26, C26, D26, E26, F26, G26, H26), '[m]')
print('The deflection at the tip is for 7:', test_function_2(Var.Span / 2, A27, B27, C27, D27, E27, F27, G27, H27), '[m]')
print('The deflection at the tip is for 8:', test_function_2(Var.Span / 2, A28, B28, C28, D28, E28, F28, G28, H28), '[m]')
print('The deflection at the tip is for 9:', test_function_2(Var.Span / 2, A29, B29, C29, D29, E29, F29, G29, H29), '[m]')

#plotting everything
#plt.plot(Span_y, int_0, 'o', label='Original Data')
#plt.plot(Span_y, int_1_manual, 'o', label='First Int Data') #change here int_1_1/2/manual to choose method
#plt.plot(Span_y, int_2_manual, 'o', label='Second Int Data') #change here int_2_1/manual to choose method
#plt.plot(Span_y, Fit_0, '-', label='Fit To Original Data')
#plt.plot(Span_y, Fit_1, '-', label='Fit To First Int Data')
plt.figure()
plt.title('Critical Loading Case 1')
plt.xlabel('Spanwise Location [m]')
plt.ylabel('Deflection [m]')
plt.plot(Span_y, Fit_21, '-', label='Deflection')
plt.plot(Span_y, Fit_22, '-', label='Deflection')
plt.plot(Span_y, Fit_23, '-', label='Deflection')
plt.legend(['philosophy 1', 'philosophy 2', 'philosophy 3'], loc='upper left', frameon=True)
plt.grid()
plt.show()
plt.figure()
plt.title('Critical Loading Case 2')
plt.xlabel('Spanwise Location [m]')
plt.ylabel('Deflection [m]')
plt.plot(Span_y, Fit_24, '-', label='Deflection')
plt.plot(Span_y, Fit_25, '-', label='Deflection')
plt.plot(Span_y, Fit_26, '-', label='Deflection')
plt.legend(['philosophy 1', 'philosophy 2', 'philosophy 3'], loc='upper right', frameon=True)
plt.grid()
plt.show()
plt.figure()
plt.title('Critical Loading Case 3')
plt.xlabel('Spanwise Location [m]')
plt.ylabel('Deflection [m]')
plt.plot(Span_y, Fit_27, '-', label='Deflection')
plt.plot(Span_y, Fit_28, '-', label='Deflection')
plt.plot(Span_y, Fit_29, '-', label='Deflection')
plt.legend(['philosophy 1', 'philosophy 2', 'philosophy 3'], loc='upper left', frameon=True)
plt.grid()
plt.show()

#printing coefficients from deflection formula
#form of bending deflection formula is a 7th order polynomial 
#print('The coefficients of the bending deflection equation are;', A2, B2, C2, D2, E2, F2, G2, H2)

