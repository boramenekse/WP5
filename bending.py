#import Moment_of_inertia_x_span_2 as momx
from Moment_of_inertia_x_span_2 import Fit_A, Fit_B, Fit_C, Fit_D
from math import *
import Variables_for_crossection_geometry as Var
from Variables_for_crossection_geometry import TnB_A, TnB_B, TnB_C, TnB_D, TnB_E, TnB_F
import numpy as np
import scipy as sp
import sympy as smp
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

#Span in a list
Span_y = np.arange(0.0, (Var.Span / 2) + 1, 1.0)

#CHECK:
#print(Span_y)

#METHOD 1 OF GETTING THE NON INTEGRATED FUNCTION
#function f(y) = -(M/E*MOI) split up in separate fractions, and put in a list 
def f_1(y):
    Mx_1 =  TnB_A*y**5
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_1 = -1*((Mx_1)/(Var.e_mod * Moixx))
    return int_fun_1

def f_2(y):
    Mx_2 = TnB_B*y**4
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_2 = -1*((Mx_2)/(Var.e_mod * Moixx))
    return int_fun_2

def f_3(y):
    Mx_3 = TnB_C*y**3
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_3 = -1*((Mx_3)/(Var.e_mod * Moixx))
    return int_fun_3

def f_4(y):
    Mx_4 = TnB_D*y**2
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_4 = -1*((Mx_4)/(Var.e_mod * Moixx))
    return int_fun_4

def f_5(y):
    Mx_5 = TnB_E*y
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_5 = -1*((Mx_5)/(Var.e_mod * Moixx))
    return int_fun_5

def f_6(y):
    Mx_6 = TnB_F
    Moixx = (Fit_A*y**3 + Fit_B*y**2 + Fit_C*y + Fit_D)
    int_fun_6 = -1*((Mx_6)/(Var.e_mod * Moixx))
    return int_fun_6

int_0 = []
for i in range(0, len(Span_y)):
    int_0_n = f_1(i) + f_2(i) + f_3(i) + f_4(i) + f_5(i) + f_6(i)
    int_0.append(int_0_n) 

#setting up test function and fitting
def test_function_0(x, A, B, C, D, E, F):
    y = A*x**5 + B*x**4 + C*x**3 + D*x**2 + E*x + F
    return y

parameters0, covariance0 = curve_fit(test_function_0, Span_y, int_0)

A0 = parameters0[0]
B0 = parameters0[1]
C0 = parameters0[2]
D0  =parameters0[3]
E0 = parameters0[4]
F0 = parameters0[5]

Fit_0 = test_function_0(Span_y, A0, B0, C0, D0, E0, F0)

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
def j(y):
    int_fun_0 = A0*y**5 + B0*y**4 + C0*y**3 + D0*y**2 + E0*y + F0
    return int_fun_0

int_1_2 = []
for i in range(0, len(Span_y / 2)):
    estimatef1_1_1, errorf1_1_1 = sp.integrate.quad(j, 0, i)
    int_1_2.append(estimatef1_1_1)

#METHOD 2.2 INTEGRATING WITH FITTED FUNCTION Fit_0 MANUALLY, REMEMBER CONSTANT IS 0
def k(y):
    int_1_manual_n = (1/6)*A0*y**6 + (1/5)*B0*y**5 + (1/4)*C0*y**4 + (1/3)*D0*y**3 + (1/2)*E0*y**2 + F0*y
    return int_1_manual_n

int_1_manual = k(Span_y)

#making the test function and fitting it for first integration:
def test_function_1(x, A, B, C, D, E, F, G):
    y = A*x**6 + B*x**5 + C*x**4 + D*x**3 + E*x**2 + F*x + G
    return y

parameters1, covariance1= curve_fit(test_function_1, Span_y, int_1_2) #change here int_1_(1/)2/manual to choose the used method

A1 = parameters1[0]
B1 = parameters1[1]
C1 = parameters1[2]
D1 = parameters1[3]
E1 = parameters1[4]
F1 = parameters1[5]
G1 = parameters1[6]

Fit_1 = test_function_1(Span_y, A1, B1, C1, D1, E1, F1, G1)

#CHECK:

#second integration
#METHOD 1 INTEGRATING WITH FITTED FUNCTION Fit_1 IN PYTHON
def g(y):
    int_fun_1 = A1*y**6 + B1*y**5 + C1*y**4 + D1*y**3 + E1*y**2 + F1*y + G1
    return int_fun_1

int_2_2 = []
for i in range(0, len(Span_y)):
    estimatef2, errorf2 = sp.integrate.quad(g, 0 , i)
    int_2_2.append(estimatef2)

#METHOD 2 INTEGRATING WITH FITTED FUNCTION Fit_1 MANUALLY REMEMBER CONSTANT IS 0
def l(y):
    int_2_manual_n = A1*(1/7)*y**7 + (1/6)*B1*y**6 + (1/5)*C1*y**5 + (1/4)*D1*y**4 + (1/3)*E1*y**3 + (1/2)*F1*y**2 + G1*y
    return int_2_manual_n

int_2_manual = l(Span_y)

#making the test function and fitting it for second integration:
def test_function_2(x, A, B, C, D, E, F, G, H):
    y = A*x**7 + B*x**6 + C*x**5 + D*x**4 + E*x**3 + F*x**2 + G*x + H
    return y

parameters2, covariance2 = curve_fit(test_function_2, Span_y, int_2_2) #change here int_2_2/manual to choose the used method

A2 = parameters2[0]
B2 = parameters2[1]
C2 = parameters2[2]
D2 = parameters2[3]
E2 = parameters2[4]
F2 = parameters2[5]
G2 = parameters2[6]
H2 = parameters2[7]

Fit_2 = test_function_2(Span_y, A2, B2, C2, D2, E2, F2, G2, H2)

#CHECK:

#total CHECK of derivatives
def h(y):
    orig_first_int_for_def = 7*A2*y**6 + 6*B2*y**5 + 5*C2*y**4 + 4*D2*y**3 + E2*3*y**2 + F2*y + G2
    return orig_first_int_for_def
orig_first_int_for = h(Span_y)

def i(y):
    orig_for_def = 6*7*A2*y**5 + 5*6*B2*y**4 + 4*5*C2*y**3 + 3*4*D2*y**2 + 2*E2*3*y + F2
    return orig_for_def
orig_for = i(Span_y)

#plt.plot(Span_y, orig_first_int_for, '-', label='First Deriv of defl')
#plt.plot(Span_y, orig_for, '-', label='Second Deriv of defl')

#printing deflection at the tip:
print('The deflection at the tip is:', test_function_2(Var.Span / 2, A2, B2, C2, D2, E2, F2, G2, H2), '[m]')

#plotting everything
#plt.plot(Span_y, int_0, 'o', label='Original Data')
#plt.plot(Span_y, int_1_manual, 'o', label='First Int Data') #change here int_1_1/2/manual to choose method
#plt.plot(Span_y, int_2_manual, 'o', label='Second Int Data') #change here int_2_1/manual to choose method
#plt.plot(Span_y, Fit_0, '-', label='Fit To Original Data')
#plt.plot(Span_y, Fit_1, '-', label='Fit To First Int Data')
plt.plot(Span_y, Fit_2, '-', label='Deflection')
plt.xlabel('Span y')
plt.ylabel('Deflection [m]')
plt.legend()
#plt.axis('equal')
plt.show()

#printing coefficients from deflection formula
#form of bending deflection formula is a 7th order polynomial 
#print('The coefficients of the bending deflection equation are;', A2, B2, C2, D2, E2, F2, G2, H2)