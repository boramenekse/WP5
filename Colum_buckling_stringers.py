#also check compressive failure for this component
import sympy as smp
import Variables_for_crossection_geometry as Var
import matplotlib.pyplot as plt
from math import *
import Normal_stress_distribution as str
import numpy as np

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

stress_dist_fun = smp.lambdify([y], str.stress_16_1.subs([(z, str.br_z1), (x, str.br_x1)]).simplify())

#Moment of inertia calculation X axis for stringer
#centroid z direction
Str_z_tilda_= (Var.Str_v_len*Var.Str_v_th*Var.Str_v_len/2)/Var.Str_A
Str_MOI_X = (Var.Str_h_len*Var.Str_h_th**3/12 + Var.Str_h_len*Var.Str_h_th*(Str_z_tilda_ - Var.Str_h_th/2)**2) + (Var.Str_v_len**3*Var.Str_v_th/12 + Var.Str_v_len*Var.Str_v_th * (Var.Str_v_len/2-Str_z_tilda_)**2) 

#Calculating critical stress
def critical_stress_b_Str(L):
    sigma = (Var.k * pi**2 *Var.e_mod *Str_MOI_X)/(L**2*Var.Str_A)
    return sigma

#Check
#print(critical_stress_b_Str(Var.Span/(2*30))/1e6, "MPa")

#Ribs spacing and ribs location along half-span
ribs_spacing_list = [0.534261602642555, 0.525745792613091, 0.517948192133515, 0.510807981088931, 0.504278027391229, 0.498322565969259, 0.492915713953116, 0.488040637415200, 0.483689285827142, 0.479862693176665, 0.476571931811814, 0.473839921959174, 0.471704487447047, 0.470223387542136, 0.469482727462029, 0.469611605477857, 0.470809336470809, 0.707830742064340, 0.719502430751495, 1.23647866009611, 0.644328916142815, 0.666243412755113, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043, 1.22935110850043]
ribs_location_list = [0.543573620188936, 1.07783522283149, 1.60358101544458, 2.12152920757810, 2.63233718866703, 3.13661521605826, 3.63493778202752, 4.12785349598063, 4.61589413339583, 5.09958341922297, 5.57944611239964, 6.05601804421145, 6.52985796617063, 7.00156245361767, 7.47178584115981, 7.94126856862184, 8.41088017409970, 8.88168951057050, 9.58952025263484, 10.3090226833863, 11.5455013434825, 12.1898302596253, 12.8560736723804, 14.0854247808808, 15.3147758893812, 16.5441269978817, 17.7734781063821, 19.0028292148825, 20.2321803233829, 21.4615314318834, 22.6908825403838, 23.9202336488842, 25.1495847573847, 26.3789358658851, 27.6082869743855, 28.8376380828859, 30.0669891913864, 31.2963402998868, 32.5256914083872]
span = np.linspace(0, Var.Span*0.5, endpoint=True)

critical_stress_list = []
critical_buck_stress=[]
for i in range(len(ribs_spacing_list)):
    sigma_buck = critical_stress_b_Str(ribs_spacing_list[i])
    critical_buck_stress.append(sigma_buck)
    critical_buck_stress.append(sigma_buck)
    critical_stress_list.append(sigma_buck)

final_position_list = []
for i in range(len(ribs_location_list)-1):
    final_position_list.append(ribs_location_list[i])
    final_position_list.append(ribs_location_list[i+1])

#print(final_position_list)

#Check
#print(critical_buck_stress)
plt.figure()
#plt.plot(ribs_location_list[:38], critical_buck_stress, marker='o')
# plt.plot(final_position_list, critical_buck_stress)
plt.plot(span, -stress_dist_fun(span[0:]))
plt.xlabel("Spanwise location")
plt.ylabel("Critical buckling stress [MPA]")
plt.show()

margin_of_safety_list = []
list_position = []
for i in range(len(ribs_location_list)-1):
    start = ribs_location_list[i]
    end = ribs_location_list[i+1]
    average = 0.5*(start+end)
    critical_stress_value = critical_stress_list[i]
    ratio_list = []
    ratio_start = -critical_stress_value/stress_dist_fun(start)
    ratio_end = -critical_stress_value/stress_dist_fun(end)
    ratio_average = -critical_stress_value/stress_dist_fun(average)
    ratio_list.append(ratio_start)
    ratio_list.append(ratio_average)
    ratio_list.append(ratio_end)
    ratio = min(ratio_list)
    if ratio == ratio_start:
        list_position.append(start)
    elif ratio == ratio_average:
        list_position.append(average)
    else:
        list_position.append(end)
    margin_of_safety_list.append(ratio)

print(margin_of_safety_list)
print(len(margin_of_safety_list))

plt.figure()
plt.plot(list_position, margin_of_safety_list)
plt.show()