import math
import numpy as np
import matplotlib.pyplot as plt
#also check compressive failure for this component

### SETTING SOME CONSTANTS VOR THE FORMULAE

pi = math.pi            # get pi from the math repository
E = 69 * (10 ** 6)      # define the young's modulus, in Pa
v = 0.33                # define the poissons ratio, no unit
k_s = 9.6               # define the k_s constant, found in figure 16 from the reader, no unit

c_root = 11.59
c_tip = 3.25
wingspan = 66.78

t_fs_root_1 = 0.05
t_rs_root_1 = 0.05
t_fs_root_2 = 0.05
t_rs_root_2 = 0.05
t_fs_root_3 = 0.05
t_rs_root_3 = 0.05

def get_chord(y):
    c = c_root - y * (c_root-c_tip) / (wingspan/2)
    return c

def get_forward_spar_lenth(y):
    l_fs = (0.045100 + 0.045200) * get_chord(y)
    return l_fs

def get_rear_spar_lenth(y):
    l_rs = (0.033900 + 0.016200) * get_chord(y)
    return l_rs

def get_thickness_forward_spar_1(y):
    t_fs_1 = t_fs_root_1 * get_chord(y) / get_chord(0)
    return t_fs_1

def get_thickness_rear_spar_1(y):
    t_rs_1 = t_rs_root_1 * get_chord(y) / get_chord(0)
    return t_rs_1

def get_thickness_forward_spar_2(y):
    t_fs_2 = t_fs_root_2 * get_chord(y) / get_chord(0)
    return t_fs_2

def get_thickness_rear_spar_2(y):
    t_rs_2 = t_rs_root_2 * get_chord(y) / get_chord(0)
    return t_rs_2

def get_thickness_forward_spar_3(y):
    t_fs_3 = t_fs_root_3 * get_chord(y) / get_chord(0)
    return t_fs_3

def get_thickness_rear_spar_3(y):
    t_rs_3 = t_rs_root_3 * get_chord(y) / get_chord(0)
    return t_rs_3

#get the ylist to do everything with
y_list = np.linspace(0, wingspan/2, num=1000, endpoint=True)

#### GET THE CRITICAL SHEAR STRESS FOR BOTH SPARS

def get_critical_shear_stress(pi, b, k_s, t, E, v):
    t_cr = (pi**2 * k_s * E) / (12 * (1 - v ** 2)) * (t/b)**2
    return t_cr

critical_shear_stress_forward_spar_list_1 = []
critical_shear_stress_rear_spar_list_1 = []
critical_shear_stress_forward_spar_list_2 = []
critical_shear_stress_rear_spar_list_2 = []
critical_shear_stress_forward_spar_list_3 = []
critical_shear_stress_rear_spar_list_3 = []
for y in y_list:
    critical_shear_stress_forward_spar_list_1.append(get_critical_shear_stress(pi, get_forward_spar_lenth(y), k_s, get_thickness_forward_spar_1(y), E, v))
    critical_shear_stress_rear_spar_list_1.append(get_critical_shear_stress(pi, get_rear_spar_lenth(y), k_s, get_thickness_rear_spar_1(y), E, v))
    critical_shear_stress_forward_spar_list_2.append(get_critical_shear_stress(pi, get_forward_spar_lenth(y), k_s, get_thickness_forward_spar_2(y), E, v))
    critical_shear_stress_rear_spar_list_2.append(get_critical_shear_stress(pi, get_rear_spar_lenth(y), k_s, get_thickness_rear_spar_2(y), E, v))
    critical_shear_stress_forward_spar_list_3.append(get_critical_shear_stress(pi, get_forward_spar_lenth(y), k_s, get_thickness_forward_spar_3(y), E, v))
    critical_shear_stress_rear_spar_list_3.append(get_critical_shear_stress(pi, get_rear_spar_lenth(y), k_s, get_thickness_rear_spar_3(y), E, v))

#check the calculations
plt.plot(y_list, critical_shear_stress_forward_spar_list_1)
plt.plot(y_list, critical_shear_stress_rear_spar_list_1)
plt.legend(["forward", "rear"])
plt.show()