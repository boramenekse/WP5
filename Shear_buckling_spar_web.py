import math
import numpy as np
import matplotlib.pyplot as plt
#also check compressive failure for this component

t_fs_root_1 = 0.05
t_rs_root_1 = 0.05
t_fs_root_2 = 0.05
t_rs_root_2 = 0.05
t_fs_root_3 = 0.05
t_rs_root_3 = 0.05

### SETTING SOME CONSTANTS VOR THE FORMULAE

c_root = 11.59
c_tip = 3.25
wingspan = 66.78

pi = math.pi            # get pi from the math repository
E = 69 * (10 ** 6)      # define the young's modulus, in Pa
v = 0.33                # define the poissons ratio, no unit
k_s = 9.6               # define the k_s constant, found in figure 16 from the reader, no unit


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
#plt.plot(y_list, critical_shear_stress_forward_spar_list_1)
#plt.plot(y_list, critical_shear_stress_rear_spar_list_1)
#plt.legend(["forward", "rear"])
#plt.show()

##### GET THE SHEAR STRESS DUE TO SHEAR FORCE

def get_shear_functions(y_list):
    LC_8_shear = []
    LC_12_shear = []
    LC_16_shear = []
    for y in y_list:
        if y < 11.69:
            LC_8_shear.append(-1.7643233643147067 * y**4 + 26.55083716939181 * y**3 + 3358.583972792644 * y**2 + -27822.798128670238 * y + -1225433.3783305446)
            LC_12_shear.append(0.9784933304341741 * y ** 4 + -26.72072199357151 * y ** 3 + -1042.9512097953173 * y ** 2 + 7702.078406464768 * y + 529912.4293065071)
            LC_16_shear.append(-2.0008036049323445 * y ** 4 + 37.59176895469123 * y ** 3 + 3142.1526531781406 * y ** 2 + -13417.741226393178 * y + -1582191.8673072504)
        else:
            LC_8_shear.append(-1.7643233643147067 * y**4 + 26.55083716939181 * y**3 + 3358.583972792644 * y**2 + -27822.798128670238 * y + -1610817.128330534)
            LC_12_shear.append(0.9784933304341741 * y ** 4 + -26.72072199357151 * y ** 3 + -1042.9512097953173 * y ** 2 + 7702.078406464768 * y + 684065.9293065054)
            LC_16_shear.append(-2.0008036049323445 * y ** 4 + 37.59176895469123 * y ** 3 + 3142.1526531781406 * y ** 2 + -13417.741226393178 * y + -1967575.617307226)
    return LC_8_shear, LC_12_shear, LC_16_shear

LC_8_shear, LC_12_shear, LC_16_shear = get_shear_functions(y_list)

#check
#plt.plot(y_list, LC_8_shear)
#plt.plot(y_list, LC_12_shear)
#plt.plot(y_list, LC_16_shear)
#plt.show()

average_shear_stress_8_1_list = []
average_shear_stress_8_2_list = []
average_shear_stress_8_3_list = []
average_shear_stress_12_1_list = []
average_shear_stress_12_2_list = []
average_shear_stress_12_3_list = []
average_shear_stress_16_1_list = []
average_shear_stress_16_2_list = []
average_shear_stress_16_3_list = []

for i in range(len(y_list)):
    y = y_list[i]
    tau_average_8_1 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    tau_average_8_2 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    tau_average_8_3 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))
    tau_average_12_1 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    tau_average_12_2 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    tau_average_12_3 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))
    tau_average_16_1 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    tau_average_16_2 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    tau_average_16_3 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))

    average_shear_stress_8_1_list.append(tau_average_8_1)
    average_shear_stress_8_2_list.append(tau_average_8_2)
    average_shear_stress_8_3_list.append(tau_average_8_3)
    average_shear_stress_12_1_list.append(tau_average_12_1)
    average_shear_stress_12_2_list.append(tau_average_12_2)
    average_shear_stress_12_3_list.append(tau_average_12_3)
    average_shear_stress_16_1_list.append(tau_average_16_1)
    average_shear_stress_16_2_list.append(tau_average_16_2)
    average_shear_stress_16_3_list.append(tau_average_16_3)

