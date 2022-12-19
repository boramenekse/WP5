import math
import numpy as np
import matplotlib.pyplot as plt
#also check compressive failure for this component

t_fs_root_1 = 0.025
t_rs_root_1 = 0.025
t_fs_root_2 = 0.04
t_rs_root_2 = 0.04
t_fs_root_3 = 0.015
t_rs_root_3 = 0.015

### SETTING SOME CONSTANTS VOR THE FORMULAE

c_root = 11.59
c_tip = 3.25
wingspan = 66.78

pi = math.pi            # get pi from the math repository
E = 69 * (10 ** 9)      # define the young's modulus, in Pa
v = 0.33                # define the poissons ratio, no unit
k_s = 9.6               # define the k_s constant, found in figure 16 from the reader, no unit
k_v = 1.3


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
plt.plot(y_list, critical_shear_stress_forward_spar_list_2)
plt.plot(y_list, critical_shear_stress_rear_spar_list_2)
plt.plot(y_list, critical_shear_stress_forward_spar_list_3)
plt.plot(y_list, critical_shear_stress_rear_spar_list_3)
plt.legend(["forward_1", "rear_1", "forward_2", "rear_2", "forward_3", "rear_3"])
plt.xlabel('Span [m]')
plt.ylabel('Critical Buckling Shear Stress [Pa]')
plt.grid()
plt.show()

##### GET THE SHEAR STRESS DUE TO SHEAR FORCE

def get_shear_functions(y_list):
    LC_8_shear = []
    LC_12_shear = []
    LC_16_shear = []
    for y in y_list:
        if y < 11.69:
            #LC_8_shear.append(-1.7643233643147067 * y**4 + 26.55083716939181 * y**3 + 3358.583972792644 * y**2 + -27822.798128670238 * y + -1225433.3783305446)
            LC_12_shear.append(0.9784933304341741 * y ** 4 + -26.72072199357151 * y ** 3 + -1042.9512097953173 * y ** 2 + 7702.078406464768 * y + 529912.4293065071)
            LC_16_shear.append(-2.0008036049323445 * y ** 4 + 37.59176895469123 * y ** 3 + 3142.1526531781406 * y ** 2 + -13417.741226393178 * y + -1582191.8673072504)
        else:
            #LC_8_shear.append(-1.7643233643147067 * y**4 + 26.55083716939181 * y**3 + 3358.583972792644 * y**2 + -27822.798128670238 * y + -1610817.128330534)
            LC_12_shear.append(0.9784933304341741 * y ** 4 + -26.72072199357151 * y ** 3 + -1042.9512097953173 * y ** 2 + 7702.078406464768 * y + 684065.9293065054)
            LC_16_shear.append(-2.0008036049323445 * y ** 4 + 37.59176895469123 * y ** 3 + 3142.1526531781406 * y ** 2 + -13417.741226393178 * y + -1967575.617307226)
    return LC_8_shear, LC_12_shear, LC_16_shear

LC_8_shear, LC_12_shear, LC_16_shear = get_shear_functions(y_list)

#check the calculations
#plt.plot(y_list, LC_8_shear)
#plt.plot(y_list, LC_12_shear)
#plt.plot(y_list, LC_16_shear)
#plt.show()

#average_shear_stress_8_1_list = []
#average_shear_stress_8_2_list = []
#average_shear_stress_8_3_list = []
average_shear_stress_12_1_list = []
average_shear_stress_12_2_list = []
average_shear_stress_12_3_list = []
average_shear_stress_16_1_list = []
average_shear_stress_16_2_list = []
average_shear_stress_16_3_list = []

for i in range(len(y_list)):
    y = y_list[i]
    #tau_average_8_1 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    #tau_average_8_2 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    #tau_average_8_3 = LC_8_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))
    tau_average_12_1 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    tau_average_12_2 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    tau_average_12_3 = LC_12_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))
    tau_average_16_1 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_1(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_1(y))
    tau_average_16_2 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_2(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_2(y))
    tau_average_16_3 = LC_16_shear[i] / (get_forward_spar_lenth(y) * get_thickness_forward_spar_3(y) + get_rear_spar_lenth(y) * get_thickness_rear_spar_3(y))

    #average_shear_stress_8_1_list.append((tau_average_8_1))
    #average_shear_stress_8_2_list.append((tau_average_8_2))
    #average_shear_stress_8_3_list.append((tau_average_8_3))
    average_shear_stress_12_1_list.append((tau_average_12_1))
    average_shear_stress_12_2_list.append((tau_average_12_2))
    average_shear_stress_12_3_list.append((tau_average_12_3))
    average_shear_stress_16_1_list.append((tau_average_16_1))
    average_shear_stress_16_2_list.append((tau_average_16_2))
    average_shear_stress_16_3_list.append((tau_average_16_3))

    #maximum_shear_stress_8_1_list = [i * k_v for i in average_shear_stress_8_1_list]
    #maximum_shear_stress_8_2_list = [i * k_v for i in average_shear_stress_8_2_list]
    #maximum_shear_stress_8_3_list = [i * k_v for i in average_shear_stress_8_3_list]
    maximum_shear_stress_12_1_list = [i * k_v for i in average_shear_stress_12_1_list]
    maximum_shear_stress_12_2_list = [i * k_v for i in average_shear_stress_12_2_list]
    maximum_shear_stress_12_3_list = [i * k_v for i in average_shear_stress_12_3_list]
    maximum_shear_stress_16_1_list = [i * k_v for i in average_shear_stress_16_1_list]
    maximum_shear_stress_16_2_list = [i * k_v for i in average_shear_stress_16_2_list]
    maximum_shear_stress_16_3_list = [i * k_v for i in average_shear_stress_16_3_list]

#check the calculations
#plt.plot(y_list, maximum_shear_stress_8_1_list)
#plt.plot(y_list, maximum_shear_stress_8_2_list)
#plt.plot(y_list, maximum_shear_stress_8_3_list)
#plt.plot(y_list, maximum_shear_stress_12_1_list)
#plt.plot(y_list, maximum_shear_stress_12_2_list)
#plt.plot(y_list, maximum_shear_stress_12_3_list)
#plt.plot(y_list, maximum_shear_stress_16_1_list)
#plt.plot(y_list, maximum_shear_stress_16_2_list)
#plt.plot(y_list, maximum_shear_stress_16_3_list)
#plt.legend(["av,tao_8_1", "av,tao_8_2", "av,tao_8_3", "av,tao_12_1", "av,tao_12_2", "av,tao_12_3", "av,tao_16_1", "av,tao_16_2", "av,tao_16_3"])
#plt.show()


###CONTRIBUTION OF TORQUE

def get_LC_8_torque(y):
    if y < 11.69:
        return -4.936650618346473 * y**4 + 335.5436542921075 * y**3 + 8009.420697124085 * y**2 + -1010097.9604714005 * y + 15577734.006808802
    else:
        return -4.936650618384433 * y**4 + 335.5436542929485 * y**3 + 8009.420697118811 * y**2 + -1010097.9604714005 * y + 18442679.692323517

def get_LC_12_torque(y):
    if y < 11.69:
        return 1.109847947180169 * y**4 + -162.7034547538029 * y**3 + 12610.129687494236 * y**2 + -481202.06432541495 * y + 9929226.911220081
    else:
        return 1.109847947180169 * y**4 + -162.7034547538029 * y**3 + 12610.129687491217 * y**2 + -481202.0643254161 * y + 6685718.310560216

def get_LC_16_torque(y):
    if y < 11.69:
        return -4.0265674601192325 * y**4 + 336.1076467822765 * y**3 + -1734.6018412792378 * y**2 + -458768.07825406635 * y + 6880148.502473444
    else:
        return -4.0265674601192325 * y**4 + 336.1076467822765 * y**3 + -1734.6018412792378 * y**2 + -458768.0782540611 * y + 9745094.187988076

def get_inclosed_area(y):
    return (get_forward_spar_lenth(y) + get_rear_spar_lenth(y)) * 0.55 * get_chord(y)

#shear_stress_torque_LC_8_1_forward = []
#shear_stress_torque_LC_8_1_rear = []
#shear_stress_torque_LC_8_2_forward = []
#shear_stress_torque_LC_8_2_rear = []
#shear_stress_torque_LC_8_3_forward = []
#shear_stress_torque_LC_8_3_rear = []
shear_stress_torque_LC_12_1_forward = []
shear_stress_torque_LC_12_1_rear = []
shear_stress_torque_LC_12_2_forward = []
shear_stress_torque_LC_12_2_rear = []
shear_stress_torque_LC_12_3_forward = []
shear_stress_torque_LC_12_3_rear = []
shear_stress_torque_LC_16_1_forward = []
shear_stress_torque_LC_16_1_rear = []
shear_stress_torque_LC_16_2_forward = []
shear_stress_torque_LC_16_2_rear = []
shear_stress_torque_LC_16_3_forward = []
shear_stress_torque_LC_16_3_rear = []

for y in y_list:
    #shear_stress_torque_LC_8_1_forward.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_1(y)))
    #shear_stress_torque_LC_8_1_rear.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_1(y)))
    #shear_stress_torque_LC_8_2_forward.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_2(y)))
    #shear_stress_torque_LC_8_2_rear.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_2(y)))
    #shear_stress_torque_LC_8_3_forward.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_3(y)))
    #shear_stress_torque_LC_8_3_rear.append(get_LC_8_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_3(y)))
    shear_stress_torque_LC_12_1_forward.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_1(y)))
    shear_stress_torque_LC_12_1_rear.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_1(y)))
    shear_stress_torque_LC_12_2_forward.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_2(y)))
    shear_stress_torque_LC_12_2_rear.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_2(y)))
    shear_stress_torque_LC_12_3_forward.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_3(y)))
    shear_stress_torque_LC_12_3_rear.append(get_LC_12_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_3(y)))
    shear_stress_torque_LC_16_1_forward.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_1(y)))
    shear_stress_torque_LC_16_1_rear.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_1(y)))
    shear_stress_torque_LC_16_2_forward.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_2(y)))
    shear_stress_torque_LC_16_2_rear.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_2(y)))
    shear_stress_torque_LC_16_3_forward.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_forward_spar_3(y)))
    shear_stress_torque_LC_16_3_rear.append(get_LC_16_torque(y) / (2 * get_inclosed_area(y) * get_thickness_rear_spar_3(y)))

#plt.plot(y_list, shear_stress_torque_LC_8_1_forward)
#plt.plot(y_list, shear_stress_torque_LC_8_1_rear)
#plt.plot(y_list, shear_stress_torque_LC_8_2_forward)
#plt.plot(y_list, shear_stress_torque_LC_8_2_rear)
#plt.plot(y_list, shear_stress_torque_LC_8_3_forward)
#plt.plot(y_list, shear_stress_torque_LC_8_3_rear)
#plt.plot(y_list, shear_stress_torque_LC_12_1_forward)
#plt.plot(y_list, shear_stress_torque_LC_12_1_rear)
#plt.plot(y_list, shear_stress_torque_LC_12_2_forward)
#plt.plot(y_list, shear_stress_torque_LC_12_2_rear)
#plt.plot(y_list, shear_stress_torque_LC_12_3_forward)
#plt.plot(y_list, shear_stress_torque_LC_12_3_rear)
#plt.plot(y_list, shear_stress_torque_LC_16_1_forward)
#plt.plot(y_list, shear_stress_torque_LC_16_1_rear)
#plt.plot(y_list, shear_stress_torque_LC_16_2_forward)
#plt.plot(y_list, shear_stress_torque_LC_16_2_rear)
#plt.plot(y_list, shear_stress_torque_LC_16_3_forward)
#plt.plot(y_list, shear_stress_torque_LC_16_3_rear)
#plt.show()


#shear_stress_LC_8_1_forward = []
#shear_stress_LC_8_1_rear = []
#shear_stress_LC_8_2_forward = []
#shear_stress_LC_8_2_rear = []
#shear_stress_LC_8_3_forward = []
#shear_stress_LC_8_3_rear = []
shear_stress_LC_12_1_forward = []
shear_stress_LC_12_1_rear = []
shear_stress_LC_12_2_forward = []
shear_stress_LC_12_2_rear = []
shear_stress_LC_12_3_forward = []
shear_stress_LC_12_3_rear = []
shear_stress_LC_16_1_forward = []
shear_stress_LC_16_1_rear = []
shear_stress_LC_16_2_forward = []
shear_stress_LC_16_2_rear = []
shear_stress_LC_16_3_forward = []
shear_stress_LC_16_3_rear = []
for i in range(len(y_list)):
    #shear_stress_LC_8_1_forward.append(average_shear_stress_8_1_list[i] - shear_stress_torque_LC_8_1_forward[i])
    #shear_stress_LC_8_1_rear.append(average_shear_stress_8_1_list[i] + shear_stress_torque_LC_8_1_rear[i])
    #shear_stress_LC_8_2_forward.append(average_shear_stress_8_2_list[i] - shear_stress_torque_LC_8_2_forward[i])
    #shear_stress_LC_8_2_rear.append(average_shear_stress_8_2_list[i] + shear_stress_torque_LC_8_2_rear[i])
    #shear_stress_LC_8_3_forward.append(average_shear_stress_8_3_list[i] - shear_stress_torque_LC_8_3_forward[i])
    #shear_stress_LC_8_3_rear.append(average_shear_stress_8_3_list[i] + shear_stress_torque_LC_8_3_rear[i])
    shear_stress_LC_12_1_forward.append(average_shear_stress_12_1_list[i] - shear_stress_torque_LC_12_1_forward[i])
    shear_stress_LC_12_1_rear.append(average_shear_stress_12_1_list[i] + shear_stress_torque_LC_12_1_rear[i])
    shear_stress_LC_12_2_forward.append(average_shear_stress_12_2_list[i] - shear_stress_torque_LC_12_2_forward[i])
    shear_stress_LC_12_2_rear.append(average_shear_stress_12_2_list[i] + shear_stress_torque_LC_12_2_rear[i])
    shear_stress_LC_12_3_forward.append(average_shear_stress_12_3_list[i] - shear_stress_torque_LC_12_3_forward[i])
    shear_stress_LC_12_3_rear.append(average_shear_stress_12_3_list[i] + shear_stress_torque_LC_12_3_rear[i])
    shear_stress_LC_16_1_forward.append(average_shear_stress_16_1_list[i] - shear_stress_torque_LC_16_1_forward[i])
    shear_stress_LC_16_1_rear.append(average_shear_stress_16_1_list[i] + shear_stress_torque_LC_16_1_rear[i])
    shear_stress_LC_16_2_forward.append(average_shear_stress_16_2_list[i] - shear_stress_torque_LC_16_2_forward[i])
    shear_stress_LC_16_2_rear.append(average_shear_stress_16_2_list[i] + shear_stress_torque_LC_16_2_rear[i])
    shear_stress_LC_16_3_forward.append(average_shear_stress_16_3_list[i] - shear_stress_torque_LC_16_3_forward[i])
    shear_stress_LC_16_3_rear.append(average_shear_stress_16_3_list[i] + shear_stress_torque_LC_16_3_rear[i])

#PLOTTING OF THE DIFFERENT PHILOSOPHYS AND FORWARD OR REAR

#Desing philosophy 1
#Front
#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
#fig.subplots_adjust(hspace=0.05)
plt.plot(y_list, critical_shear_stress_forward_spar_list_1, '-', label='Tau_Crit')
#plt.plot(y_list, shear_stress_LC_8_1_forward, '-', label='Tau_LC8')
plt.plot(y_list, shear_stress_LC_12_1_forward, '-', label='Tau_LC12')
plt.plot(y_list, shear_stress_LC_16_1_forward, '-', label='Tau_LC16')
plt.title('Shear stress in the forward spar, philosophy 1')
#ax2.plot(y_list, critical_shear_stress_forward_spar_list_1, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_1_forward, '-', label='Tau_LC8')
#ax2.plot(y_list, shear_stress_LC_12_1_forward, '-', label='Tau_LC12')
#ax2.plot(y_list, shear_stress_LC_16_1_forward, '-', label='Tau_LC16')
plt.legend(loc = 'upper right')
plt.xlabel('Span [m]')
plt.ylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_1_forward), min(shear_stress_LC_16_1_forward))
maximum = max(max(shear_stress_LC_12_1_forward), max(shear_stress_LC_16_1_forward))
plt.ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), critical_shear_stress_forward_spar_list_1[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
#ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
#ax1.spines.bottom.set_visible(False)
#ax2.spines.top.set_visible(False)
#ax1.xaxis.tick_top()
#ax1.tick_params(labeltop=False)  # don't put tick labels at the top
#ax2.xaxis.tick_bottom()
#d = .5  # proportion of vertical to horizontal extent of the slanted line
#kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
#ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
#ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_1_forward[i] < critical_shear_stress_forward_spar_list_1[i] and shear_stress_LC_16_1_forward[i] < critical_shear_stress_forward_spar_list_1[i]:
        l+=1
    else:
        print("philosophy 1, forward spar: fail")
if l == len(y_list):
    print("philosophy 1, forward spar: pass")

#Rear
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)
ax1.plot(y_list, critical_shear_stress_rear_spar_list_1, '-', label='Tau_Crit')
#ax1.plot(y_list, shear_stress_LC_8_1_rear, '-', label='Tau_LC8')
ax1.plot(y_list, shear_stress_LC_12_1_rear, '-', label='Tau_LC12')
ax1.plot(y_list, shear_stress_LC_16_1_rear, '-', label='Tau_LC16')
plt.suptitle('Shear stress in the rear spar, philosophy 1')
ax2.plot(y_list, critical_shear_stress_rear_spar_list_1, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_1_rear, '-', label='Tau_LC8')
ax2.plot(y_list, shear_stress_LC_12_1_rear, '-', label='Tau_LC12')
ax2.plot(y_list, shear_stress_LC_16_1_rear, '-', label='Tau_LC16')
ax1.legend(loc = 'upper right')
plt.xlabel('Span [m]')
fig.supylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_1_rear), min(shear_stress_LC_16_1_rear))
maximum = max(max(shear_stress_LC_12_1_rear), max(shear_stress_LC_16_1_rear))
ax1.set_ylim(critical_shear_stress_rear_spar_list_1[0] - 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))), critical_shear_stress_rear_spar_list_1[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
ax1.grid()
ax2.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_1_rear[i] < critical_shear_stress_rear_spar_list_1[i] and shear_stress_LC_16_1_rear[i] < critical_shear_stress_rear_spar_list_1[i]:
        l+=1
    else:
        print("philosophy 1, rear spar: fail")
if l == len(y_list):
    print("philosophy 1, rear spar: pass")

#Desing philosophy 2
#Front
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)
ax1.plot(y_list, critical_shear_stress_forward_spar_list_2, '-', label='Tau_Crit')
#ax1.plot(y_list, shear_stress_LC_8_2_forward, '-', label='Tau_LC8')
ax1.plot(y_list, shear_stress_LC_12_2_forward, '-', label='Tau_LC12')
ax1.plot(y_list, shear_stress_LC_16_2_forward, '-', label='Tau_LC16')
plt.suptitle('Shear stress in the forward spar, philosophy 2')
ax2.plot(y_list, critical_shear_stress_forward_spar_list_2, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_2_forward, '-', label='Tau_LC8')
ax2.plot(y_list, shear_stress_LC_12_2_forward, '-', label='Tau_LC12')
ax2.plot(y_list, shear_stress_LC_16_2_forward, '-', label='Tau_LC16')
ax1.legend(loc = 'upper right')
plt.xlabel('Span [m]')
fig.supylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_2_forward), min(shear_stress_LC_16_2_forward))
maximum = max(max(shear_stress_LC_12_2_forward), max(shear_stress_LC_16_2_forward))
ax1.set_ylim(critical_shear_stress_forward_spar_list_2[0] - 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))), critical_shear_stress_forward_spar_list_2[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
ax1.grid()
ax2.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_2_forward[i] < critical_shear_stress_forward_spar_list_2[i] and shear_stress_LC_16_2_forward[i] < critical_shear_stress_forward_spar_list_2[i]:
        l+=1
    else:
        print("philosophy 2, forward spar: fail")
if l == len(y_list):
    print("philosophy 2, forward spar: pass")

#Rear
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05)
ax1.plot(y_list, critical_shear_stress_rear_spar_list_2, '-', label='Tau_Crit')
#ax1.plot(y_list, shear_stress_LC_8_2_rear, '-', label='Tau_LC8')
ax1.plot(y_list, shear_stress_LC_12_2_rear, '-', label='Tau_LC12')
ax1.plot(y_list, shear_stress_LC_16_2_rear, '-', label='Tau_LC16')
plt.suptitle('Shear stress in the rear spar, philosophy 2')
ax2.plot(y_list, critical_shear_stress_rear_spar_list_2, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_2_rear, '-', label='Tau_LC8')
ax2.plot(y_list, shear_stress_LC_12_2_rear, '-', label='Tau_LC12')
ax2.plot(y_list, shear_stress_LC_16_2_rear, '-', label='Tau_LC16')
ax1.legend(loc = 'upper right')
plt.xlabel('Span [m]')
fig.supylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_2_rear), min(shear_stress_LC_16_2_rear))
maximum = max(max(shear_stress_LC_12_2_rear), max(shear_stress_LC_16_2_rear))
ax1.set_ylim(critical_shear_stress_rear_spar_list_2[0] - 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))), critical_shear_stress_rear_spar_list_2[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.tick_top()
ax1.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
d = .5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
ax1.grid()
ax2.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_1_rear[i] < critical_shear_stress_rear_spar_list_1[i] and shear_stress_LC_16_1_rear[i] < critical_shear_stress_rear_spar_list_1[i]:
        l+=1
    else:
        print("philosophy 2, rear spar: fail")
if l == len(y_list):
    print("philosophy 2, rear spar: pass")

#Desing philosophy 3
#Front
#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
#fig.subplots_adjust(hspace=0.05)
plt.plot(y_list, critical_shear_stress_forward_spar_list_3, '-', label='Tau_Crit')
#plt.plot(y_list, shear_stress_LC_8_3_forward, '-', label='Tau_LC8')
plt.plot(y_list, shear_stress_LC_12_3_forward, '-', label='Tau_LC12')
plt.plot(y_list, shear_stress_LC_16_3_forward, '-', label='Tau_LC16')
plt.title('Shear stress in the forward spar, philosophy 3')
#ax2.plot(y_list, critical_shear_stress_forward_spar_list_3, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_3_forward, '-', label='Tau_LC8')
#ax2.plot(y_list, shear_stress_LC_12_3_forward, '-', label='Tau_LC12')
#ax2.plot(y_list, shear_stress_LC_16_3_forward, '-', label='Tau_LC16')
plt.legend(loc = 'upper right')
plt.xlabel('Span [m]')
plt.ylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_3_forward), min(shear_stress_LC_16_3_forward))
maximum = max(max(shear_stress_LC_12_3_forward), max(shear_stress_LC_16_3_forward))
plt.ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), critical_shear_stress_forward_spar_list_3[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
#ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
#ax1.spines.bottom.set_visible(False)
#ax2.spines.top.set_visible(False)
#ax1.xaxis.tick_top()
#ax1.tick_params(labeltop=False)  # don't put tick labels at the top
#ax2.xaxis.tick_bottom()
#d = .5  # proportion of vertical to horizontal extent of the slanted line
#kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
#ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
#ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_3_forward[i] < critical_shear_stress_forward_spar_list_3[i] and shear_stress_LC_16_3_forward[i] < critical_shear_stress_forward_spar_list_3[i]:
        l+=1
    else:
        print("philosophy 3, forward spar: fail")
if l == len(y_list):
    print("philosophy 3, forward spar: pass")

#Rear
#fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
#fig.subplots_adjust(hspace=0.05)
plt.plot(y_list, critical_shear_stress_rear_spar_list_3, '-', label='Tau_Crit')
#plt.plot(y_list, shear_stress_LC_8_3_rear, '-', label='Tau_LC8')
plt.plot(y_list, shear_stress_LC_12_3_rear, '-', label='Tau_LC12')
plt.plot(y_list, shear_stress_LC_16_3_rear, '-', label='Tau_LC16')
plt.title('Shear stress in the rear spar, philosophy 3')
#ax2.plot(y_list, critical_shear_stress_rear_spar_list_3, '-', label='Tau_Crit')
#ax2.plot(y_list, shear_stress_LC_8_3_rear, '-', label='Tau_LC8')
#ax2.plot(y_list, shear_stress_LC_12_3_rear, '-', label='Tau_LC12')
#ax2.plot(y_list, shear_stress_LC_16_3_rear, '-', label='Tau_LC16')
plt.legend(loc = 'upper right')
plt.xlabel('Span [m]')
plt.ylabel('Shear stress [Pa]')
minimum = min(min(shear_stress_LC_12_3_rear), min(shear_stress_LC_16_3_rear))
maximum = max(max(shear_stress_LC_12_3_rear), max(shear_stress_LC_16_3_rear))
plt.ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), critical_shear_stress_rear_spar_list_3[0] + 0.5 * (maximum - minimum + max(abs(minimum), abs(maximum))))  # outliers only
#ax2.set_ylim(minimum - 0.5 * max(abs(minimum), abs(maximum)), maximum + 0.5 * max(abs(minimum), abs(maximum)))
#ax1.spines.bottom.set_visible(False)
#ax2.spines.top.set_visible(False)
#ax1.xaxis.tick_top()
#ax1.tick_params(labeltop=False)  # don't put tick labels at the top
#ax2.xaxis.tick_bottom()
#d = .5  # proportion of vertical to horizontal extent of the slanted line
#kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
#ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
#ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
plt.grid()
plt.show()
l=0
for i in range(len(y_list)):
    if shear_stress_LC_12_3_rear[i] < critical_shear_stress_rear_spar_list_3[i] and shear_stress_LC_16_3_rear[i] < critical_shear_stress_rear_spar_list_3[i]:
        l+=1
    else:
        print("philosophy 3, rear spar: fail")
if l == len(y_list):
    print("philosophy 3, rear spar: pass")

