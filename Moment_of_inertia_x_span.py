import matplotlib.pyplot as plt
import math
import numpy as np
import Variables_for_crossection_geometry as Var
from scipy.optimize import curve_fit

#Creating x and y axis:
moment_of_inertia_x_span_1 = []
moment_of_inertia_x_span_2 = []
moment_of_inertia_x_span_3 = []

no_list = [0, 18, 34]
fr_t_list = [0.045, 0.121, 0.100]
re_t_list = [0.045, 0.121, 0.100]
top_t_list = [0.063, 0.056, 0.057]
bottom_t_list = [0.063, 0.056, 0.057]

#Creating the function which calculates the moment of inertia at a certain spanwise location:
def Moi_x_wingbox(p):
    #input parameters

    #Taper ratio and span
    Taper_ratio = Var.Taper_ratio
    Span = Var.Span

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

    #Distances from center of part to x axis: 
    dx_Spar_fr = abs(Centroid_z - z_tilda_Spar_fr)
    dx_Spar_re = abs(Centroid_z - z_tilda_Spar_re)
    dx_Sheet_top = abs(Centroid_z - z_tilda_Sheet_top)
    dx_Sheet_bottom = abs(Centroid_z - z_tilda_Sheet_bottom)

    dx_Str_top = []
    for i in range(int(Str_N / 2)):
        dx_Str_top_n = abs(Centroid_z - z_tilda_Str_top[i])
        dx_Str_top.append(dx_Str_top_n)

    dx_Str_bottom = []
    for i in range(int(Str_N / 2)):
        dx_Str_bottom_n = abs(Centroid_z - z_tilda_Str_bottom[i])
        dx_Str_bottom.append(dx_Str_bottom_n)

    #Check:
    #print(dx_Spar_fr)
    #print(dx_Spar_re)
    #print(dx_Sheet_top)
    #print(dx_Sheet_bottom)
    #print(dx_Str_top)
    #print(dx_Str_bottom)

    #I_x' calculations, where "I_x' = twisted formula for moi" if appropriate:
    Ix_pr_Spar_fr = (1/12)*(Spar_fr_th)*((Spar_fr_len)**3)
    Ix_pr_Spar_re = (1/12)*(Spar_re_th)*((Spar_re_len)**3)

    Ix_tw_Sheet_top = (1/12)*(Sheet_top_len)*((Sheet_top_th)**3)
    Iz_tw_Sheet_top = (1/12)*(Sheet_top_th)*((Sheet_top_len)**3)
    Ix_pr_Sheet_top = ((Ix_tw_Sheet_top + Iz_tw_Sheet_top)/2) + (((Ix_tw_Sheet_top - Iz_tw_Sheet_top)/2) * math.cos(2 * Sheet_top_angle))

    Ix_tw_Sheet_bottom = (1/12)*(Sheet_bottom_len)*((Sheet_bottom_th)**3)
    Iz_tw_Sheet_bottom = (1/12)*(Sheet_bottom_th)*((Sheet_bottom_len)**3)
    Ix_pr_Sheet_bottom = ((Ix_tw_Sheet_bottom + Iz_tw_Sheet_bottom)/2) + (((Ix_tw_Sheet_bottom - Iz_tw_Sheet_bottom)/2) * math.cos(2 * Sheet_bottom_angle))

    Ix_pr_Str_top = []
    for i in range(int(Str_N / 2)):
        Ix_pr_Str_top.append(0)

    Ix_pr_Str_bottom = []
    for i in range(int(Str_N / 2)):
        Ix_pr_Str_bottom.append(0)

    #Check:
    #print(Ix_pr_Spar_fr)
    #print(Ix_pr_Spar_re)
    #print(Ix_pr_Sheet_top)
    #print(Ix_pr_Sheet_bottom)
    #print(Ix_pr_Str_top)
    #print(Ix_pr_Str_bottom)

    #Calculation of steiner terms "A*(dx**2)":
    stei_Spar_fr_x = Spar_fr_A * (dx_Spar_fr**2)
    stei_Spar_re_x = Spar_re_A * (dx_Spar_re**2)
    stei_Sheet_top_x = Sheet_top_A * (dx_Sheet_top**2)
    stei_Sheet_bottom_x = Sheet_bottom_A * (dx_Sheet_bottom**2)

    stei_stringer_top_x = []
    for i in range(int(Str_N / 2)):
        stei_stringer_top_n = Str_A * (dx_Str_top[i]**2)
        stei_stringer_top_x.append(stei_stringer_top_n)

    stei_stringer_bottom_x = []
    for i in range(int(Str_N / 2)):
        stei_stringer_bottom_n = Str_A * (dx_Str_bottom[i]**2)
        stei_stringer_bottom_x.append(stei_stringer_bottom_n)

    #Check:
    #print(stei_Spar_fr_x)
    #print(stei_Spar_re_x)
    #print(stei_Sheet_top_x)
    #print(stei_Sheet_bottom_x)
    #print(stei_stringer_top_x)
    #print(stei_stringer_bottom_x)

    #Moi Ix calculations seperate parts:
    Ix_Spar_fr = Ix_pr_Spar_fr + stei_Spar_fr_x
    Ix_Spar_re = Ix_pr_Spar_re + stei_Spar_re_x
    Ix_Sheet_top = Ix_pr_Sheet_top + stei_Sheet_top_x
    Ix_Sheet_bottom = Ix_pr_Sheet_bottom + stei_Sheet_bottom_x

    Ix_Str_top = []
    for i in range(int(Str_N / 2)):
        Ix_Str_top_n = Ix_pr_Str_top[i] + stei_stringer_top_x[i]
        Ix_Str_top.append(Ix_Str_top_n)
    Ix_Str_top_tot = sum(Ix_Str_top)

    Ix_Str_bottom = []
    for i in range(int(Str_N / 2)):
        Ix_Str_bottom_n = Ix_pr_Str_bottom[i] + stei_stringer_bottom_x[i]
        Ix_Str_bottom.append(Ix_Str_bottom_n)
    Ix_Str_bottom_tot = sum(Ix_Str_bottom)

    #Check:
    #print(Ix_Spar_fr)
    #print(Ix_Spar_re)
    #print(Ix_Sheet_top)
    #print(Ix_SHeet_bottom)
    #print(Ix_Str_top_tot)
    #print(Ix_Str_bottom_tot)

    #Total moment of inertia wing box about x axis:
    Ix_Wingbox = Ix_Spar_fr + Ix_Spar_re + Ix_Sheet_top + Ix_Sheet_bottom + Ix_Str_top_tot + Ix_Str_bottom_tot
    #print("Area moment of inertia of the wingbox about the x-axis is:", Ix_Wingbox, "[m^4]")

    return Ix_Wingbox

#Check:
#print(Moi_x_wingbox(0))

#Calculating all the different values along the span for both x and y values
Span_y_x = np.arange(0.0, (Var.Span / 2), 1.0)

#Check:
#print(Span_y_1)
index = 0
for i in range(0, len(Span_y_x)):
    Moi_x_wingbox_y = Moi_x_wingbox(Span_y_x[i])
    moment_of_inertia_x_span_1.append(Moi_x_wingbox_y)
index += 1
for i in range(0, len(Span_y_x)):
    Moi_x_wingbox_y = Moi_x_wingbox(Span_y_x[i])
    moment_of_inertia_x_span_2.append(Moi_x_wingbox_y)
index += 1
for i in range(0, len(Span_y_x)):
    Moi_x_wingbox_y = Moi_x_wingbox(Span_y_x[i])
    moment_of_inertia_x_span_3.append(Moi_x_wingbox_y)

#Check:
#print(moment_of_inertia_x_span)

#Set up for test function:
def test_function(x, A, B, C, D):
    y = A*(x**3) + B*(x**2) + C*x + D 
    return y 

Parameters1, covariance = curve_fit(test_function, Span_y_x, moment_of_inertia_x_span_1)

Fit_A_1 = Parameters1[0]
Fit_B_1 = Parameters1[1]
Fit_C_1 = Parameters1[2]
Fit_D_1 = Parameters1[3]

def print_fit():
    return Parameters1

Fit_y_1 = test_function(Span_y_x, Fit_A_1, Fit_B_1, Fit_C_1, Fit_D_1)
#print('The values for A, B and C in the function Ax^3 + Bx^2 + Cx + D are:', Fit_A, Fit_B, Fit_C, Fit_D)

#Plotting the results 
# plt.plot(Span_y_x, moment_of_inertia_x_span_1, 'o', label='Data')
plt.plot(Span_y_x, Fit_y_1, '-', label='Fit')
plt.xlabel("Spanwise location")
plt.ylabel("Moment of inertia about x axis")
plt.show()











