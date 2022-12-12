import math
import numpy as np
import Centroid_x_direction as Cx
import Centroid_z_direction as Cz
import Variables_for_crossection_geometry as Var

#Distances from center of part to x axis: 
dx_Spar_fr = abs(Cz.Centroid_z - Cz.z_tilda_Spar_fr)
dx_Spar_re = abs(Cz.Centroid_z - Cz.z_tilda_Spar_re)
dx_Sheet_top = abs(Cz.Centroid_z - Cz.z_tilda_Sheet_top)
dx_Sheet_bottom = abs(Cz.Centroid_z - Cz.z_tilda_Sheet_bottom)

dx_Str_top = []
for i in range(int(Var.Str_N / 2)):
    dx_Str_top_n = abs(Cz.Centroid_z - Cz.z_tilda_Str_top[i])
    dx_Str_top.append(dx_Str_top_n)

dx_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    dx_Str_bottom_n = abs(Cz.Centroid_z - Cz.z_tilda_Str_bottom[i])
    dx_Str_bottom.append(dx_Str_bottom_n)

#Check:
#print(dx_Spar_fr)
#print(dx_Spar_re)
#print(dx_Sheet_top)
#print(dx_Sheet_bottom)
#print(dx_Str_top)
#print(dx_Str_bottom)

#I_x' calculations, where "I_x' = twisted formula for moi" if appropriate:
Ix_pr_Spar_fr = (1/12)*(Var.Spar_fr_th)*((Var.Spar_fr_len)**3)
Ix_pr_Spar_re = (1/12)*(Var.Spar_re_th)*((Var.Spar_re_len)**3)

Ix_tw_Sheet_top = (1/12)*(Var.Sheet_top_len)*((Var.Sheet_top_th)**3)
Iz_tw_Sheet_top = (1/12)*(Var.Sheet_top_th)*((Var.Sheet_top_len)**3)
Ix_pr_Sheet_top = ((Ix_tw_Sheet_top + Iz_tw_Sheet_top)/2) + (((Ix_tw_Sheet_top - Iz_tw_Sheet_top)/2) * math.cos(2 * Var.Sheet_top_angle))

Ix_tw_Sheet_bottom = (1/12)*(Var.Sheet_bottom_len)*((Var.Sheet_bottom_th)**3)
Iz_tw_Sheet_bottom = (1/12)*(Var.Sheet_bottom_th)*((Var.Sheet_bottom_len)**3)
Ix_pr_Sheet_bottom = ((Ix_tw_Sheet_bottom + Iz_tw_Sheet_bottom)/2) + (((Ix_tw_Sheet_bottom - Iz_tw_Sheet_bottom)/2) * math.cos(2 * Var.Sheet_bottom_angle))

Ix_pr_Str_top = []
for i in range(int(Var.Str_N / 2)):
    Ix_pr_Str_top.append(0)

Ix_pr_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    Ix_pr_Str_bottom.append(0)

#Check:
#print(Ix_pr_Spar_fr)
#print(Ix_pr_Spar_re)
#print(Ix_pr_Sheet_top)
#print(Ix_pr_Sheet_bottom)
#print(Ix_pr_Str_top)
#print(Ix_pr_Str_bottom)

#Calculation of steiner terms "A*(dx**2)":
stei_Spar_fr_x = Cx.Spar_fr_A * (dx_Spar_fr**2)
stei_Spar_re_x = Cx.Spar_re_A * (dx_Spar_re**2)
stei_Sheet_top_x = Cx.Sheet_top_A * (dx_Sheet_top**2)
stei_Sheet_bottom_x = Cx.Sheet_bottom_A * (dx_Sheet_bottom**2)

stei_stringer_top_x = []
for i in range(int(Var.Str_N / 2)):
    stei_stringer_top_n = Var.Str_A * (dx_Str_top[i]**2)
    stei_stringer_top_x.append(stei_stringer_top_n)

stei_stringer_bottom_x = []
for i in range(int(Var.Str_N / 2)):
    stei_stringer_bottom_n = Var.Str_A * (dx_Str_bottom[i]**2)
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
for i in range(int(Var.Str_N / 2)):
    Ix_Str_top_n = Ix_pr_Str_top[i] + stei_stringer_top_x[i]
    Ix_Str_top.append(Ix_Str_top_n)
Ix_Str_top_tot = sum(Ix_Str_top)

Ix_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
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