import math
import numpy as np
import Centroid_x_direction as Cx
import Centroid_z_direction as Cz
import Variables_for_crossection_geometry as Var

#Distances from center of part to z axis: 
dz_Spar_fr = abs(Cx.Centroid_x)
dz_Spar_re = abs(Cx.Centroid_x - Cx.x_tilda_Spar_re)
dz_Sheet_top = abs(Cx.Centroid_x - Cx.x_tilda_Sheet_top)
dz_Sheet_bottom = abs(Cx.Centroid_x - Cx.x_tilda_Sheet_bottom)

dz_Str_top = []
for i in range(int(Var.Str_N / 2)):
    dz_Str_top_n = abs(Cx.Centroid_x - Cx.x_tilda_Str_top[i])
    dz_Str_top.append(dz_Str_top_n)

dz_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    dz_Str_bottom_n = abs(Cx.Centroid_x - Cx.x_tilda_Str_bottom[i])
    dz_Str_bottom.append(dz_Str_bottom_n)

#Check:
#print(dz_Spar_fr)
#print(dz_Spar_re)
#print(dz_Sheet_top)
#print(dz_Sheet_bottom)
#print(dz_Str_top)
#print(dz_Str_bottom)

#I_z' calculations, where "I_z' = twisted formula for moi" if appropriate:
Iz_pr_Spar_fr = (1/12)*(Var.Spar_fr_len)*((Var.Spar_fr_th)**3)
Iz_pr_Spar_re = (1/12)*(Var.Spar_re_len)*((Var.Spar_re_th)**3)

Ix_tw_Sheet_top = (1/12)*(Var.Sheet_top_len)*((Var.Sheet_top_th)**3)
Iz_tw_Sheet_top = (1/12)*(Var.Sheet_top_th)*((Var.Sheet_top_len)**3)
Iz_pr_Sheet_top = ((Ix_tw_Sheet_top + Iz_tw_Sheet_top)/2) - (((Ix_tw_Sheet_top - Iz_tw_Sheet_top)/2) * math.cos(2 * Var.Sheet_top_angle))

Ix_tw_Sheet_bottom = (1/12)*(Var.Sheet_bottom_len)*((Var.Sheet_bottom_th)**3)
Iz_tw_Sheet_bottom = (1/12)*(Var.Sheet_bottom_th)*((Var.Sheet_bottom_len)**3)
Iz_pr_Sheet_bottom = ((Ix_tw_Sheet_bottom + Iz_tw_Sheet_bottom)/2) - (((Ix_tw_Sheet_bottom - Iz_tw_Sheet_bottom)/2) * math.cos(2 * Var.Sheet_bottom_angle))

Iz_pr_Str_top = []
for i in range(int(Var.Str_N / 2)):
    Iz_pr_Str_top.append(0)

Iz_pr_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    Iz_pr_Str_bottom.append(0)

#Check:
#print(Iz_pr_Spar_fr)
#print(Iz_pr_Spar_re)
#print(Iz_pr_Sheet_top)
#print(Iz_pr_Sheet_bottom)
#print(Iz_pr_Str_top)
#print(Iz_pr_Str_bottom)

#Calculation of steiner terms "A*(dz**2)":
stei_Spar_fr_z = Cz.Spar_fr_A * (dz_Spar_fr**2)
stei_Spar_re_z = Cz.Spar_re_A * (dz_Spar_re**2)
stei_Sheet_top_z = Cz.Sheet_top_A * (dz_Sheet_top**2)
stei_Sheet_bottom_z = Cz.Sheet_bottom_A * (dz_Sheet_bottom**2)

stei_stringer_top_z = []
for i in range(int(Var.Str_N / 2)):
    stei_stringer_top_n = Var.Str_A * (dz_Str_top[i]**2)
    stei_stringer_top_z.append(stei_stringer_top_n)

stei_stringer_bottom_z = []
for i in range(int(Var.Str_N / 2)):
    stei_stringer_bottom_n = Var.Str_A * (dz_Str_bottom[i]**2)
    stei_stringer_bottom_z.append(stei_stringer_bottom_n)

#Check:
#print(stei_Spar_fr_z)
#print(stei_Spar_re_z)
#print(stei_Sheet_top_z)
#print(stei_Sheet_bottom_z)
#print(stei_stringer_top_z)
#print(stei_stringer_bottom_Z)

#Moi Iz calculations seperate parts:
Iz_Spar_fr = Iz_pr_Spar_fr + stei_Spar_fr_z
Iz_Spar_re = Iz_pr_Spar_re + stei_Spar_re_z
Iz_Sheet_top = Iz_pr_Sheet_top + stei_Sheet_top_z
Iz_Sheet_bottom = Iz_pr_Sheet_bottom + stei_Sheet_bottom_z

Iz_Str_top = []
for i in range(int(Var.Str_N / 2)):
    Iz_Str_top_n = Iz_pr_Str_top[i] + stei_stringer_top_z[i]
    Iz_Str_top.append(Iz_Str_top_n)
Iz_Str_top_tot = sum(Iz_Str_top)

Iz_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    Iz_Str_bottom_n = Iz_pr_Str_bottom[i] + stei_stringer_bottom_z[i]
    Iz_Str_bottom.append(Iz_Str_bottom_n)
Iz_Str_bottom_tot = sum(Iz_Str_bottom)

#Check:
#print(Iz_Spar_fr_z)
#print(Iz_Spar_re_z)
#print(Iz_Sheet_top_z)
#print(Iz_SHeet_bottom_z)
#print(Iz_Str_top_tot_z)
#print(Iz_Str_bottom_tot_z)

#Total moment of inertia wing box about z axis:
Iz_Wingbox = Iz_Spar_fr + Iz_Spar_re + Iz_Sheet_top + Iz_Sheet_bottom + Iz_Str_top_tot + Iz_Str_bottom_tot
#print("Area moment of inertia of the wingbox about the z-axis is:", Iz_Wingbox, "[m^4]")



