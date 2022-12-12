import math
import Variables_for_crossection_geometry as Var

#Total area and distance to stringer from start of sheet
total_area = Var.Spar_fr_len * Var.Spar_fr_th + Var.Spar_re_len * Var.Spar_re_th + Var.Sheet_bottom_len * Var.Sheet_bottom_th + Var.Sheet_top_len * Var.Sheet_top_th + Var.Str_A * Var.Str_N
Str_dis_top = Var.Sheet_top_len / ((Var.Str_N / 2) + 1)
Str_dis_bottom = Var.Sheet_bottom_len / ((Var.Str_N / 2) + 1)

#Z tilda calculations for every part
z_tilda_Spar_fr = Var.Spar_fr_len / 2
z_tilda_Spar_re = (Var.Sheet_top_len) * math.sin(Var.Sheet_top_angle) + Var.Spar_re_len / 2
z_tilda_Sheet_top = (Var.Sheet_top_len / 2) * math.sin(Var.Sheet_top_angle)
z_tilda_Sheet_bottom = (Var.Spar_fr_len) - ((Var.Sheet_bottom_len / 2) * math.sin(Var.Sheet_bottom_angle))

z_tilda_Str_bottom = []
for i in range(int(Var.Str_N / 2)):
    z_tilda_Str_N = Var.Spar_fr_len - (Str_dis_bottom * (i+1) * math.sin(Var.Sheet_bottom_angle))
    z_tilda_Str_bottom.append(z_tilda_Str_N)

z_tilda_Str_top = []
for i in range(int(Var.Str_N / 2)):
    z_tilda_Str_P = Str_dis_top * (i+1) * math.sin(Var.Sheet_top_angle)
    z_tilda_Str_top.append(z_tilda_Str_P)

#Area's of seperate parts
Spar_fr_A = Var.Spar_fr_len*Var.Spar_fr_th
Spar_re_A = Var.Spar_re_len*Var.Spar_re_th
Sheet_top_A = Var.Sheet_top_len*Var.Sheet_top_th
Sheet_bottom_A = Var.Sheet_bottom_len*Var.Sheet_bottom_th

#Area multiplied by z tilda for every part
sum_of_products_Str_top = 0
for i in range(int(Var.Str_N / 2)):
    sum_of_products_Str_top += z_tilda_Str_top[i] * Var.Str_A

sum_of_products_Str_bottom = 0
for i in range(int(Var.Str_N / 2)):
    sum_of_products_Str_bottom += z_tilda_Str_bottom[i] * Var.Str_A

sum_of_products_Spar = z_tilda_Spar_fr * Spar_fr_A + z_tilda_Spar_re * Spar_re_A

sum_of_products_Sheet = z_tilda_Sheet_top * Sheet_top_A + z_tilda_Sheet_bottom * Sheet_bottom_A

sum_of_products = sum_of_products_Sheet + sum_of_products_Spar + sum_of_products_Str_bottom + sum_of_products_Str_top

#Calculation of centroid
Centroid_z = sum_of_products / total_area



