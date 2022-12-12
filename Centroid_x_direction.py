import Variables_for_crossection_geometry  as Var
from math import cos

#Total area and distance to stringer from start of sheet
total_area = Var.Spar_fr_len*Var.Spar_fr_th + Var.Spar_re_len*Var.Spar_re_th + Var.Sheet_bottom_len*Var.Sheet_bottom_th + Var.Sheet_top_len*Var.Sheet_top_th + Var.Str_A*Var.Str_N
Str_dis_top = Var.Sheet_top_len / (Var.Str_N/2 + 1)
Str_dis_bottom = Var.Sheet_bottom_len / (Var.Str_N/2 + 1)

#Z tilda calculations for every part
x_tilda_Spar_fr = 0
x_tilda_Spar_re = Var.Sheet_top_len * cos(Var.Sheet_top_angle)
x_tilda_Sheet_top = (Var.Sheet_top_len/2) * cos(Var.Sheet_top_angle)
x_tilda_Sheet_bottom = (Var.Sheet_bottom_len/2) * cos(Var.Sheet_bottom_angle)

x_tilda_Str_bottom = []
for i in range(int(Var.Str_N/2)):
    x_tilda_Str_bottom.append( (Str_dis_bottom * (i+1)) * cos(Var.Sheet_bottom_angle))

x_tilda_Str_top = []
for i in range(int(Var.Str_N/2)):
    x_tilda_Str_top.append( (Str_dis_top * (i+1)* cos(Var.Sheet_top_angle)) )

#Area's of seperate parts
Spar_fr_A = Var.Spar_fr_len*Var.Spar_fr_th
Spar_re_A = Var.Spar_re_len*Var.Spar_re_th
Sheet_top_A = Var.Sheet_top_len*Var.Sheet_top_th
Sheet_bottom_A = Var.Sheet_bottom_len*Var.Sheet_bottom_th

#Area multiplied by z tilda for every part
sum_of_products_Str = 0
for i in range(int(Var.Str_N/2)):
    sum_of_products_Str = sum_of_products_Str + Var.Str_A * x_tilda_Str_bottom[i]
    sum_of_products_Str = sum_of_products_Str + Var.Str_A * x_tilda_Str_top[i]

sum_of_products = sum_of_products_Str + Spar_fr_A*x_tilda_Spar_fr + Spar_re_A*x_tilda_Spar_re + Sheet_top_A*x_tilda_Sheet_top + Sheet_bottom_A*x_tilda_Sheet_bottom

#Calculation of centroid
Centroid_x = sum_of_products/total_area

#check:
#print(x_tilda_Str_bottom)
#print(x_tilda_Str_top)
#print(sum_of_products)
#print(Centroid_x, "yeyyyy")