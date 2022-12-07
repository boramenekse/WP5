#also check compressive failure for this component

import Copy_Variables_for_crossection as Var
from math import *



#Moment of inertia calculation X axis for stringer
#centroid z direction
Str_z_tilda_= (Var.Str_v_len*Var.Str_v_th*Var.Str_v_len/2)/Var.Str_A
Str_MOI_X = (Var.Str_h_len*Var.Str_h_th**3/12 + Var.Str_h_len*Var.Str_h_th*(Str_z_tilda_ - Var.Str_h_th/2)**2) + (Var.Str_v_len**3*Var.Str_v_th/12 + Var.Str_v_len*Var.Str_v_th * (Var.Str_v_len/2-Str_z_tilda_)**2) 

def stress_b_Str():
    sigma = (Var.k * pi)/(1)
    return 

print(pi)