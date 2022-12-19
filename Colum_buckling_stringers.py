#also check compressive failure for this component

import Variables_for_crossection_geometry as Var
from math import *



#Moment of inertia calculation X axis for stringer
#centroid z direction
Str_z_tilda_= (Var.Str_v_len*Var.Str_v_th*Var.Str_v_len/2)/Var.Str_A
Str_MOI_X = (Var.Str_h_len*Var.Str_h_th**3/12 + Var.Str_h_len*Var.Str_h_th*(Str_z_tilda_ - Var.Str_h_th/2)**2) + (Var.Str_v_len**3*Var.Str_v_th/12 + Var.Str_v_len*Var.Str_v_th * (Var.Str_v_len/2-Str_z_tilda_)**2) 

def critical_stress_b_Str(L):
    sigma = (Var.k * pi**2 *Var.e_mod *Str_MOI_X)/(L**2*Var.Str_A)
    return sigma

print(critical_stress_b_Str(Var.Span/(2*40))/1e6, "MPa")             #-67 MPa

ribs_spacing_list = [0.534261602642555, 0.525745792613091, 0.517948192133515, 0.510807981088931, 0.504278027391229, 0.498322565969259, 0.492915713953116, 0.488040637415200, 0.483689285827142, 0.479862693176665, 0.476571931811814, 0.473839921959174, 0.471704487447047, 0.470223387542136, 0.469482727462029, 0.469611605477857, 0.470809336470809, 0.707830742064340, 0.719502430751495, 1.23647866009611, 0.644328916142815, 0.666243412755113, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710, 0.903329851150710]