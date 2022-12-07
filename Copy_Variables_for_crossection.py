#copy of the variable list used for the previous workpackage

import math

#input parameters

#Taper ratio and span and chord length at the root
Taper_ratio = 0.279881175
Span = 66.77675839 #[m]
Chord_root = 11.59426354 #[m]

#material properties for AL6061-T6
e_mod = 69*10**9 #[pa]
g_mod = 26*10**9 #[pa]

#L stringer parameters

Str_h_len = 1
Str_v_len = 1
Str_h_th = 0.1
Str_v_th = 0.1

#Str_A = 0.001 #[m^2]
Str_A = Str_h_len*Str_h_th + Str_v_len*Str_v_th
Str_N = 0 #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...




#spanwise location
y = 0

#GEOMETRIES AT THE ROOT

#spars root (and tip for front length) #spar is located at x=0.2c
Spar_fr_len_root = (0.045100 + 0.045200)*Chord_root #data from https://m-selig.ae.illinois.edu/ads/coord/sc20710.dat
Spar_fr_len_tip = Spar_fr_len_root * Taper_ratio
Spar_fr_th_root= 0.045

#spar rear is at x=0.75c
Spar_re_len_root = (0.033900 + 0.016200)*Chord_root #data from https://m-selig.ae.illinois.edu/ads/coord/sc20710.dat
Spar_re_th_root = 0.045

#sheets root 
Sheet_top_len_root = math.sqrt((0.045100 - 0.033900)**2 + 0.55**2)*Chord_root #simply using some Pythagoras.
Sheet_top_th_root = 0.063 # Need more or less 0.063 to meet the requirement

Sheet_bottom_len_root = math.sqrt((0.045200 - 0.016200)**2 + 0.55**2)*Chord_root  #simply using some Pythagoras.
Sheet_bottom_th_root = 0.063 # Need more or less 0.063 to meet the requirement

Sheet_top_angle = math.atan((0.045100-0.033900)/(0.75-0.2))
Sheet_bottom_angle = math.atan((0.045200-0.016200)/(0.75-0.2))

#Ratio's
Constant = (Spar_fr_len_tip - Spar_fr_len_root)/(Span/2)

Ratio_Sheet_top_len = Spar_fr_len_root / Sheet_top_len_root
Ratio_Sheet_bottom_len = Spar_fr_len_root / Sheet_bottom_len_root
Ratio_Spar_re_len = Spar_fr_len_root / Spar_re_len_root
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


