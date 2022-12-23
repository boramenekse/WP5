import math
#from wp4_1.xflr5 import M_distr_list

#input parameters

#Taper ratio and span and chord length at the root
Taper_ratio = 0.279881175
Span = 66.77675839 #[m]
Chord_root = 11.59426354 #[m]

#material properties for AL6061-T6
e_mod = 69*10**9 #[pa]
g_mod = 26*10**9 #[pa]

#stringers
Str_A = 0.001 #[m^2]
Str_N = 34 #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...

#spanwise location
y = 0

#GEOMETRIES AT THE ROOT

#spars root (and tip for front length) #spar is located at x=0.2c
Spar_fr_len_root = (0.045100 + 0.045200)*Chord_root #data from https://m-selig.ae.illinois.edu/ads/coord/sc20710.dat
Spar_fr_len_tip = Spar_fr_len_root * Taper_ratio
Spar_fr_th_root= 0.015

#spar rear is at x=0.75c
Spar_re_len_root = (0.033900 + 0.016200)*Chord_root #data from https://m-selig.ae.illinois.edu/ads/coord/sc20710.dat
Spar_re_th_root = 0.015

#sheets root 
Sheet_top_len_root = math.sqrt((0.045100 - 0.033900)**2 + 0.55**2)*Chord_root #simply using some Pythagoras.
Sheet_top_th_root = 0.013 # Need more or less 0.063 to meet the requirement

Sheet_bottom_len_root = math.sqrt((0.045200 - 0.016200)**2 + 0.55**2)*Chord_root  #simply using some Pythagoras.
Sheet_bottom_th_root = 0.013 # Need more or less 0.063 to meet the requirement

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


#Check:
#print(Spar_fr_len)
#print(Spar_fr_th)
#print(Spar_re_len)
#print(Spar_re_th)
#print(Sheet_top_len)
#print(Sheet_top_th)
#print(Sheet_bottom_len)
#print(Sheet_bottom_th)


#Function for moment

#Moment_distribution_constants = [-0.09613131517333155, 5.6005659862894595, -227.62613584140453, 24235.558185874754, -1093618.4267466187, 14997899.497052994]  #LC 1
#Moment_distribution_constants = [-0.09619501984986228, 5.602851321582722, -227.6887446894729, 24254.589352459083, -1094624.2999532684, 15012401.913543124] #LC 2
#Moment_distribution_constants = [-0.17102654074339463, 11.266598408064446, -642.7225762637838, 59457.406498209646, -2435602.742200649, 32056283.381364055] #LC 3
#Moment_distribution_constants = [-0.1710367457937722, 11.266968241633654, -642.732699863465, 59460.432103918596, -2435762.578457649, 32058587.79909822] #LC 4
#Moment_distribution_constants = [0.19889633611036153, -7.8835636493059225, 132.34764504574704, -45011.27746777236, 2500969.825372003, -36707168.522779234] #LC 5
#Moment_distribution_constants = [0.19882762656682307, -7.879646173703972, 132.24544996799534, -45000.338448958115, 2500424.8199410927, -36699379.969619706] #LC 6
#Moment_distribution_constants = [0.3527726381564994, -12.641417719819973, -4.3652103102061535, -63701.46469707234, 3958499.772459415, -59919672.9227723] #LC 7
Moment_distribution_constants_8 = [0.3528646728626992, -6.637709292337059, -1119.527990931244, 13911.399064342344, 1610817.1283304724, -34013528.36512487] #LC 8
#Moment_distribution_constants = [-0.114076587665573, 6.351119188973497, -246.9625831161021, 28609.600548619503, -1321275.4667740301, 18254467.348068513] #LC 9
#Moment_distribution_constants = [-0.1140870375866153, 6.351812169868608, -246.97908746198178, 28610.171555441848, -1321296.6397836932, 18254724.439961758] #LC 10
#Moment_distribution_constants = [-0.19568452667000605, 12.142188968395022, -666.2415582285619, 66716.72227515523, -2818848.8123028614, 37570080.41352775] #LC 11
Moment_distribution_constants_12 = [-0.19569866608662626, 6.680180498386035, 347.65040326524246, -3851.039203234463, -684065.9293064886, 14011505.048648307] #LC 12
#Moment_distribution_constants = [0.2285183814038126, -9.627584515401974, 178.12005938690191, -49496.032616846955, 2722962.572266368, -39885043.86687112] #LC 13
Moment_distribution_constants_14 = [0.2303774463057893, -8.32782475226923, -78.97409340148076, -31866.982573840316, 2200496.6737145814, -34216483.98319514] #LC 14
#Moment_distribution_constants = [0.3725902160334711, -13.778330812362093, 25.29532125913686, -66816.08003782954, 4113416.239124321, -62133352.62347586] #LC 15
Moment_distribution_constants_16 = [0.40016072098555516, -9.397942238637903, -1047.3842177269048, 6708.870613212022, 1967575.6173070923, -39113307.35963459] #LC 16

def print_load_cases():
  return Moment_distribution_constants_8, Moment_distribution_constants_12, Moment_distribution_constants_16

TnB_A = Moment_distribution_constants_16[0]
TnB_B = Moment_distribution_constants_16[1]
TnB_C = Moment_distribution_constants_16[2]
TnB_D = Moment_distribution_constants_16[3]
TnB_E = Moment_distribution_constants_16[4]
TnB_F = Moment_distribution_constants_16[5]

#Function for torque


#maximum deflection and twist
max_defl = 0
max_twist = 0

#L stringer parameters

Str_h_len = 0.04
Str_v_len = 0.07
Str_h_th = 0.01
Str_v_th = 0.01

#Str_A = 0.001 #[m^2]
Str_A = Str_h_len*Str_h_th + Str_v_len*Str_v_th         #should be 10 cm^2 from last report
#Str_N = 0 #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...
k = 4   #clamping conditions
L = 1   #spacing between the ribs