import Copy_Variables_for_crossection as Var





def Moi_z_wingbox(p, no_str, fr_t_root, re_t_root, top_sheet_t_root, bottom_sheet_t_root):
        #input parameters

    #Taper ratio and span
    Taper_ratio = Var.Taper_ratio
    Span = Var.Span

    #material properties for AL6061-T6
    e_mod = Var.e_mod #[gpa]
    g_mod = Var.g_mod #[gpa]

    #stringers
    Str_A = Var.Str_A
    Str_N = no_str #the number of stringers has to be 2, 6, 10, 14, 18, 22, etc...

    #spanwise location
    y = p

    #GEOMETRIES AT THE ROOT

    #spars root (and tip for front length)
    Spar_fr_len_root = Var.Spar_fr_len_root
    Spar_fr_len_tip = Var.Spar_fr_len_tip
    Spar_fr_th_root= fr_t_root

    Spar_re_len_root = Var.Spar_re_len_root
    Spar_re_th_root = re_t_root

    #sheets root 
    Sheet_top_len_root = Var.Sheet_top_len_root
    Sheet_top_th_root = top_sheet_t_root

    Sheet_bottom_len_root = Var.Sheet_bottom_len_root
    Sheet_bottom_th_root = bottom_sheet_t_root

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
    total_area = Spar_fr_len*Spar_fr_th + Spar_re_len*Spar_re_th + Sheet_bottom_len*Sheet_bottom_th + Sheet_top_len*Sheet_top_th + Str_A*Str_N
    Str_dis_top = Sheet_top_len / (Str_N/2 + 1)
    Str_dis_bottom = Sheet_bottom_len / (Str_N/2 + 1)

    #Z tilda calculations for every part
    x_tilda_Spar_fr = 0
    x_tilda_Spar_re = Sheet_top_len * math.cos(Sheet_top_angle)
    x_tilda_Sheet_top = (Sheet_top_len/2) * math.cos(Sheet_top_angle)
    x_tilda_Sheet_bottom = (Sheet_bottom_len/2) * math.cos(Sheet_bottom_angle)

    x_tilda_Str_bottom = []
    for i in range(int(Str_N/2)):
        x_tilda_Str_bottom.append( (Str_dis_bottom * (i+1)) * math.cos(Sheet_bottom_angle))

    x_tilda_Str_top = []
    for i in range(int(Str_N/2)):
        x_tilda_Str_top.append( (Str_dis_top * (i+1)* math.cos(Sheet_top_angle)) )

    #Area's of seperate parts
    Spar_fr_A = Spar_fr_len*Spar_fr_th
    Spar_re_A = Spar_re_len*Spar_re_th
    Sheet_top_A = Sheet_top_len*Sheet_top_th
    Sheet_bottom_A = Sheet_bottom_len*Sheet_bottom_th

    #Area multiplied by z tilda for every part
    sum_of_products_Str = 0
    for i in range(int(Str_N/2)):
        sum_of_products_Str = sum_of_products_Str + Str_A * x_tilda_Str_bottom[i]
        sum_of_products_Str = sum_of_products_Str + Str_A * x_tilda_Str_top[i]

    sum_of_products = sum_of_products_Str + Spar_fr_A*x_tilda_Spar_fr + Spar_re_A*x_tilda_Spar_re + Sheet_top_A*x_tilda_Sheet_top + Sheet_bottom_A*x_tilda_Sheet_bottom

    #Calculation of centroid
    Centroid_x = sum_of_products/total_area

    #check:
    #print(x_tilda_Str_bottom)
    #print(x_tilda_Str_top)
    #print(sum_of_products)
    #print(Centroid_x, "yeyyyy")

    #Distances from center of part to z axis: 
    dz_Spar_fr = abs(Centroid_x)
    dz_Spar_re = abs(Centroid_x - x_tilda_Spar_re)
    dz_Sheet_top = abs(Centroid_x - x_tilda_Sheet_top)
    dz_Sheet_bottom = abs(Centroid_x - x_tilda_Sheet_bottom)

    dz_Str_top = []
    for i in range(int(Str_N / 2)):
        dz_Str_top_n = abs(Centroid_x - x_tilda_Str_top[i])
        dz_Str_top.append(dz_Str_top_n)

    dz_Str_bottom = []
    for i in range(int(Str_N / 2)):
        dz_Str_bottom_n = abs(Centroid_x - x_tilda_Str_bottom[i])
        dz_Str_bottom.append(dz_Str_bottom_n)

    #Check:
    #print(dz_Spar_fr)
    #print(dz_Spar_re)
    #print(dz_Sheet_top)
    #print(dz_Sheet_bottom)
    #print(dz_Str_top)
    #print(dz_Str_bottom)

    #I_z' calculations, where "I_z' = twisted formula for moi" if appropriate:
    Iz_pr_Spar_fr = (1/12)*(Spar_fr_len)*((Spar_fr_th)**3)
    Iz_pr_Spar_re = (1/12)*(Spar_re_len)*((Spar_re_th)**3)

    Ix_tw_Sheet_top = (1/12)*(Sheet_top_len)*((Sheet_top_th)**3)
    Iz_tw_Sheet_top = (1/12)*(Sheet_top_th)*((Sheet_top_len)**3)
    Iz_pr_Sheet_top = ((Ix_tw_Sheet_top + Iz_tw_Sheet_top)/2) - (((Ix_tw_Sheet_top - Iz_tw_Sheet_top)/2) * math.cos(2 * Sheet_top_angle))

    Ix_tw_Sheet_bottom = (1/12)*(Sheet_bottom_len)*((Sheet_bottom_th)**3)
    Iz_tw_Sheet_bottom = (1/12)*(Sheet_bottom_th)*((Sheet_bottom_len)**3)
    Iz_pr_Sheet_bottom = ((Ix_tw_Sheet_bottom + Iz_tw_Sheet_bottom)/2) - (((Ix_tw_Sheet_bottom - Iz_tw_Sheet_bottom)/2) * math.cos(2 * Sheet_bottom_angle))

    Iz_pr_Str_top = []
    for i in range(int(Str_N / 2)):
        Iz_pr_Str_top.append(0)

    Iz_pr_Str_bottom = []
    for i in range(int(Str_N / 2)):
        Iz_pr_Str_bottom.append(0)

    #Check:
    #print(Iz_pr_Spar_fr)
    #print(Iz_pr_Spar_re)
    #print(Iz_pr_Sheet_top)
    #print(Iz_pr_Sheet_bottom)
    #print(Iz_pr_Str_top)
    #print(Iz_pr_Str_bottom)

    #Calculation of steiner terms "A*(dz**2)":
    stei_Spar_fr_z = Spar_fr_A * (dz_Spar_fr**2)
    stei_Spar_re_z = Spar_re_A * (dz_Spar_re**2)
    stei_Sheet_top_z = Sheet_top_A * (dz_Sheet_top**2)
    stei_Sheet_bottom_z = Sheet_bottom_A * (dz_Sheet_bottom**2)

    stei_stringer_top_z = []
    for i in range(int(Str_N / 2)):
        stei_stringer_top_n = Str_A * (dz_Str_top[i]**2)
        stei_stringer_top_z.append(stei_stringer_top_n)

    stei_stringer_bottom_z = []
    for i in range(int(Str_N / 2)):
        stei_stringer_bottom_n = Str_A * (dz_Str_bottom[i]**2)
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
    for i in range(int(Str_N / 2)):
        Iz_Str_top_n = Iz_pr_Str_top[i] + stei_stringer_top_z[i]
        Iz_Str_top.append(Iz_Str_top_n)
    Iz_Str_top_tot = sum(Iz_Str_top)

    Iz_Str_bottom = []
    for i in range(int(Str_N / 2)):
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

    return Iz_Wingbox