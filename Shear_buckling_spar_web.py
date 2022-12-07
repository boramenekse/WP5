import math
#also check compressive failure for this component

### SETTING SOME CONSTANTS VOR THE FORMULAE

pi = math.pi            # get pi from the math repository
E = 69 * (10 ** 6)      # define the young's modulus, in Pa
v = 0.33                # define the poissons ratio, no unit
k_s = 9.6               # define the k_s constant, found in figure 16 from the reader, no unit

def get_chord(y):
    c = 11.59 - y * (11.59-3.25) / 33.39
    return c

def get_forward_spar_lenth(y):
    l_fs = (0.045100 + 0.045200) * get_chord(y)
    return l_fs

def get_rear_spar_lenth(y):
    l_rs = (0.033900 + 0.016200) * get_chord(y)
    return l_rs

def get_thickness_forward_spar_1(y):
    t_fs_1 =
    return t_fs_1

def get_thickness_rear_spar_1(y):
    t_rs_1 =
    return t_rs_1

def get_thickness_forward_spar_2(y):
    t_fs_2 =
    return t_fs_2

def get_thickness_rear_spar_2(y):
    t_fr_2 =
    return t_fr_2

def get_thickness_forward_spar_3(y):
    t_fs_3 =
    return t_fs_3

def get_thickness_rear_spar_3(y):
    t_fr_3 =
    return t_fr_3
