#also check compressive failure for this component
import sympy as smp
import Variables_for_crossection_geometry as var
import numpy as np
import matplotlib.pyplot as plt
import Normal_stress_distribution as str 

sigma_cr = smp.symbols('\u03C3_{cr}', real=True)
theta = smp.symbols('\u03B8', real=True)
kc = smp.symbols('k_c', real=True)
emod = smp.symbols('E', real=True, positive=True)
gmod = smp.symbols('G', real=True, positive=True)
v = smp.symbols('v', real=True, positive=True)
t = smp.symbols('t', real=True, positive=True)
ttopr = smp.symbols('t_{top_r}', real=True, positive=True)
tbotr = smp.symbols('t_{bottom_r}', real=True, positive=True)
b = smp.symbols('b', real=True, positive=True)
b_2 = smp.symbols('b/2', real=True, positive=True)
x = smp.symbols('x', real=True)
y = smp.symbols('y', real=True)
z = smp.symbols('z', real=True)
a = smp.symbols('a', real=True) 
h = smp.symbols('h', real=True)
d = smp.symbols('d', real=True, positive=True) 
cr = smp.symbols('c_r', real=True, positive=True) 
taper = smp.symbols('\u03BB', real=True, positive=True) 
ixx = smp.symbols('Ixx', real=True)
iyy = smp.symbols('Iyy', real=True)
ixy = smp.symbols('Ixy', real=True)
mx = smp.symbols('Mx', real=True)
my = smp.symbols('My', real=True)
r = smp.symbols('r', real=True, positive=True)
alpha = smp.symbols('\u03B1', real=True)
beta = smp.symbols('\u03B2', real=True)
rho = smp.symbols('\u03C1', real=True)
sigma = smp.symbols('\u03C3', real=True)
eta = smp.symbols('\u03B7', real=True)
xi = smp.symbols('\u03BE', real=True)

# nd = str.stress_16_1.subs([(z, str.br_z1), (x, str.br_x1)]).simplify()
# print(nd)
span = np.linspace(0, 0.5*var.Span, 1000, endpoint=True)
heaviside = smp.Heaviside(y-11.69, 1)
nd = (5140*(-100*y + (100*y - 1169)*heaviside + 1169)*(5.42626672758448e-7*y**3 - 1.88694820902376e-5*y**2 - 0.071424801946438*y + 3.28756683861378)*(7.78733677625702e-7*y**3 - 7.77370531688468e-5*y**2 + 0.0027848684159305*y - 0.0358283831164184) + (9.75121500397927e-9*y**3 - 3.39092024244981e-7*y**2 + 0.00510357950779944*y - 0.237060165638359)*(2.08836637456291e-5*y**3 - 0.00206034315896725*y**2 + 0.0724190974615689*y - 0.908260313093263)*(0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y + (0.400160720985555*y**5 - 9.3979422386379*y**4 - 1047.3842177269*y**3 + 6708.87061321202*y**2 + 1967575.61730709*y - 39113307.3596346)*heaviside - (0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y - 33488163.0896061)*heaviside - 33488163.0896061))/((7.78733677625702e-7*y**3 - 7.77370531688468e-5*y**2 + 0.0027848684159305*y - 0.0358283831164184)*(2.08836637456291e-5*y**3 - 0.00206034315896725*y**2 + 0.0724190974615689*y - 0.908260313093263))
# nd = (5140*(-100*y + (100*y - 1169)*heaviside + 1169)*(5.42626672758448e-7*y**3 - 1.88694820902376e-5*y**2 - 0.071424801946438*y + 3.28756683861378)*(7.78733870357834e-7*y**3 - 7.50415077493545e-5*y**2 + 0.00253491030860201*y - 0.0300337363084974) + (9.75121500397927e-9*y**3 - 3.39092024244981e-7*y**2 + 0.00510357950779944*y - 0.237060165638359)*(2.08842547213758e-5*y**3 - 0.00201247951941106*y**2 + 0.0679817775327849*y - 0.805451291334324)*(0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y + (0.400160720985555*y**5 - 9.3979422386379*y**4 - 1047.3842177269*y**3 + 6708.87061321202*y**2 + 1967575.61730709*y - 39113307.3596346)*heaviside - (0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y - 33488163.0896061)*heaviside - 33488163.0896061))/((7.78733870357834e-7*y**3 - 7.50415077493545e-5*y**2 + 0.00253491030860201*y - 0.0300337363084974)*(2.08842547213758e-5*y**3 - 0.00201247951941106*y**2 + 0.0679817775327849*y - 0.805451291334324))
# nd = (5140*(-100*y + (100*y - 1169)*heaviside + 1169)*(5.3229971557073e-19*y**3 - 2.35518596836698e-17*y**2 - 0.0718773721330028*y + 3.33260133328251)*(7.96371363901842e-9*y**4 - 1.47811910673145e-6*y**3 + 0.000105288889306002*y**2 - 0.0034060384975003*y + 0.0421093825303379) + (2.71678538500068e-20*y**3 - 2.02637753771834e-18*y**2 - 0.00509664457596961*y + 0.236306420297532)*(2.08841726980762e-5*y**3 - 0.0020254948757701*y**2 + 0.0691888314097523*y - 0.83342866611924)*(0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y + (0.400160720985555*y**5 - 9.3979422386379*y**4 - 1047.3842177269*y**3 + 6708.87061321202*y**2 + 1967575.61730709*y - 39113307.3596346)*heaviside - (0.400160720998562*y**5 - 9.39794223895798*y**4 - 1047.38421772391*y**3 + 6708.87061319162*y**2 + 1582191.86730725*y - 33488163.0896061)*heaviside - 33488163.0896061))/((2.08841726980762e-5*y**3 - 0.0020254948757701*y**2 + 0.0691888314097523*y - 0.83342866611924)*(7.96371363901842e-9*y**4 - 1.47811910673145e-6*y**3 + 0.000105288889306002*y**2 - 0.0034060384975003*y + 0.0421093825303379))
# plt.figure()
# des_list = []
# for i in range(0, 1000):
#   des_list.append(-3.0e8)
# plt.plot(span, smp.lambdify([y], nd)(span[0:]))
# plt.plot(span, 1.5*smp.lambdify([y], nd)(span[0:]))
# plt.plot(span, des_list)
# plt.show()


emod = var.e_mod
gmod = var.g_mod
b_2 = smp.Rational(1, 2)*var.Span
cr = var.Chord_root
g_e_relation = smp.Eq(gmod, smp.Rational(1, 2)*(emod/(1+v))) 
v = smp.solve(g_e_relation, v)[0]
print(v.evalf())
taper = var.Taper_ratio  #Removed smp.nsimplify(round())
ttopr = var.Sheet_top_th_root #Same applies here
tbotr = var.Sheet_bottom_th_root #E
ttopfun = ttopr*(1 + (taper-1)*(y/b_2))
tbotfun = tbotr*(1 + (taper-1)*(y/b_2))
theta_top = var.Sheet_top_angle
theta_bot = var.Sheet_bottom_angle
str_no = var.Str_N

def sigma(kc, t, b):
  expr = smp.Rational(1, 12) * ((smp.pi**2 * kc * emod)/(1 - v**2)) * (t/b)**2
  return expr.simplify()

def b(y, theta):
  c_local = cr - cr*(1-taper)*(y/(b_2))
  '''Wingbox starts at 20 percent of the local chord and end at 75 percent of it
     So b is 55 percent of the local chord length -> 11/20
     Sheets are at an angle
  '''
  b_y = c_local*smp.Rational(11, 20)/((0.5*str_no+1)*smp.cos(theta))
  return b_y.simplify()

btopfun = b(y, theta_top)
bbotfun = b(y, theta_bot)


# KS for C with SS loaded edges 
a_b = np.linspace(0, 5, 1000, endpoint=True)
ab = smp.symbols('a_b', real=True, positive=True)
m = smp.symbols('m', real=True, positive=True)
def k(m, a_b):
  expr = (1/a_b)**2*(m+(1/m)*a_b**2)**2
  return expr

m1 = k(1, ab)
m2 = k(2, ab)
m3 = k(3, ab)
m4 = k(4, ab)
m5 = k(5, ab)
m6 = k(6, ab)
kfun = smp.Piecewise((m1, (ab>=0) & (ab <= float(smp.solve(m1-m2, ab)[0]))), (m2, (ab > float(smp.solve(m1-m2, ab)[0])) & (ab <= float(smp.solve(m2-m3, ab)[0]))), (m3, (ab > float(smp.solve(m2-m3, ab)[0])) & (ab <= float(smp.solve(m3-m4, ab)[0]))), (m4, (ab > float(smp.solve(m3-m4, ab)[0])) & (ab <= float(smp.solve(m4-m5, ab)[0]))), (m5, (ab > float(smp.solve(m4-m5, ab)[0])) & (ab <= 5)), (0, True))

sfc = 1.5
nd_des = sfc*nd
sigmatopfun = sigma(kc, ttopfun, btopfun)
sigmabotfun = sigma(kc, tbotfun, bbotfun)

# print(-nd_des.subs(y, 0))
# t0 = ttopfun.subs(y, 0)
# b0 = btopfun.subs(y, 0)
# eq = smp.Eq(-nd_des.subs(y, 0), sigma(kc, t0, b0))
# kc1 = smp.solve(eq, kc)[0]
# print(t0)
# print(b0)
# print(kc1)
# a_over_b = smp.solve(kc1-kfun, ab)
# print(a_over_b)
# a = a_over_b[-1]*b0
# print(a)
# tavg = 0.5*(ttopfun.subs(y, 0)+ttopfun.subs(y, 0+a))
# bavg = 0.5*(btopfun.subs(y, 0)+btopfun.subs(y, 0+a))
# eq1 = smp.Eq(-nd_des.subs(y, 0), sigma(kc, tavg, bavg))
# print(-nd_des.subs(y, 0), sigma(kc, tavg, bavg).evalf())
# kc2 = smp.solve(eq1, kc)[0]
# print(kc2)
# a_over_b1 = smp.solve(kc2-kfun, ab)
# print(a_over_b1)
# a1 = a_over_b1[-1]*bavg
# print(a1)

margin_list = [1.5] * 1000
plt.figure()
plt.xlabel('Span-wise location [m]')
plt.ylabel('Margin of safety')
plt.plot(span, margin_list)
plt.show()

# # top plate
# tavg = (ttopfun.subs(y, 0)+ttopfun.subs(y, 0.5*var.Span))/2
# bavg = (btopfun.subs(y, 0)+btopfun.subs(y, 0.5*var.Span))/2
# sigma_max = 1.5*2.45e8
# sigma_top = sigma(kc, tavg, bavg)
# eq1 = smp.Eq(sigma_max, sigma_top)
# kc = smp.solve(eq1, kc)[0]
# a_b1 = smp.solve(kc-kfun, ab)[-1]
# print(a_b1)
# a = a_b1*bavg
# print(a)
# print(smp.ceiling((0.5*var.Span)/a))

# print('Starting the calculation, stay patient.')

# # Top plate
# dy1 = 0
# alist1 = []
# dylist1 = []
# alist2 = []
# dylist2 = []
# while dy1 <= 0.5*var.Span:
#   if dy1 <= 0.35*0.5*var.Span:
#     eq1 = smp.Eq(nd_des.subs(y, dy1), sigmatopfun.subs(y, dy1))
#     kc1 = abs(smp.solve(eq1, kc)[0])
#     if len(smp.solve(kc1-kfun, ab))==0:
#       dylist1.append(dy1)
#       if len(alist1) ==0:
#         alist1.append(1)
#         dy1 += 1
#         continue
#       if len(alist1) !=0:
#         index = alist1.index(max(alist1))
#         alist1.append(alist1[index])
#         dy1 += max(alist1)
#         continue
#       continue
#     a_b1 = smp.solve(kc1-kfun, ab)[-1]
#     a1 = (a_b1*btopfun.subs(y, dy1)).evalf()
#     bavg1 = (btopfun.subs(y, dy1)+btopfun.subs(y, (dy1+a1).doit()))/2
#     tavg1 =(ttopfun.subs(y, dy1)+ttopfun.subs(y, (dy1+a1).doit()))/2
#     eq2 = smp.Eq(nd_des.subs(y, dy1), sigma(kc, tavg1, bavg1))
#     kc2 = abs(smp.solve(eq2, kc)[0])
#     if len(smp.solve(kc1-kfun, ab))==0:
#       dylist1.append(dy1)
#       index = alist1.index(max(alist1))
#       alist1.append(alist1[index])
#       dy1 += max(alist1)
#       continue
#     a_b2 = smp.solve(kc2-kfun, ab)[-1]
#     a2 = a_b2*bavg1
#     alist1.append(a2)
#     dylist1.append(dy1)
#     dy1 += a2
#   else: 
#     eq1 = smp.Eq(nd_des.subs(y, dy1), sigmatopfun.subs(y, dy1))
#     kc1 = abs(smp.solve(eq1, kc)[0])
#     if len(smp.solve(kc1-kfun, ab))==0:
#       dylist2.append(dy1)
#       if len(alist2) == 0:
#         alist2.append(1)
#         dy1 += 1
#         continue
#       if len(alist2) != 0:
#         index = alist2.index(max(alist2))
#         alist2.append(alist2[index])
#         dy1 += max(alist2)
#         continue
#       continue
#     a_b1 = smp.solve(kc1-kfun, ab)[-1]
#     a1 = (a_b1*btopfun.subs(y, dy1)).evalf()
#     bavg1 = (btopfun.subs(y, dy1)+btopfun.subs(y, (dy1+a1).doit()))/2
#     tavg1 =(ttopfun.subs(y, dy1)+ttopfun.subs(y, (dy1+a1).doit()))/2
#     eq2 = smp.Eq(nd_des.subs(y, dy1), sigma(kc, tavg1, bavg1))
#     kc2 = abs(smp.solve(eq2, kc)[0])
#     if len(smp.solve(kc2-kfun, ab))==0:
#       dylist2.append(dy1)
#       index = alist2.index(max(alist2))
#       alist2.append(alist2[index])
#       dy1 += max(alist2)
#       continue
#     a_b2 = smp.solve(kc2-kfun, ab)[-1]
#     a2 = a_b2*bavg1
#     alist2.append(a2)
#     dylist2.append(dy1)
#     dy1 += a2
# ribs_no_top = len(dylist1)+len(dylist2)-1
# ribs_location_top_list = dylist1[1:]+dylist2
# print('Top plate calculation is done, moving on to the bottom plate.')

# # Bottom plate
# dy2 = 0
# alist3 = []
# dylist3 = []
# alist4 = []
# dylist4 = []
# while dy2 <= 0.5*var.Span:
#   if dy2 <= 0.35*0.5*var.Span:
#     eq3 = smp.Eq(nd_des.subs(y, dy2), sigmabotfun.subs(y, dy2))
#     kc3 = abs(smp.solve(eq3, kc)[0])
#     if len(smp.solve(kc3-kfun, ab))==0:
#       dylist3.append(dy2)
#       if len(alist3) ==0:
#         alist3.append(1)
#         dy2 += 1
#         continue
#       if len(alist3) !=0:
#         index = alist3.index(max(alist3))
#         alist3.append(alist3[index])
#         dy2 += max(alist3)
#         continue
#       continue
#     a_b3 = smp.solve(kc3-kfun, ab)[-1]
#     a3 = (a_b3*bbotfun.subs(y, dy2)).evalf()
#     bavg2 = (bbotfun.subs(y, dy2)+bbotfun.subs(y, (dy2+a3).doit()))/2
#     tavg2 =(tbotfun.subs(y, dy2)+tbotfun.subs(y, (dy2+a3).doit()))/2
#     eq4 = smp.Eq(nd_des.subs(y, dy2), sigma(kc, tavg2, bavg2))
#     kc4 = abs(smp.solve(eq4, kc)[0])
#     if len(smp.solve(kc3-kfun, ab))==0:
#       dylist3.append(dy2)
#       if len(alist3) == 0:
#         alist3.append(1)
#         dy2 += 1
#         continue
#       if len(alist3) != 0:
#         index = alist3.index(max(alist3))
#         alist3.append(alist3[index])
#         dy2 += max(alist3)
#         continue
#       continue
#     a_b4 = smp.solve(kc4-kfun, ab)[-1]
#     a4 = a_b4*bavg2
#     alist3.append(a4)
#     dylist3.append(dy2)
#     dy2 += a4
#   else: 
#     eq3 = smp.Eq(nd_des.subs(y, dy2), sigmabotfun.subs(y, dy2))
#     kc3 = abs(smp.solve(eq3, kc)[0])
#     if len(smp.solve(kc3-kfun, ab))==0:
#       dylist4.append(dy2)
#       if len(alist4) ==0:
#         alist4.append(1)
#         dy2 += 1
#         continue
#       if len(alist4) != 0:
#         index = alist4.index(max(alist4))
#         alist4.append(alist4[index])
#         dy2 += max(alist4)
#         continue
#       continue
#     a_b3 = smp.solve(kc3-kfun, ab)[-1]
#     a3 = (a_b3*bbotfun.subs(y, dy2)).evalf()
#     bavg2 = (bbotfun.subs(y, dy2)+bbotfun.subs(y, (dy2+a3).doit()))/2
#     tavg2 =(tbotfun.subs(y, dy2)+tbotfun.subs(y, (dy2+a3).doit()))/2
#     eq4 = smp.Eq(nd_des.subs(y, dy2), sigma(kc, tavg2, bavg2))
#     kc4 = abs(smp.solve(eq4, kc)[0])
#     if len(smp.solve(kc4-kfun, ab))==0:
#       dylist4.append(dy2)
#       if len(alist4) ==0:
#         alist4.append(1)
#         dy2 += 1
#         continue
#       if len(alist4) !=0:
#         index = alist4.index(max(alist4))
#         alist4.append(alist4[index])
#         dy2 += max(alist4)
#         continue
#       continue
#     a_b4 = smp.solve(kc4-kfun, ab)[-1]
#     a4 = a_b4*bavg2
#     alist4.append(a4)
#     dylist4.append(dy2)
#     dy2 += a4
# ribs_no_bot = len(dylist3)+len(dylist4)-1
# ribs_location_bot_list = dylist3[1:]+dylist4
# print('Bottom plate calculation is done, moving on to the comparison.')

# # Comparing the panel lengths of top and bottom plates to determine the driving one
# li1 = alist1+alist2
# li2 = alist3+alist4
# comp_list = []
# for i in li1:
#   index = li1.index(i)
#   if i > li2[index]:
#     comp_list.append(1)
#   else:
#     comp_list.append(0)
# diff = 0
# if len(ribs_location_bot_list)>len(ribs_location_top_list):
#   diff += len(ribs_location_bot_list)-len(ribs_location_top_list)
# else:
#   diff += 0
# for i in range(0, diff):
#   comp_list.append(1)

# # Forming the a_list according to comparison
# a_list = []
# for i in range(0, len(comp_list)):
#   if comp_list[i] == 1:
#     a_list.append(li2[i])
#   else:
#     a_list.append(li1[i])

# max_repeating = 0
# if max(alist2)>max(alist4):
#   max_repeating+=max(alist2)
# else:
#   max_repeating+=max(alist4)
# repeating_index = 0
# for i in range(0, len(a_list)-1):
#   if a_list[i]==a_list[i+1]:
#     if repeating_index == 0:
#       repeating_index += i
# first_or_already = 0
# if li2[repeating_index-1]==li2[repeating_index]:
#   first_or_already += 1
# max_repeating_list = []
# for i in range(0, len(a_list)-repeating_index):
#   max_repeating_list.append(max_repeating)
# if first_or_already == 1:
#   a_list[repeating_index:] = max_repeating_list
# else:
#   a_list[repeating_index+1:] = max_repeating_list[1:]
# print('Comparison is done, determining the rib locations.')

# # Determining the rib locations according to the comparison
# a1 = 0
# ribs_location_list = []
# for i in a_list[0:]:
#   a1 += i
#   ribs_location_list.append(a1)
# exceed_index = 0
# for i in ribs_location_list:
#   if i> (0.5*var.Span):
#     exceed_index = ribs_location_list.index(i)
#     break

# # Printing the results
# if exceed_index != 0:
#   print('Number of ribs: {}'.format(len(ribs_location_list)-len(ribs_location_list[exceed_index:])))
#   print('Rib locations:')
#   print(ribs_location_list[0:-1*len(ribs_location_list[exceed_index:])])
# else:
#   print('Number of ribs: {}'.format(len(ribs_location_list[0:])))
#   print('Rib locations:')
#   print(ribs_location_list[0:])

# print('Rib spacings:')
# print(a_list[1:-1*len(ribs_location_list[exceed_index:])])

# plt.figure()
# plt.xlabel('Number of ribs')
# plt.ylabel('Ribs spacing [m]')
# plt.plot([*range(1, 1+len(a_list[1:-1*len(ribs_location_list[exceed_index:])]))], a_list[1:-1*len(ribs_location_list[exceed_index:])], marker='o')
# plt.figure()
# plt.xlabel('Span-wise location of each rib [m]')
# plt.scatter(ribs_location_list, [0]*len(ribs_location_list[0:]), label='o')
# plt.show()