#yeah
import sympy as smp

#Moment Distribution
#loading case 16 and loading case 12
#jump in the moment diagram 11.69m

x = smp.symbols('x', real=True)
M1_12 = [-0.19569866609295541, 6.6801804985865445, 347.650403262873, -3851.0392032213094, -529912.4293065271, 13329458.866177564] 
M2_12 = [-0.19569866608662626, 6.680180498386035, 347.65040326524246, -3851.039203234463, -684065.9293064886, 14011505.048648307]
#Moment distribution constants for LC-16: [0.40016072099856226, -9.397942238957981, -1047.384217723905, 6708.870613191617, 1582191.8673072485, -33488163.089606084] Part 2:[0.40016072098555516, -9.397942238637903, -1047.3842177269048, 6708.870613212022, 1967575.6173070923, -39113307.35963459] 



#Normal force diagram (until 11.69)
#force is -500000 N (compression)
md_fun1_12 = smp.nsimplify(round(M1_12[0], 6))*x**5+smp.nsimplify(round(M1_12[1], 6))*x**4+smp.nsimplify(round(M1_12[2], 6))*x**3+smp.nsimplify(round(M1_12[3], 6))*x**2+smp.nsimplify(round(M1_12[4], 6))*x+smp.nsimplify(round(M1_12[5], 6))
md_fun2_12 = smp.nsimplify(round(M2_12[0], 6))*x**5+smp.nsimplify(round(M2_12[1], 6))*x**4+smp.nsimplify(round(M2_12[2], 6))*x**3+smp.nsimplify(round(M2_12[3], 6))*x**2+smp.nsimplify(round(M2_12[4], 6))*x+smp.nsimplify(round(M2_12[5], 6))

heaviside = smp.Heaviside(x-11.69, 0)
m_fun1 = md_fun1_12-md_fun1_12*heaviside+md_fun2_12*heaviside
print(m_fun1.simplify())