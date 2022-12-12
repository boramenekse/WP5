#Solver figuring out
import pandas as pd
from pulp import * 
import scipy as sc
from scipy.optimize import leastsq, curve_fit
import numpy as np
import sympy as smp

#
#df = pd.DataFrame({'Variable': ['x1', 'x2'],
#                   'Price': [45, 20],
#                   'Cost': [30, 10],
#                   'Demand': [2000, 8000]
#                    })
#
#df.eval('Margin=Price-Cost', inplace=True)
#
## set the dictionary for each feature
#prob = LpProblem('Sell', LpMaximize) # Objective function
#
#inv_item = list(df['Variable']) # Variable name
#
#margin = dict(zip(inv_items, df['Margin'])) # Function 
#
#demand = dict(zip(inv_items, df['Demand'])) # Function
#
## next, we are defining our decision variables as investments and are adding a few parameters to it
#inv_vars = LpVariable.dicts('Variable', inv_items, lowBound=0, cat='Integar')
#
## set the decision variables
## all add in the problem setting
#prob += lpSum([inv_vars[i] * margin[i] for i in inv_items])
#
## Constraint
#prob += lpSum([inv_vars[i] for i in inv_items]) <= 8000, 'Total Demand'
#
## Constraint
#prob += inv_vars['x1'] <= 2000, 'x1 Demand'
#
## Constraint
#prob += inv_vars['x2'] <= 8000, 'x2 Demand'
#



##lele = np.linspace(0, 9, 10)
#mierda = 5
#def puta(g):
#    g**2
##print(lele)
#STAN = sc.optimize.leastsq(lele, puta(mierda))


# y = np.array([12, 8, 11, 7, 5, 2, 3, 5, 6, 4, 5, 7, 8, 13, 19, 22, 25])
# x = np.array(range(len(y)))

# def func1(params, x, y):
#     a, b, c = params[0], params[1], params[2]
#     residual = y - (a*x**2+b*x+c)
#     return residual

# params=[0, 0, 0]

# result = leastsq(func1, params, (x, y))  # first one is residual, which is minimised, second one is parameters to change, third one are
# a, b, c = result[0][0], result[0][1], result[0][2]
# yfit1 = a*x**2+b*x+c

y = np.arange(0.0, (6), 1.0)
z = [3,5,7,6,5,4]
def test_function(x, A, B, C, D):
    y = A*(x**3) + B*(x**2) + C*x + D
    return y 

Parameters, covariance = curve_fit(test_function, y, z)
#print(Parameters)
# print(covariance)

def trial(a, b, c, d, e, f, g, h, i, j, x):
  return (a* x**5+b* x**4+c* x**3+d* x**2+e* x+f)/(g* x**3+h* x**2+i* x+j)

x = smp.symbols('x', real = True)

result = smp.integrate(trial(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, x), x)
# print(result)
result2 = smp.integrate(result, x)
print(smp.lambdify([x], result2))
print(smp.lambdify([x], result2)(1))



