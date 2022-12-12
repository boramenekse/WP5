import sympy as smp # that's all you need :)
# Sympy stands for symbolic python (numpy -> numerical python, scipy -> scientific python)
# This library allows you to do all the math in symbolic way in python as if you were doing it on a paper -> Helpful for detecting errors since there is no need to plug numbers in
# It would be better to start by showing an example
# First, always define your symbols. These are basically the independent variables we use while doing math on paper. It's an useful habit to give as much arguments as possible to your symbols
x = smp.symbols('x')
y = smp.symbols('y')
z = smp.symbols('z')
a = smp.symbols('a', real=True) # First string determines how the symbol will look like
h = smp.symbols('h', real=True)
d = smp.symbols('d', real=True, positive=True) 
r = smp.symbols('r', real=True, positive=True) # Radius cannot be negative so, I specify it to make sympy's work easier. Those arguments help sympy to do the math (makes it easier or faster)
alpha = smp.symbols('\u03B1', real=True) # use unicode for greek letters
beta = smp.symbols('\u03B2', real=True)
rho = smp.symbols('\u03C1', real=True)
sigma = smp.symbols('\u03C3', real=True)
theta = smp.symbols('\u03B8', real=True)
eta = smp.symbols('\u03B7', real=True)
xi = smp.symbols('\u03BE', real=True)
mx = smp.symbols('Mx', real=True)
my = smp.symbols('My', real=True)
ixx = smp.symbols('Ixx', real=True)
iyy = smp.symbols('Ixy', real=True)
ixy = smp.symbols('Iyy', real=True)
t = smp.symbols('t', real=True)
t1 = smp.symbols('t1', real=True)
t2 = smp.symbols('t2', real=True)
t3 = smp.symbols('t3', real=True)
vx = smp.symbols('Vx', real=True)
vy = smp.symbols('Vy', real=True)
s = smp.symbols('s', real=True)
s1 = smp.symbols('s1', real=True)
s2 = smp.symbols('s2', real=True)
s3 = smp.symbols('s3', real=True)
s4 = smp.symbols('s4', real=True)
s5 = smp.symbols('s5', real=True)
s6 = smp.symbols('s6', real=True)

def qs(vx, vy, ixx, iyy, ixy, t, x, y, s):
  # sympy does't take integration constants into account so, while doing an indefinite integral always think about your integration constants
  # smp.integrate(function, (variable to integrate for, lower bound, upper bound))
  # note that you can put variables to bounds just as I did here
  # for an indefinite integral -> smp.integrate(function, variable to integrate for)
  exp = -((vy*iyy-vx*ixy)/(ixx*iyy-ixy**2))*smp.integrate(t*y, (s, 0, s)) - ((vx*ixx-vy*ixy)/(ixx*iyy-ixy**2))*smp.integrate(t*x, (s, 0, s))
  return exp

# It is so useful to work with rationals instead of floating numbers in sympy, especially while doing a polynomial division. Otherwise, for so complex functions, sympy can crash
# You can understand it's crashed if the programs works more than 30 seconds 
h = smp.Rational(5, 100) # 5/100
t1 = smp.Rational(2, 1000)
t2 = smp.Rational(3, 1000)
t3 = smp.Rational(4, 1000)
ybar = (h*2*h*t2+2*h*2*h*t1)/(h*t3+2*h*t2+2*h*t1)
xbar = (h*2*h*t1+h*2*h*t2+smp.Rational(3,2)*h*h*t3)/(h*t3+2*h*t2+2*h*t1) # Use smp.Rational even for simple cases
ixx = smp.Rational(8,12)*t2*h**3+2*h*t2*(h-ybar)**2+2*h*t1*(2*h-ybar)**2+h*t3*ybar**2
iyy = smp.Rational(8,12)*t1*h**3+smp.Rational(1,12)*t3*h**3+2*h*t1*(h-xbar)**2+2*h*t2*(h-xbar)**2+h*t3*(smp.Rational(3,2)*h-xbar)**2
ixy = 2*h*t1*(2*h-ybar)*(h-xbar)+2*h*t2*(h-ybar)*(h-xbar)+h*t3*(-ybar)*(smp.Rational(3,2)*h-xbar)
qs54 = qs(vx, vy, ixx.simplify(), iyy.simplify(), ixy.simplify(), t3, 2*h-xbar-s4, -ybar, s4).simplify() # It should be clear what .simplify() does -> sympy_expression.simplify() use like this
m = -smp.integrate(qs54*2*h, (s4, 0, h)) # Output: 17*Vx/520 - 7*Vy/520 
print(m.expand()) # sympy_expression.expand expands the expression :o
# Let's try to solve the following equation: vy*ξ + vx*η + m for ξ and η 
print(m.args) # sympy_expression.args breaks down an expression into seperate parts -> Returns a tuple
# m is -> 17*Vx/520 - 7*Vy/520 So, m.args would be (-7*Vy/520, 17*Vx/520)  
print(m.args[0].args) # output: (-7/520, Vy) 
print(m.args[1].args) # output: (17/520, Vx)
print('Eta is: {}'.format(m.args[0].args[0]))
print('Xi is: {}'.format(m.args[1].args[0]))
print('Xi is: {1} and Eta is:{0}'.format(m.args[0].args[0], m.args[1].args[0])) #more than 1 bracked inside the string
# Or I could just use sympy_expression.coeff(the variable to get the coefficient for)
print(m.coeff(vy)) # Output: -7/520
print(m.coeff(vx)) # Output: 17/520
# I can get them in floating numbers by using sympy_expression.evalf() -> Being used for evaluating a numerical expression into a floating number
print(m.coeff(vy).evalf()) # Output: -0.0134615384615385
print(m.coeff(vx).evalf()) # Output: 0.0326923076923077
# If I didn't enter values for h, t1, t2, and t3, the simplified output will be this:
# 2*h*t3*(8*Vx*t1*t2 + 3*Vx*t1*t3 + 2*Vx*t2**2 + Vx*t2*t3 - 6*Vy*t1**2 - 3*Vy*t1*t2)/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2)
# From that the coefficients of Vx and Vy can be obtained as:
# print(m.simplify().expand().coeff(vx).simplify())
# print(m.simplify().expand().coeff(vy).simplify())
# First, I simplify it, and then, expand the simplified version, which makes expression less complex than just saying m.expand()
# m.simplify().expand() gives: 16*Vx*h*t1*t2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 6*Vx*h*t1*t3**2/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 4*Vx*h*t2**2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 2*Vx*h*t2*t3**2/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) - 12*Vy*h*t1**2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) - 6*Vy*h*t1*t2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) 
# Then, with .coeff(vx), I get the summation of all coefficients of vx, and following is the output:
# 16*h*t1*t2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 6*h*t1*t3**2/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 4*h*t2**2*t3/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2) + 2*h*t2*t3**2/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2)
# Since it looks ugly, I use .simplify() and get the following:
# 2*h*t3*(8*t1*t2 + 3*t1*t3 + 2*t2**2 + t2*t3)/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2)
# print(m.simplify().expand().coeff(vy).simplify()) gives: 6*h*t1*t3*(-2*t1 - t2)/(16*t1**2*t2 + 24*t1**2*t3 + 4*t1*t2**2 + 16*t1*t2*t3 + 3*t1*t3**2 + 2*t2**2*t3 + t2*t3**2)
# By using .subs(), results can be proven such as:
# print(m.simplify().expand().coeff(vx).simplify().subs([(h, smp.Rational(5, 100)), (t1, smp.Rational(2, 1000)), (t2, smp.Rational(3, 1000)), (t3, smp.Rational(4, 1000))])) -> 17/520
# print(m.simplify().expand().coeff(vy).simplify().subs([(h, smp.Rational(5, 100)), (t1, smp.Rational(2, 1000)), (t2, smp.Rational(3, 1000)), (t3, smp.Rational(4, 1000))])) -> -7/520
# THIS IS THE POWER OF SYMPY!!!

# Expression manipulation
# Let's start with the most basic and commonly used one -> substitution
expr1 = smp.cos(x) + 1
expr2 = expr1.subs(x, y) # First the variable you wanna change, then the variable(or this can be a number as well) you wanna substitute to the first variable
# Important thing to note is that .subs is an immutable method, which means subs does not modify the expression in-place -> Therefore, I defined a new variable to store the result
print(expr2) # output: cos(y) + 1

expr3 = x**3 + 4*x*y - z
expr4 = expr3.subs([(x, 2), (y, 4), (z, 0)]) # For multiple substitutions -> Output: 40

# Rewriting
print(smp.tan(theta).rewrite(smp.cos(theta))) # Output: cos(θ - pi/2)/cos(θ) -> .rewrite(a function)

# Creating an equation
print(smp.Eq(sigma, (mx*y)/ixx)) # .Eq(lhs, rhs) -> sigma = (mx*y)/ixx -> however, when you print it, you'll see it is not shown like that -> I recommend not using this actually -> it might be easier to work with expressions such as: expr_for_sigma = (mx*y)/ixx 

# Converting strings to sympy expressions
str_expr = "x**2 + 3*x - 1/2"
expr = smp.sympify(str_expr)

# Converting sympy expressions to python function with the same functionality of                                   def any_function(parameters): return output
m_fun = smp.lambdify([vx, vy], m) # smp.lambdify([list of parameters], expression to convert)
# Now, you can use m_fun to plot the function or calculate it at any value -> you can do the second thing in sympy as well, your choice -> You cannot print m_fun, but you can print the sympy expression of it, which is m

# Polynomial division
f = s6*x**5 + s5*x**4 + s4*x**3 + s3*x**2 + s2*x+ s1
g = t3*x**3 + t2*x**2 + t1*x + t # t + x**3/250 + 3*x**2/1000 + x/500
print(smp.div(f, g, domain='QQ')) # Returns a tuple in the form: (quotient, remainder) Output: (250*s4 + 250*s5*x - 375*s5/2 + 250*s6*x**2 - 375*s6*x/2 + 125*s6/8, s1 + s2*x + s3*x**2 - 250*s4*t - 3*s4*x**2/4 - s4*x/2 - 250*s5*t*x + 375*s5*t/2 + s5*x**2/16 + 3*s5*x/8 - 250*s6*t*x**2 + 375*s6*t*x/2 - 125*s6*t/8 + 21*s6*x**2/64 - s6*x/32) 

# Solving equations
eq1  = smp.Eq(x + y - 2*z, 0) # x+y-2*z = 0
print(smp.solve(eq1, x)) # [-y+2*z] -> solved eq1 for x -> if dict=True is not included, output is given as a list
print(smp.solve(eq1, x, dict=True)) # -> [{x: -y + 2*z}]
eq2  = smp.Eq(y+4*z, 0) # y+4*z = 0
print(smp.solve([eq1, eq2], [x, y], dict=True)) # [{x: 6*z, y: -4*z}]
print(smp.solve([x + y - 2*z, y + 4*z], [x, y], dict=True)) # [{x: 6*z, y: -4*z}]
print(smp.solve([x + y - 2*z, y + 4*z], [x, y])) # {x: 6*z, y: -4*z}
# These two are the same -> an example of solving systems of equations
# Another way of solving systems of equations is using matrices
eq3 = smp.Eq(x+y+z, 1)
eq4 = smp.Eq(x+y+2*z, 3)
matrix = smp.Matrix(([1, 1, 1, 1], [1, 1, 2, 3])) # Augmented matrix form
print(smp.linsolve(matrix, (x, y, z))) # Output: {(-y - 1, y, 2)} -> x = -y-1, y=y, z = 2 -> this means y is a free variable, and equations are satisfied as long as x=-y-1 and z=2 
# In order to get the solutions:
results = smp.linsolve(matrix, (x, y, z)) # returns a sympy set
print(list(results)) # Convert results sympy set to a list Output: [(-y - 1, y, 2)] 
print(list(results)[0]) # getting the first thing in the list, and since it is a tuple, get the first thing in that tuple

# Step functions
# Because of a convention in electronics and signal processing, when you put the value that makes the function equal to zero, heaviside returns 1/2. However, you can change it by adding an argument as I did. 
step1 = smp.Heaviside(x-2) # at 2 -> gives 1/2
step2 = smp.Heaviside(x-2, 0) # at 2 -> gives 0
print(step2.subs(x, 3)) # Output: 1
print(step2.subs(x, 1)) # Output: 0
print(step2.subs(x, 2)) # Output: 0
step3 = smp.Heaviside(x-2, 1)
print(step3.subs(x, 2)) # Output: 1
fun1 = (x**2+1)*step2
print(fun1)