import sympy as sym

# Creacion de symbols para uso de sym.integrate

t = sym.Symbol('t',Real=True)
h = sym.Symbol('h',Real=True)

# Solucion para los 3 puntos Adams-Bashforth

F1 = (t-(-h))*(t-(-2*h))/(2*h**2)
F2 = -(t-(0))*(t-(-2*h))/(h**2)
F3 = (t-(-h))*(t-(0))/(2*h**2)

C1 = sym.integrate(F1, (t,0,h))
C2 = sym.integrate(F2, (t,0,h))
C3 = sym.integrate(F3, (t,0,h))

# donde C1,C2,C3 son los coeficientes del metodo de Adams-Bashforth

# Solucion para los 4 coeficientes Adams-Bashforth

F1 = (t-(-h))*(t-(-2*h))*(t-(-3*h))/(6*h**3)
F2 = -(t-(0))*(t-(-2*h))*(t-(-3*h))/(2*h**3)
F3 = (t-(-h))*(t-(0))*(t-(-3*h))/(2*h**3)
F4 = -(t-(-h))*(t-(0))*(t-(-2*h))/(6*h**3)

C1= sym.integrate(F1,(t,0,h))
C2= sym.integrate(F2,(t,0,h))
C3= sym.integrate(F3,(t,0,h))
C4= sym.integrate(F4,(t,0,h))

# donde C1,C2,C3 son los coeficientes del metodo de Adams-Bashforth

# Solucion para los 3 puntos de Adams-Moulton

F1 = (t-(0))*(t-(-h))/(2*h**2)
F2 = -(t-(h))*(t-(-h))/(h**2)
F3 = (t-(0))*(t-(h))/(2*h**2)

C1 = sym.integrate(F1, (t,0,h))
C2 = sym.integrate(F2, (t,0,h))
C3 = sym.integrate(F3, (t,0,h))

# donde C1,C2,C3 son los coeficientes del metodo de Adams-Moulton

# Solucion para los 4 puntos de Adams-Moulton

F1 = (t-(0))*(t-(-h))*(t-(-2*h))/(6*h**3)
F2 = -(t-(h))*(t-(-h))*(t-(-2*h))/(2*h**3)
F3 = (t-(0))*(t-(h))*(t-(-2*h))/(2*h**3)
F4 = -(t-(0))*(t-(h))*(t-(-h))/(6*h**3)

C1= sym.integrate(F1,(t,0,h))
C2= sym.integrate(F2,(t,0,h))
C3= sym.integrate(F3,(t,0,h))
C4= sym.integrate(F4,(t,0,h))

# donde C1,C2,C3,C4 son los coeficientes del metodo de Adams-Moulton