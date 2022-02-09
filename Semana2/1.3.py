import numpy as np
import scipy as sp
from scipy import integrate

def function_Riemann(t):
    f = 1/np.pi * (t/12*(t**2-np.pi**2))**2
    return f

def f_Riemann():
    est,error = integrate.quad(lambda t: function_Riemann(t),-np.pi,np.pi)
    return "Para dseta igual a 6 el valor es " + str(est) + " y su error es +/- " + str(error)

print(f_Riemann())