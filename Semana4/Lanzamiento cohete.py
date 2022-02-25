import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import integrate

G = 6.6710e-11 
mT = 5.9736e24 
rT = 6.3781e6 
mL = 0.07349e24 
rL = 1.7374e6 
d = 3.844e8 
w = 2.6617e-6

v0 = 11190/d

phi0 = 0/180*np.pi
theta0 = 25/180*np.pi
r0 = rT/d

delta = (G*mT)/(d**3)
miu = mL/mT