import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

# discretizacion 

n = 51
x = np.linspace(0,1,n)
h = x[1]-x[0]
v = 1.

def initU(n,v,h):
    
    t = np.zeros((n,n))
    
    for i in range(n):
        for j in range(n):
            t[i,j] = j*v
    for i in range(n):
        t[i,-1] = t[i,-2] + v*h
    for j in range(n):
        t[1][j] = t[0][j]
    for i in range(n):
        if i <= 5 or i >= 30:
            t[i,0] = 0.
    for j in range(1,n-1):
        t[-1][j] = t[-2][j]

    return t

def initW(n):

    t = np.zeros((n,n))

    t[0,:] = 0.
    t[-1:] = 0.
    t[:0] = 0.
    t[:,-1] = 0.
    
    return t

# creacion variables
u =  initU(n,v,h)
w = initW(n)

def choque(u,w,h):

    for i in range(5,30):
        u[i,5] = 0
        w[i,5] = -2*(u[i-1,5] - u[i,5])/(h**2)
        u[i,25] = 0 
        w[i,25] = -2*(u[i+1,25] - u[i,25])/(h**2)
    for j in range(5,25):
        u[5,j] = 0 
        w[5,j] = -2*(u[5,j-1] - u[5,j])/(h**2)
        u[30,j] = 0 
        w[30,j] = -2*(u[30,j+1] - u[30,j])/(h**2)

    return u,w

u,w = choque(u,w,h)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
X,Y = np.meshgrid(x,x)