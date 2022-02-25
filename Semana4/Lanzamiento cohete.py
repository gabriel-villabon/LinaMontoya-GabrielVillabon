import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from scipy import integrate

G = 6.67e-11 
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

def rpri (r,phi,t):
    return (1+r**2-2*r*np.cos(phi-w*t))**(1/2)
    
def pr0_ (v,theta,phi):
    return v*np.cos(theta - phi)

def pphi0_ (r,v,theta,phi):
    return r*v*np.cos(theta - phi)

pr0 = pr0_(v0,theta0,phi0)
pphi0 = pphi0_ (r0,v0,theta0,phi0)

def f1(r,phi,pr,pphi,t):
    return pr

def f2(r,phi,pr,pphi,t):
    return pphi/(r**2)

def f3(r,phi,pr,pphi,t):
    return (pphi**2)/(r**3) - delta*(1/(r**2)+miu/(rpri(r,phi,t)**3)(r-np.cos(phi-w*t)))

def f4(r,phi,pr,pphi,t):
    return - ((delta*miu*r)/(rpri(r,phi,t)))*np.sin(phi-w*t)

def GetRK4(r0,phi0,pr0,pphi0,t0,tf,h):
    
    N = int((tf-t0)/h) +1
    
    t = np.linspace(t0,tf,N)
    
    r = np.zeros([N])
    phi = np.zeros([N])
    pr = np.zeros([N])
    pphi = np.zeros([N])
    
    k1 = np.zeros(N)
    k2 = np.zeros(N)
    k3 = np.zeros(N)
    k4 = np.zeros(N)
    
    r[0] = r0
    phi[0] = phi0
    pr[0] = pr0
    pphi[0] = pphi0
    
    for i in range(1,N):
        
        k1[0] = h*f1(r[i-1],phi[i-1],pr[i-1],pphi[i-1],t[i-1])
        k1[1] = h*f2(r[i-1],phi[i-1],pr[i-1],pphi[i-1],t[i-1])
        k1[2] = h*f3(r[i-1],phi[i-1],pr[i-1],pphi[i-1],t[i-1])
        k1[3] = h*f4(r[i-1],phi[i-1],pr[i-1],pphi[i-1],t[i-1])
        
        k2[0] = f1(r[i-1]+0.5*k1[0],phi[i-1]+0.5*k1[1],pr[i-1]+0.5*k1[2],pphi[i-1]+0.5*k1[3],t[i-1]+0.5*h)
        k2[1] = f2(r[i-1]+0.5*k1[0],phi[i-1]+0.5*k1[1],pr[i-1]+0.5*k1[2],pphi[i-1]+0.5*k1[3],t[i-1]+0.5*h)
        k2[2] = f3(r[i-1]+0.5*k1[0],phi[i-1]+0.5*k1[1],pr[i-1]+0.5*k1[2],pphi[i-1]+0.5*k1[3],t[i-1]+0.5*h)
        k2[3] = f4(r[i-1]+0.5*k1[0],phi[i-1]+0.5*k1[1],pr[i-1]+0.5*k1[2],pphi[i-1]+0.5*k1[3],t[i-1]+0.5*h)
        
        k3[0] = f1(r[i-1]+0.5*k2[0],phi[i-1]+0.5*k2[1],pr[i-1]+0.5*k2[2],pphi[i-1]+0.5*k2[3],t[i-1]+0.5*h)
        k3[1] = f2(r[i-1]+0.5*k2[0],phi[i-1]+0.5*k2[1],pr[i-1]+0.5*k2[2],pphi[i-1]+0.5*k2[3],t[i-1]+0.5*h)
        k3[2] = f3(r[i-1]+0.5*k2[0],phi[i-1]+0.5*k2[1],pr[i-1]+0.5*k2[2],pphi[i-1]+0.5*k2[3],t[i-1]+0.5*h)
        k3[3] = f4(r[i-1]+0.5*k2[0],phi[i-1]+0.5*k2[1],pr[i-1]+0.5*k2[2],pphi[i-1]+0.5*k2[3],t[i-1]+0.5*h)
        
        k4[0] = h*f1(r[i-1]+k3[0],phi[i-1]+k3[1],pr[i-1]+k3[2],pphi[i-1]+k3[3],t[i-1]+h)
        k4[1] = h*f2(r[i-1]+k3[0],phi[i-1]+k3[1],pr[i-1]+k3[2],pphi[i-1]+k3[3],t[i-1]+h)
        k4[2] = h*f3(r[i-1]+k3[0],phi[i-1]+k3[1],pr[i-1]+k3[2],pphi[i-1]+k3[3],t[i-1]+h)
        k4[3] = h*f4(r[i-1]+k3[0],phi[i-1]+k3[1],pr[i-1]+k3[2],pphi[i-1]+k3[3],t[i-1]+h)
        
        r[i] = r[i-1] + (k1[0]+2*k2[0]+2*k3[0]+k4[0])/6
        phi[i] = phi[i-1] + (k1[1]+2*k2[1]+2*k3[1]+k4[1])/6
        pr[i] = pr[i-1] + (k1[2]+2*k2[2]+2*k3[2]+k4[2])/6
        pphi[i] = pphi[i-1] + (k1[3]+2*k2[3]+2*k3[3]+k4[3])/6
            
    return t,r,phi

t,r,phi = GetRK4(r0,phi0,pr0,pphi0,0,7*10**5,10)

t_red = []
r_red = []
phi_red = []

for i in range(len(t)):
    if i%1000 == 0:
        t_red.append(t)
        r_red.append(r)
        phi_red.append(phi)
        
t_red = np.array(t_red)
r_red = np.array(r_red)
phi_red = np.array(phi_red)

x,y = r_red*np.cos(phi_red), r_red*np.sin(phi_red)
xl,yl = r_red*np.cos(w*t_red), r_red*np.sin(w*t_red)

fig = plt.figure()
plt.scatter(x,y)
plt.scatter(xl,yl)
plt.xlim(-1,1)
plt.ylim(-1,1)

fig1 = plt.figure()
ax = fig.add_subplot(1,1,1)

def init ():
    
    ax.set_xlim (-1,1)
    ax.set_ylim (-1,1)
    
def animate():
    plot = ax.clear()
    init()
    
    Tierra = plt.Circle((0,0), rT/d, lable='Tierra')
    Luna = plt.Circle((xl[i],yl[i]), rL/d, lable='Luna')
    
    plot = ax.add_patch(Tierra)
    plot = ax.add_patch(Luna)
    
    plot = ax.scatter(x[i],y[i], lable='Rocket')
    
    plto = plt.legend()
    
    return plot

Animatio = animation.FuncAnimation(fig1, animate, )
