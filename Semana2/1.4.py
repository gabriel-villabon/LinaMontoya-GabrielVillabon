import numpy as np
import matplotlib.pyplot as plt

t=np.linspace(-2*np.pi,2*np.pi,100)
h = 100/(4*np.pi)

def function(t):
    f = np.exp(-0.1*t)*np.sin(t)
    return f

def central_Derivate(f,t,h):
    d = 0
    if h != 0:
        d = (f(t+h)-f(t-h))/(2*h)
    return d

def right_Derivate(f,t,h):
    d = 0
    if h != 0:
        d = (f(t+h)-f(t))/h
    return d

def exact_Derivate(t):
    f = -0.1*np.exp(-0.1*t)*np.sin(t) + np.exp(-0.1*t)*np.cos(t)
    return f

X = np.fft.fft(function(t))
centrals_Derivate = central_Derivate(function,t,h)
rights_Derivate = right_Derivate(function,t,h)
exacts_Derivate = (exact_Derivate(t))/10

frequecy =np.fft.fftfreq(np.size(t),np.max(exacts_Derivate))
X *=1j*frequecy
spectral_Derivate = np.fft.ifft(X)

plt.figure()
plt.plot(t,spectral_Derivate,color='g',label='Spectral Derivate')
plt.plot(t,rights_Derivate,color='b',label='Right Derivate')
plt.plot(t,exacts_Derivate,color='r',label='Exact Derivate')
plt.plot(t,centrals_Derivate,color='k',label='Central Derivate')
plt.legend()
plt.show()
plt.savefig('Derivative superposition.png')