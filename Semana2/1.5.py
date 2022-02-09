import numpy as np
import matplotlib.pyplot as plt
import scipy.io.wavfile as wav
import os.path as path
import wget

data = np.loadtxt('ManchasSolares.dat.txt')

#a) Filtrar datos

mask1= data[:,0] >=1900
data = data[mask1]
print(data)

#b) Quitar valor medio (centrar los datos)
mean = np.mean(data[:,3])
spots = data[:,3] - mean 
print(spots)

#c) fft y fftfreq

fft = np.fft.fft(spots)

Npoints = len(spots)
Frequencies = np.fft.fftfreq(Npoints,1)
print(fft)
print(Frequencies)

#d) Encontrar la frecuencia dominante por año
ffta = np.abs(fft)

mask2 = Frequencies < 0 
fft[mask2] = 0

NewFreq = np.arange(len(fft))

i = np.argmax(ffta)
FreqMax = Frequencies[i]
print(Frequencies[i])

#e) Periodo por año
T = 1/(Frequencies[i]*12)
print (T)

#f) Figura
x = data[:,0] + data[:,1]/12

y=np.cos((1/T)*x)

plt.figure()
plt.plot(x,spots, label='Datos')
plt.plot(x, y, label='Frecuencia dominante', color='r')
plt.title('Periodo 10.67 en años')
plt.legend()
plt.show()
plt.savefig('ManchasSolares.png')