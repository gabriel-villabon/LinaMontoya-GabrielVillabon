import numpy as np
from celluloid import Camera
import matplotlib.pyplot as plt

def function(time,armonic):
    f = (2/armonic)*(-1)**(armonic-1)*np.sin(armonic*time)
    return f

L = np.pi
cicles = 2
time = np.linspace(0,2*L*cicles,100)

armonics_init = 1
armonics_end = 50
Function = 0

figure = plt.figure(figsize=(6,6))
camera = Camera(figure)

while armonics_init < armonics_end+1:

    Function += function(time,armonics_init)
    plt.plot(time,Function,c="g")
    camera.snap()
    armonics_init += 1

animation = camera.animate()
animation.save("fourierFunction.gif")