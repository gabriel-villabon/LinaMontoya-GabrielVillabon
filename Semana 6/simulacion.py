import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

n = 51
x = np.linspace(0.,1.,n)
y = x[:]
h = x[1]-x[0]

v = 1.
nu = 0.2