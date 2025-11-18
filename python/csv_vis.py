import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt("output/pressure_final.csv", delimiter=",")
plt.imshow(p, cmap="seismic", origin="lower")
plt.colorbar()
plt.show()
