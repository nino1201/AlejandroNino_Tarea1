import numpy as np
import matplotlib.pyplot as plt

datos=np.genfromtxt('tiempos.dat')
plt.scatter(datos[:,1],datos[:,0])
plt.savefig("tiempos.pdf")
