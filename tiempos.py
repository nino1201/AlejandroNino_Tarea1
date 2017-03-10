import numpy as np
import matplotlib.pyplot as plt

datos=np.genfromtxt('tiempos.dat')
plt.scatter(datos[:,0],datos[:,1])
plt.safe("tiempos.pdf")
