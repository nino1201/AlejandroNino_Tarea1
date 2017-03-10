import numpy as np
import matplotlib.pyplot as plt

datos=np.genfromtxt('tiempos.dat')
n=np.array([1,2,4])

plt.figure()
plt.scatter(n,datos[:3])
plt.xlabel("Procesadores")
plt.ylabel("Tiempo")
plt.savefig("tiempos.pdf")
plt.show()
