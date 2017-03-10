import numpy as np
import matplotlib.pyplot as plt

E=np.genfromtxt("datos.dat")
X=np.linspace(0,1000,len(E[:1000,0]))

plt.figure()
plt.plot(X,E[:1000,0], label="Modo 1")
plt.plot(X,E[:1000,1], label="Modo 2")
plt.plot(X,E[:1000,2], label="Modo 3")
plt.legend()
plt.title("Energias")
plt.ylabel("E")
plt.xlabel("Tiempo")
plt.savefig("energias.pdf")
plt.show()
plt.close()
