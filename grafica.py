import numpy as np
import matplotlib.pyplot as plt

E=np.genfromtxt("a.dat")
X=np.linspace(0,999,1000)

plt.figure()
plt.plot(X,E[:,0])
plt.plto(X,E[:,1])
plt.plot(X,E[:,2])
plt.savefig("energias.pdf")
plt.show()
plt.close()
