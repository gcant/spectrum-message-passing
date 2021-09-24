import numpy as np
import matplotlib.pyplot as plt
import os

# compute spectral density with broadening parameter eta, from
# the true eigenvalues, l
def rho(x, eta, l):
    z = x+(0+1j)*eta
    return np.imag(-np.mean(1/(z-l))/np.pi)

# eigenvalues computed directly
l = np.loadtxt("eigenvalues.txt")

## run message passing for r=0,1,2
os.system("../spectrum  laplacian.txt  0 0.05 0 10 100 > out0.txt")
os.system("../spectrum  laplacian.txt  1 0.05 0 10 100 > out1.txt")
os.system("../spectrum  laplacian.txt  2 0.05 0 10 100 > out2.txt")

X0 = np.loadtxt("out0.txt")
X1 = np.loadtxt("out1.txt")
X2 = np.loadtxt("out2.txt")
y = np.array([rho(x,0.05,l) for x in X0[:,0]])


fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(3,7))

ax1.plot(X0[:,0], X0[:,1], color='black')
ax1.set_xlim(0,10)
ax1.set_ylim(0,1)
ax1.fill_between(X0[:,0],0,y,alpha=1.0)
ax1.set_ylabel(r'$\rho(x)$, $r=0$')

ax2.plot(X1[:,0], X1[:,1], color='black')
ax2.set_xlim(0,10)
ax2.set_ylim(0,1)
ax2.fill_between(X0[:,0],0,y,alpha=1.0)
ax2.set_ylabel(r'$\rho(x)$, $r=1$')

ax3.plot(X2[:,0], X2[:,1], color='black')
ax3.set_xlim(0,10)
ax3.set_ylim(0,1)
ax3.fill_between(X0[:,0],0,y,alpha=1.0)
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$\rho(x)$, $r=2$')

plt.tight_layout()
plt.savefig('figure1.pdf')

plt.show()
