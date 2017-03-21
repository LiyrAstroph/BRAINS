# check Appendix A in Pancoast et al. 2014a

import numpy as np
import matplotlib.pyplot as plt

Ec = -0.5

Lc = 1.0

lam = 1.5
n = 10000
chi = np.random.randn(n) * lam

E =(chi+ 1.0) * Ec

#plt.hist(E)

Lmax = np.sqrt(2.0 * (E + 1.0))
L = Lmax * lam * np.log( (np.exp(1.0/lam) - 1.0) * np.random.random(n)  + 1.0)

vr2 = 2.0*(E + 1.0) - L*L
u = np.random.rand(n)
sig = [1 if ut > 0.5 else -1 for ut in u]

idx = np.where(vr2 < 0)
vr2[idx[0]] = 0
vr = np.sqrt(vr2) * sig

u = np.random.rand(n)
sig = [1 if ut > 0.5 else -1 for ut in u]
vphi = L * sig


plt.plot(vr, vphi, ls='none', marker='.', markersize=1)

r= 1.0
theta = np.linspace(0, 2.0*np.pi, 100)
x = r * np.cos(theta)
y = r * np.sin(theta)

plt.plot(x, y)

plt.show()
