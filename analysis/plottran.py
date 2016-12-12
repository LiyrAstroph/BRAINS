##
#  \file plottran.py
#  \brief plot transfer function.
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter

#redshift = 0.0165
#ccf = np.loadtxt("../data/1044_vr_lags.txt", skiprows=1)
#print(ccf)
#ccf[:, 0:3] = (ccf[:, 0:3]/(1.0+redshift) - 4861.0)/4861.0 * 3e5/1.0e3

nt = 200
nv = 25

fp=open("../data/tran2d.txt", "r")

grid_tau=np.zeros(nt)
grid_vel=np.zeros(nv)
tran = np.zeros((nt, nv))
tran1d=np.zeros(nv)
fp.readline()
for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel[i], grid_tau[j], tran[j, i] = line.split()
  
  line=fp.readline()

grid_vel /= 1.0e3

for i in range(0, nv):
  tran1d[i] = np.sum(tran[:, i] * grid_tau)/(np.sum(tran[:, i])+1.0e-10)

print(tran1d)
fig = plt.figure()
ax = fig.add_axes((0.1, 0.1, 0.5, 0.8))
plt.plot(grid_vel, tran1d, lw=2, color='w')
#plt.errorbar(ccf[:, 2], ccf[:, 4], xerr=[ccf[:, 2]-ccf[:, 0], ccf[:, 1]-ccf[:, 2]], yerr=[ccf[:, 4]-ccf[:, 5], ccf[:, 6]-ccf[:, 4]], ls='none', marker='o', color='r')
plt.imshow(tran, origin="low", aspect='auto', extent=[grid_vel[0], grid_vel[-1], grid_tau[0], grid_tau[-1]])

ax.set_xlabel('Velocity')
ax.set_ylabel("Time Lag")

ax = fig.add_axes((0.62, 0.1, 0.3, 0.8))
tran1d = np.sum(tran, axis=1)

plt.plot(tran1d, grid_tau, color='k')
ax.set_xlabel("Transfer Function")
[yt.set_visible(False) for yt in ax.get_yticklabels()]

xlim = ax.get_xlim()
dxlim = (xlim[1] - xlim[0])/5
ax.xaxis.set_major_locator(MultipleLocator(dxlim))

plt.savefig("tran.pdf", bbox_inches='tight')
plt.show()
