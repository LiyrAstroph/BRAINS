##
# \file plot2dlc.py
# \brief plot light curves of emission line for 2d analysis.
#

import numpy as np 
import matplotlib.pyplot as plt 

VelUnit = np.sqrt( 6.672e-8 * 1.0e6 * 1.989e33 / (2.9979e10*8.64e4)) / 1.0e5

obj = "mrk493"

fp=open("../data/pline2d_data.txt", "r")
nv = 25
nt = 25

grid_tau=np.zeros(nt)
grid_vel=np.zeros(nv)
tran2d = np.zeros((nt, nv))

for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel[i], grid_tau[j], tran2d[j, i] = line.split()
  
  line=fp.readline()

fp.close()

tran1d = np.sum(tran2d, axis= 1)

plt.plot(grid_tau, tran1d)

plt.show()


