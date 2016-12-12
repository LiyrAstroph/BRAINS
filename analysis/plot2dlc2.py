##
# \file plot2dlc.py
# \brief plot light curves of emission line for 2d analysis.
#

import numpy as np 
import matplotlib.pyplot as plt 

VelUnit = np.sqrt( 6.672e-8 * 1.0e6 * 1.989e33 / (2.9979e10*8.64e4)) / 1.0e5

obj = "mrk493"

fp=open("../data/sim_hb2d.txt", "r")
line=fp.readline()
text=line.split()
nt=int(text[1])
nv=int(text[2])

print(nv, nt)

date_hb=np.zeros(nt)
grid_vel=np.zeros(nv)
prof = np.zeros((nt, nv))
prof_err = np.zeros((nt, nv))

for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel[i], date_hb[j], prof[j, i], prof_err[j, i] = line.split()
  
  line=fp.readline()

fp.close()

lc = np.sum(prof, axis=1)
plt.plot(date_hb, lc, marker='o')
plt.show()


