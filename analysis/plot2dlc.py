##
# \file plot2dlc.py
# \brief plot light curves of emission line for 2d analysis.
#

import numpy as np 
import matplotlib.pyplot as plt 

VelUnit = np.sqrt( 6.672e-8 * 1.0e6 * 1.989e33 / (2.9979e10*8.64e4)) / 1.0e5

obj = "mrk493"

fp=open("../data/"+ obj + "_hb2d.txt", "r")
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

grid_wave = grid_vel / 3e5 * 4861.0 + 4861.0
dV = (grid_vel[1]-grid_vel[0]) / VelUnit

logP = np.loadtxt("../data/posterior_sample_info_2d.txt", skiprows=1)
nmax = np.argmax(logP)
print(nmax, logP[nmax])

prof_rec = np.zeros((nt, nv))
prof_rec_max = np.zeros((nt, nv))
line_rec = np.zeros(nt)
line_rec_max = np.zeros(nt)
fp = open("line2d_rec.txt", "r")

for i in range(len(logP)):
  for j in range(nt):
    line = fp.readline()
    prof_rec[j, :]=line.split()
  
  fp.readline()
  line_rec = np.sum(prof_rec, axis=1)
  plt.plot(grid_vel, line_rec, color='grey', alpha=0.01)

  if i == 500:
    plt.plot(grid_vel, line_rec, color='g', lw = 2)

fp.close()
plt.show()


