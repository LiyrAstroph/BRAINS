##
#  \file plottran.py
#  \brief plot transfer function.
#

import numpy as np
import copy as copy
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import sys

red ={"mrk335":0.0258, "mrk1044":0.0165, "iras":0.0889, "mrk382":0.0337, "mrk142":0.0449, "mcg06":0.0328,
      "irasf":0.0435,  "mrk486":0.0389,  "mrk493":0.0313, "1h0323":0.061}
def plottran(ichimin):
  print(red)
  obj = "mrk493"
  redshift = red[obj]
  ccf = np.loadtxt("../data/"+ obj +"_vr_lags.txt", skiprows=1)
  print(ccf)
  ccf[:, 0:3] = (ccf[:, 0:3]/(1.0+redshift) - 4861.0)/4861.0 * 3e5/1.0e3

  logP = np.loadtxt("../data/posterior_sample_info2d.txt", skiprows=1)
  nmax = np.argmax(logP)
  print(nmax, logP[nmax])


  fp=open("../data/tran2d_data.txt", "r")

  line=fp.readline()
  text=line.split()
  nt=int(text[1])
  nv=int(text[2])

  grid_tau=np.zeros(nt)
  grid_vel=np.zeros(nv)
  tran = np.zeros((nt, nv))
  tran_max = np.zeros((nt, nv))
  tran1d=np.zeros(nv)
  tran1d_rec=np.zeros((len(logP), nv))
  tran1d_ccf=np.zeros((len(logP), len(ccf[:, 0])))

  for j in range(0, nt):
    for i in range(0, nv):
      line=fp.readline()
      grid_vel[i], grid_tau[j], tran[j, i] = line.split()
  
    line=fp.readline()

  fp.close()
  grid_vel /= 1.0e3

  
  fp = open("../data/tran2d_rec.txt", "r")
  for i in range(len(logP)):
    for j in range(0, nt):
        line=fp.readline()
        tran[j, :] = line.split()
    
    if i == ichimin:
      tran_max = copy.copy(tran)

    fp.readline()
    for j in range(nv):
      tran1d_rec[i, j] = np.sum(tran[:, j] * grid_tau)/(np.sum(tran[:, j]) + 1.0e-10)

    for j in range(len(ccf[:, 0])):
      idx = np.where( (grid_vel < ccf[j, 1]) & (grid_vel >= ccf[j, 0]) )
      tran1d_ccf[i, j] = np.sum(tran1d_rec[i, idx[0]] * grid_tau[idx[0]])/(np.sum(tran1d_rec[i, idx[0]]) + 1.0e-10)

  fp.close()

  tran1d_rec_mean = np.mean(tran1d_rec, axis=0)
  tran1d_rec_std = np.std(tran1d_rec, axis=0)

  tran = copy.copy(tran_max)
  for i in range(0, nv):
    tran1d[i] = np.sum(tran[:, i] * grid_tau)/(np.sum(tran[:, i]) + 1.0e-10)

  figtran = plt.figure(1)

  ax = figtran.add_axes((0.1, 0.1, 0.5, 0.8))
  ax.plot(grid_vel, tran1d, lw=2, color='w')
  plt.errorbar(ccf[:, 2], ccf[:, 4], xerr=[ccf[:, 2]-ccf[:, 0], ccf[:, 1]-ccf[:, 2]], yerr=[ccf[:, 4]-ccf[:, 5], ccf[:, 6]-ccf[:, 4]], 
    ls='none', marker='o', color='r')
  ax.imshow(tran, origin="low", aspect='auto', extent=[grid_vel[0], grid_vel[-1], grid_tau[0], grid_tau[-1]])

  ax.set_xlabel('Velocity')
  ax.set_ylabel("Time Lag")

  ax = figtran.add_axes((0.62, 0.1, 0.3, 0.8))
  tran1d = np.sum(tran, axis=1)

  ax.plot(tran1d, grid_tau, color='k')
  ax.set_xlabel("Transfer Function")
  [yt.set_visible(False) for yt in ax.get_yticklabels()]

  xlim = ax.get_xlim()
  dxlim = (xlim[1] - xlim[0])/5
  ax.xaxis.set_major_locator(MultipleLocator(dxlim))

  figtran.savefig("tran.pdf", bbox_inches='tight')

  fig = plt.figure(2)
  ax1 = fig.add_axes((0.1, 0.1, 0.5, 0.8))
  
  #for i in range(0,len(logP)):
    #ax1.plot(grid_vel, tran1d_rec[i, :], color='grey')

  ax1.errorbar(ccf[:, 2], ccf[:, 4], xerr=[ccf[:, 2]-ccf[:, 0], ccf[:, 1]-ccf[:, 2]], yerr=[ccf[:, 4]-ccf[:, 5], ccf[:, 6]-ccf[:, 4]], 
    ls='none', marker='o', color='r')
  ax1.errorbar(ccf[:, 2], np.mean(tran1d_ccf, axis=0), xerr=[ccf[:, 2]-ccf[:, 0], ccf[:, 1]-ccf[:, 2]], yerr=np.std(tran1d_ccf, axis=0), 
    ls='none', marker='o', color='b')

  plt.show()

if __name__ == '__main__':
  print(sys.argv[1])
  plottran(int(sys.argv[1]))

