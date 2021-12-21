##
# \file plotlimits.py
# \brief plot the limits of parameters at each level.
#

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import sys

def plotlimits(ndim):
  if ndim == 0:
    str_dim = '_con'
    npar = 4+1
  elif ndim == 1:
    str_dim = '_1d'
    npar =  8+4+1
  elif ndim == 2:
    str_dim = '_2d'
    npar = 12+4+1
  else:
    print('incorrect dimension.')
    exit()

  pdf = PdfPages("limits"+str_dim+".pdf")

  limits = np.loadtxt("../data/limits"+str_dim+".txt", comments='#')

  for i in range(npar):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(limits[:, 0], limits[:, (i+1)*2-1], marker='o')
    ax.plot(limits[:, 0], limits[:, (i+1)*2], marker='o')
    ax.set_xlabel("Level")
    ax.set_ylabel("Limits")
    pdf.savefig()
    plt.close()

  pdf.close()



if __name__ == '__main__':

  if(len(sys.argv) < 2):
    print("No dimension specified.")
    exit(0)

  plotlimits(int(sys.argv[1]))
