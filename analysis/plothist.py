##
# \file plothist.py
# \brief plot histogram of the posterior sample.
#        this procedure plots histogram of the posterior sample 
#        according to the given dimension
#

import numpy as np
import corner
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages


def plothist(ndim):
  pdf = PdfPages('hist.pdf')
  if(ndim == 0):
    samples = np.loadtxt("../data/posterior_sample.txt", skiprows=1)
  elif(ndim == 1):
    samples = np.loadtxt("../data/posterior_sample1d.txt", skiprows=1)
  elif(ndim == 2):
    samples = np.loadtxt("../data/posterior_sample2d.txt", skiprows=1)
  else:
    print("Incorrect dimension.")
    exit(0)

  for i in range(samples.shape[1]):
    plt.hist(samples[:, i])
    plt.axvline(x=np.mean(samples[:, i]))
    pdf.savefig()
    plt.close()
  
  pdf.close()


if __name__ == "__main__":
  if(len(sys.argv) < 2):
    print("No dimension specified.")
    exit(0);

  plothist(int(sys.argv[1]))

