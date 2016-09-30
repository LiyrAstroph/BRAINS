import numpy as np
import corner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


pdf = PdfPages('hist.pdf')

samples = np.loadtxt("../data/posterior_sample1d.txt", skiprows=1)


for i in range(samples.shape[1]):
  plt.hist(samples[:, i])
  plt.axvline(x=np.mean(samples[:, i]))
  pdf.savefig()
  plt.close()

pdf.close()
