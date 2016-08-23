import numpy as np
import matplotlib.pyplot as plt

sample = np.loadtxt("../data/sample2d.txt", skiprows=1)
sample_info = np.loadtxt("../data/sample_info2d.txt", skiprows=1)

plt.plot(sample_info[:, 0], sample[:, 8]/np.log(10.0)+6.0, ls='none', marker='o')

plt.show()
