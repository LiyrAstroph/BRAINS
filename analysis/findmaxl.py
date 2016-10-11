import numpy as np

info = np.loadtxt("../data/sample_info2d.txt", skiprows=1)
sample = np.loadtxt("../data/sample2d.txt", skiprows=1)

idx = np.argmax(info[:, 1])

print(idx, info[idx, 1])
print(sample[idx, :])
