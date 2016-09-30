import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter

nt = 100
nv = 100

fp=open("../data/tran2d.txt", "r")

date=np.zeros(nt)
grid_vel=np.zeros(nv)
tran = np.zeros((nt, nv))

for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel[i], date[j], tran[j, i] = line.split()
  
  line=fp.readline()

fig = plt.figure()

plt.imshow(tran, origin="low")

fig = plt.figure()

tran1d = np.sum(tran, axis=1)

plt.plot(date, tran1d)

plt.show()
