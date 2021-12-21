import BBackend as bb 
import matplotlib.pyplot as plt 

bplot = bb.Bbackend()

bplot.load_results()

line2d_rec = bplot.results["line2d_rec"]

for i in range(line2d_rec.shape[0]):
    plt.plot(line2d_rec[i, 0, :])

plt.show()