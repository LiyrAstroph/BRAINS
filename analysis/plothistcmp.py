import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

d=np.loadtxt("../data/posterior_sample2d.txt", skiprows=1)

fig = plt.figure(figsize=(10, 10))

ax1 = fig.add_axes((0.1, 0.7, 0.2, 0.25))

ax1.plot(d[:, 4], d[:, 8]/np.log(10.0)+6.0, ls='none', marker='.', markersize=3)

[t.set_visible(False) for t in ax1.get_xticklabels()]
ax1.xaxis.set_major_locator(MultipleLocator(50))
ax1.yaxis.set_major_locator(MultipleLocator(0.5))

ax2 = fig.add_axes((0.3, 0.7, 0.2, 0.25))

ax2.plot(d[:, 3], d[:, 8]/np.log(10.0)+6.0, ls='none', marker='.', markersize=3)
[t.set_visible(False) for t in ax2.get_yticklabels()]
[t.set_visible(False) for t in ax2.get_xticklabels()]
ax2.yaxis.set_major_locator(MultipleLocator(0.5))

ax3 = fig.add_axes((0.5, 0.7, 0.2, 0.25))
ax3.hist(d[:, 8]/np.log(10.0)+6.0, bins = 10, normed=True)
[t.set_visible(False) for t in ax3.get_yticklabels()]
ax3.xaxis.set_major_locator(MultipleLocator(1.0))


ax4 = fig.add_axes((0.1, 0.45, 0.2, 0.25))
ax4.plot(d[:, 4], d[:, 3], ls='none', marker='.', markersize=3)
[t.set_visible(False) for t in ax4.get_xticklabels()]
ax4.xaxis.set_major_locator(MultipleLocator(50))
ax4.yaxis.set_major_locator(MultipleLocator(30))


ax5 = fig.add_axes((0.3, 0.45, 0.2, 0.25))
ax5.hist(d[:, 3], bins = 10, normed=True)
[t.set_visible(False) for t in ax5.get_yticklabels()]
ax5.xaxis.set_major_locator(MultipleLocator(30))


ax6 = fig.add_axes((0.1, 0.20, 0.2, 0.25))
ax6.hist(d[:, 4], bins = 10, normed=True)
[t.set_visible(False) for t in ax6.get_yticklabels()]
ax6.xaxis.set_major_locator(MultipleLocator(50))

ax7 = fig.add_axes((0.78, 0.72, 0.2, 0.22))
ax7.hist(d[:, 0], bins = 10, normed=True)

ax8 = fig.add_axes((0.78, 0.45, 0.2, 0.22))
ax8.hist(d[:, 1], bins = 10, normed=True)

ax9 = fig.add_axes((0.78, 0.2, 0.2, 0.22))
ax9.hist(d[:, 10], bins = 10, normed=True)

plt.savefig("histcmp.pdf", bbox_inches='tight')

#plt.show()
