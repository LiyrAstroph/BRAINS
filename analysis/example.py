import matplotlib.pyplot as plt 

from plotbackend import plotbackend 

pb = plotbackend("../src/param")

pb.load_results()

fig = pb.plot_drw_parameters()
plt.show()