import matplotlib.pyplot as plt 

from plotbackend import plotbackend 

pb = plotbackend("../param/param")

pb.load_results()

fig = pb.plot_drw_parameters()
plt.show()