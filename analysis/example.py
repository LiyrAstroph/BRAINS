import matplotlib.pyplot as plt 

from BBackend import BBackend 

pb = BBackend("../param/param")

# plot continuum DRW parameters
fig = pb.plot_drw_parameters()

# plot RM 2D results
# 
# pb.plot_results_2d_style2018()
pb.plot_results_2d_style2022()

