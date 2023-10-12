import matplotlib.pyplot as plt 

from bbackend import bplotlib 

pb = bplotlib("../param/param")
flagdim = int(pb.param["flagdim"])

#===============================================
# continuum results, flagdim == 0
if flagdim == 0:
  pb.plot_results_con()

#===============================================
# RM 2D results, flagdim == 2
# 
if flagdim == 2:
  # print parameter names
  pb.print_blrmodel_para_names()
  
  # plot continuum DRW parameters
  fig = pb.plot_drw_parameters()
  
  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature)
  
  # get continuum data 
  con_data = pb.get_con_data()
  
  # get line 2D data 
  line2d_data = pb.get_line2d_data()
  
  # pb.plot_results_2d_style2018()
  pb.plot_results_2d_style2022()
  
  # plot 2d tranfer function with the maximum prob
  # pb.plot_tran2d(tau_range=[a, b])
  # pb.plot_tran2d(vel_range=[a, b])
  # pb.plot_tran2d(tau_range=[a, b], vel_range=[c, d])
  pb.plot_tran2d()
  
  # plot 1d transfer function
  # pb.plot_tran1d(tau_range=[a, b])
  pb.plot_tran1d()
  
  # plot histograms of BLR model parameters
  pb.plot_blrmodel_para_hist()