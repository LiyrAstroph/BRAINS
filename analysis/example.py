# -*- coding: utf-8 -*-

#
# BRAINS
# (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
# Yan-Rong Li, liyanrong@ihep.ac.cn
# Thu, Aug 4, 2016
#
# an example script for plotting BRAINS results using the backend bplotlib

import matplotlib.pyplot as plt 
import os.path as path
from bbackend import bplotlib 

# load the param file
pb = bplotlib("../param/param")

# get flagdim and file directory
flagdim = int(pb.param["flagdim"])
fdir = pb.param['filedir']

#===============================================
# continuum results, flagdim == 0
if flagdim == 0:
  # plot continuum DRW parameters
  pb.plot_drw_parameters(doshow=False)
  
  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature, doshow=False)

  pb.plot_results_con()

#===============================================
# RM 1D results, flagdim == 1
if flagdim == 1:
  # plot continuum DRW parameters
  pb.plot_drw_parameters(doshow=False)
  
  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature, doshow=False)

  pb.plot_results_1d(doshow=False)

  # plot 1d transfer function
  # pb.plot_tran1d(tau_range=[a, b])
  pb.plot_tran1d(doshow=False)

  # plot histograms of BLR model parameters
  pb.plot_blrmodel_para_hist(doshow=False)

  # plot clouds' distribution
  pb.plot_clouds(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=False, doshow=False)

  # plot clouds' distribution viewd from line of sight
  pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=False, doshow=False)


#===============================================
# RM 2D results, flagdim == 2
# 
if flagdim == 2:
  # print parameter names
  pb.print_blrmodel_para_names()
  
  # plot continuum DRW parameters
  pb.plot_drw_parameters(doshow=False)
  
  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature, doshow=False)
  
  # get continuum data 
  con_data = pb.get_con_data()
  
  # get line 2D data 
  line2d_data = pb.get_line2d_data()
  
  # pb.plot_results_2d_style2018(doshow=False)
  pb.plot_results_2d_style2022(doshow=False)
  
  # plot 2d tranfer function with the maximum prob
  # pb.plot_tran2d(tau_range=[a, b])
  # pb.plot_tran2d(vel_range=[a, b])
  # pb.plot_tran2d(tau_range=[a, b], vel_range=[c, d], doshow=False)
  pb.plot_tran2d(doshow=False)

  # plot mean and rms spectra and their reconstruction
  # ci: confidence interval
  pb.plot_mean_rms(doshow=False, ci=0.95)
  
  # plot 1d transfer function
  # pb.plot_tran1d(tau_range=[a, b], doshow=False)
  pb.plot_tran1d(doshow=False)
  
  # plot histograms of BLR model parameters
  pb.plot_blrmodel_para_hist(doshow=False)

  # plot clouds' distribution
  pb.plot_clouds(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)
  
  
  # plot clouds' distribution viewd from line of sight
  pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)

#===============================================
# Line profile fitting, flagdim == 3
# 
if flagdim == 3:

  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature, doshow=False)

  # plot line profile fitting results
  pb.plot_results_lp(doshow=False)

  # plot histograms of BLR model parameters
  pb.plot_blrmodel_para_hist(doshow=False)

  # plot clouds' distribution
  pb.plot_clouds(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)
  
  # plot clouds' distribution viewd from line of sight
  pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), doshow=False)
  #pb.plot_clouds_los(path.join(fdir, "data/clouds.txt"), range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)


#===============================================
# RM SA results, flagdim == 4
# 
if flagdim == 4:

  # plot CDNest diagnoistics
  temperature = 1
  pb.postprocess(temperature, doshow=False)
  
  # plot SA results
  pb.plot_results_sa(show_offset=True, subtract_offset=True, phase_limit=[-0.9, 0.9])
  # full arguments are:
  # pb.plot_results_sa(show_offset=True, subtract_offset=True, phase_limit=[-0.9, 0.9], column_first=False, average_baseline=3)
  # column_first: when plotting multiple columns (baselines>6), the baselines are arranged by column
  # average_baseline: the number of baselines to average, counting from the baseline with the largest blr signal.
  