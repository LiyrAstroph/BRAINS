# -*- coding: utf-8 -*-

#
# BRAINS
# (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
# Yan-Rong Li, liyanrong@ihep.ac.cn
# Thu, Aug 4, 2016
#

import os
import sys
import corner
import numpy as np 
import configparser as cp
import matplotlib.pyplot as plt

__all__ = ['plotbackend']

class plotbackend:
  """
  plot backend for BRAINS
  """
  def __init__(self, fname="../src/param", fopt=""):
    
    # setup param
    self.param_file = fname
    self._param_parser(self.param_file)
    self.file_dir = self.param['filedir']+"/"

    # setup options
    if fopt == "":
      if self.param['flagdim'] == '-2':
        self.option_file = ""
      elif self.param['flagdim'] == '-1':
        self.option_file = self.file_dir + "/src/OPTIONSCON"
      elif self.param['flagdim'] == '0':
        self.option_file = self.file_dir + "/src/OPTIONSCON"
      elif self.param['flagdim'] == '1':
        self.option_file = self.file_dir + "/src/OPTIONS1D"
      elif self.param['flagdim'] == '2':
        self.option_file = self.file_dir + "/src/OPTIONS2D"
      elif self.param['flagdim'] == '3':
        self.option_file = self.file_dir + "/src/OPTIONSSA"
      elif self.param['flagdim'] == '4':
        self.option_file = self.file_dir + "/src/OPTIONSSA1D"
      elif self.param['flagdim'] == '5':
        self.option_file = self.file_dir + "/src/OPTIONSSA2D"
      
      if self.option_file != "":
        self._option_load(self.option_file)
        self._get_sample_size()
      else:
        self.sample_size = 0
    
    else:
      self.set_option_file(fopt)

    

  def _option_load(self, fname):
    """
    load option file
    """
    with open(fname, "r") as f:
      lines = f.readlines()
    
    # negect comments
    i=0
    for line in lines:
      if line[0] == '#' or len(line.strip()) == 0:
        i+=1

    option={} 

    option['num_particles'] = int(lines[i].split()[0])
    i+=1
    option['new_level_interval'] = int(lines[i].split()[0])
    i+=1
    option['save_interval'] = int(lines[i].split()[0])
    i+=1
    option['thread_step'] = int(lines[i].split()[0])
    i+=1
    option['num_levels'] = int(lines[i].split()[0])
    i+=1
    option['lambda'] = int(lines[i].split()[0])
    i+=1
    option['beta'] = int(lines[i].split()[0])
    i+=1
    option['num_saves'] = int(lines[i].split()[0])
    i+=1
    option['file_sample'] = lines[i].split()[0]
    i+=1
    option['file_sample_info'] = lines[i].split()[0]
    i+=1
    option['file_levels'] = lines[i].split()[0]
    i+=1
    option['file_sampler_state'] = lines[i].split()[0]
    i+=1
    option['file_post_sample'] = lines[i].split()[0]
    i+=1
    option['file_post_sample_info'] = lines[i].split()[0]
    i+=1
    option['file_limits'] = lines[i].split()[0]
  
    self.option = option

  def _param_parser(self, fname):
    """
    parse parameter file
    """
    config = cp.RawConfigParser(delimiters=' ', comment_prefixes='%', inline_comment_prefixes='%', 
    default_section=cp.DEFAULTSECT, empty_lines_in_values=False)
  
    with open(fname) as f:
      file_content = '[dump]\n' + f.read()

    config.read_string(file_content)
    
    # check the absolute path
    if os.path.isabs(config['dump']['filedir']) == False:
      raise Exception("FileDir in %s is not an absoulte path.\n"%self.param_file)
    
    self.param = config['dump']
  
  def _get_sample_size(self):
    """
    load results
    """
    with open(self.file_dir+self.option['file_post_sample']) as f:
      self.sample_size = int(f.readline().split()[1])
    
  
  def set_param_file(self, fname):
    """
    set parameter file
    """
    self.param_file = fname
    self._param_parser(fname)
    self.file_dir = self.param['filedir']+"/"
    return
  
  def set_option_file(self, fname):
    """
    set option file
    """
    self.option_file = fname
    self._option_load(fname)
    self._get_sample_size()
    return

  def load_results(self):
    """
    load results
    """
    self.results={}
    
    if self.param['flagdim'] == '-2':
      self.results['con_sim'] = np.loadtxt(self.file_dir+"data/sim_con.txt")
      self.results['line_sim'] = np.loadtxt(self.file_dir+"data/sim_hb.txt")
      self.results['line2d_sim'] = np.loadtxt(self.file_dir+"data/sim_hb2d.txt")

    elif self.param['flagdim'] == '-1':
      self.results['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      self.results['con_sim'] = np.loadtxt(self.file_dir+"data/sim_con.txt")
      self.results['line_sim'] = np.loadtxt(self.file_dir+"data/sim_hb.txt")
      self.results['line2d_sim'] = np.loadtxt(self.file_dir+"data/sim_hb2d.txt")
    
    elif self.param['flagdim'] == '0':
      self.results['sample'] = np.loadtxt(self.file_dir + self.option['file_post_sample'])
      self.results['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      self.results['con_rec'] = np.loadtxt(self.file_dir+"data/con_rec.txt")
      
    return
  
  def plot_drw_parameters(self):
    if self.param['flagdim'] == '3' or self.param['flagdim'] == '-2':
      raise Exception("FlagDim=%d, no DRW parameters.\n"%self.param['flagdim'])
    
    sample = self.results['sample']
    fig = corner.corner(sample[:, 1:3], smooth=True, smooth1d=True, labels=[r"$\ln(\hat\sigma)$", r"$\ln(\tau)$"])

    return fig
    

  def plot_con_rec(self):

    if self.param['flagdim'] == '3' or self.param['flagdim'] == '-2':
      raise Exception("FlagDim=%d, no continuum reconstruction.\n"%self.param['flagdim'])

    con_data = self.results['con_data']
    con = self.results['con_rec']
    offset = int(con.shape[0]/self.sample_size)

    fig, ax = plt.subplots(1,1)
    ax.errorbar(con_data[:, 0], con_data[:, 1], yerr = con_data[:, 2], ls='none', marker='o', label='Data')
    for i in np.random.randint(self.sample_size, size=np.min((100, self.sample_size))):
      plt.plot(con[i*offset:(i+1)*offset, 0], con[i*offset:(i+1)*offset, 1], lw=0.2)
    
    ax.set_xlabel('Time')
    ax.set_ylabel('Flux')
    ax.legend()

    return fig
