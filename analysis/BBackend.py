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

class Param:
  """
  load param file
  """
  def __init__(self, fname):
    self.param_file = fname
    self._param_parser(self.param_file)
  
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
  
  def set_param_file(self, fname):
    """
    set parameter file
    """
    self.param_file = fname
    self._param_parser(fname)
    return

class Options:
  """
  load options file
  """
  def __init__(self, file_dir, flag_dim):
     self._options_parser(file_dir, flag_dim)

  def _options_parser(self, file_dir, flag_dim):
    config = cp.RawConfigParser(delimiters=' ', comment_prefixes='#', inline_comment_prefixes='#', 
    default_section=cp.DEFAULTSECT, empty_lines_in_values=False)
    
    if flag_dim == '-2':
      self.option_file = ""
    elif flag_dim == '-1':
      self.option_file = file_dir + "/src/OPTIONSCON"
    elif flag_dim == '0':
      self.option_file = file_dir + "/src/OPTIONSCON"
    elif flag_dim == '1':
      self.option_file = file_dir + "/src/OPTIONS1D"
    elif flag_dim == '2':
      self.option_file = file_dir + "/src/OPTIONS2D"
    elif flag_dim == '3':
      self.option_file = file_dir + "/src/OPTIONSSA"
    elif flag_dim == '4':
      self.option_file = file_dir + "/src/OPTIONSSA1D"
    elif flag_dim == '5':
      self.option_file = file_dir + "/src/OPTIONSSA2D"

    if self.option_file != "":
      with open(self.option_file) as f:
        file_content = '[dump]\n' + f.read()
  
      config.read_string(file_content)
      self.options = config['dump']
    else:
      self.options = None


class ParaName:
  """
  load parameter names
  """
  def __init__(self, file_dir, flag_dim):
    self._load_para_names(file_dir, flag_dim)

  def _load_para_names(self, file_dir, flag_dim):

    if flag_dim == '-2':
      self.para_names_file = ""
    elif flag_dim == '-1':
      self.para_names_file = ""
    elif flag_dim == '0':
      self.para_names_file = file_dir + "/data/para_names_con.txt"
    elif flag_dim == '1':
      self.para_names_file = file_dir + "/data/para_names_1d.txt"
    elif flag_dim == '2':
      self.para_names_file = file_dir + "/data/para_names_2d.txt"
    elif flag_dim == '3':
      self.para_names_file = file_dir + "/data/para_names_sa.txt"
    elif flag_dim == '4':
      self.para_names_file = file_dir + "/data/para_names_sa1d.txt"
    elif flag_dim == '5':
      self.para_names_file = file_dir + "/data/para_names_sa2d.txt"
    else:
      raise Exception("Incorrect FlagDim.")

    names = np.genfromtxt(self.para_names_file, skip_header=7, \
      dtype=[int, 'S16', float, float, int, int, float], delimiter=[4, 16, 12, 12, 4, 5, 15])
    
    self.para_names = {}
    self.para_names['name'] = np.array([str(f.strip(), encoding='utf8') for f in names['f1']])
    self.para_names['min']  = names['f2']
    self.para_names['max']  = names['f3']
    self.para_names['prior']  = names['f4']
    self.para_names['fix']  = names['f5']
    self.para_names['val']  = names['f6']

    self.num_param_blrmodel_rm = 0
    self.num_param_blr_rm = 0
    self.num_param_rm_extra = 0
    self.num_param_sa = 0
    self.num_param_blrmodel_sa = 0
    self.num_param_sa_extra = 0
    
    if 'BLR model' in self.para_names['name']:
      self.num_param_blrmodel_rm = np.count_nonzero(self.para_names['name'] == 'BLR model')
    else:
      self.num_param_blrmodel_rm = 0
      
    if 'SA BLR model' in self.para_names['name']:
      idx = np.nonzero(self.para_names['name'] == 'SA BLR model')
      self.num_param_rm_extra = idx[0][0] - self.num_param_blrmodel_rm
      self.num_param_blr_rm = self.num_param_rm_extra + self.num_param_blrmodel_rm
      self.num_param_blrmodel_sa = np.count_nonzero(self.para_names['name'] == 'SA BLR model')
      self.num_param_sa_extra = np.count_nonzero(self.para_names['name'] ==  'SA Extra Par')
      self.num_param_sa = self.num_param_blrmodel_sa + self.num_param_sa_extra
    else:
      self.num_param_sa = 0
      self.num_param_blrmodel_sa = 0
      self.num_param_sa_extra = 0
      idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')
      self.num_param_rm_extra = idx_con[0][0] - self.num_param_blrmodel_rm
      self.num_param_blr_rm = self.num_param_rm_extra + self.num_param_blrmodel_rm

    if 'time series' in self.para_names['name']:
      idx = np.nonzero(self.para_names['name'] == 'time series')
      self.num_param_con = idx[0][0] - self.num_param_sa - self.num_param_blr_rm
    else:
      self.num_param_con = 0

    #print(self.num_param_blr_rm, self.num_param_sa, self.num_param_con)
    


class BBackend(Param, Options, ParaName):
  """
  backend for ploting BRAINS results
  """
  def __init__(self, fname="./param"):
    Param.__init__(self, fname)
    Options.__init__(self, self.param['filedir'], self.param['flagdim'])
    ParaName.__init__(self, self.param['filedir'], self.param['flagdim'])
    
    self.file_dir = self.param['filedir']+"/"
    self._load_data()
    self._load_results() 

  def _load_line2d_data(self):
    """
    load 2d line data
    """
    fp = open(self.file_dir+self.param['line2dfile'])
    line = fp.readline()
    nt = int(line[1:].split()[0])
    nv = int(line[1:].split()[1])
    line2d_time = np.empty(nt)
    line2d_profile = np.empty((nt, nv, 3))
    for i in range(nt):
      line = fp.readline()
      line2d_time[i] = float(line[1:])
      for j in range(nv):
        line = fp.readline()
        line2d_profile[i, j, :] = line.split()
      
      fp.readline()

    fp.close()
    return {"time":line2d_time, "profile":line2d_profile}
  
  def _load_sa_data(self):
    """
    load sa data
    """
    fp = open(self.file_dir+self.param['safile'])
    line = fp.readline()
    ns = int(line[1:].split()[0]) # number of line profiles
    nv = int(line[1:].split()[1]) # number of wavelengths
    nb = int(line[1:].split()[2]) # number of baseline
    sa_profiles = np.empty((ns, nv, 3))
    sa_baselines = np.empty((nb, 2))
    sa_phases = np.empty((nb, nv, 3))

    # first read line profiles
    for i in range(ns):
      for j in range(nv):
        line = fp.readline()
        sa_profiles[i, j, :] = line.split()
      
      fp.readline()

    # then read phase data
    for i in range(nb):
      line = fp.readline()
      sa_baselines[i, :] = line[1:].split()  
      for j in range(nv):
        line = fp.readline()
        sa_phases[i, j, :] = line.split()
      fp.readline()

    fp.close()
    
    return {"profile":sa_profiles, "baseline":sa_baselines, "phase":sa_phases}

  def _load_data(self):
    """
    load data
    """
    self.data={}

    if self.param['flagdim'] == '0':
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
    
    elif self.param['flagdim'] == '1':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line data
      self.data['line_data'] = np.loadtxt(self.file_dir+self.param['linefile'])

    elif self.param['flagdim'] == '2':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line2d data
      self.data['line2d_data'] = self._load_line2d_data()
    
    elif self.param['flagdim'] == '3':
      self.data["sa_data"] = self._load_sa_data()

    elif self.param['flagdim'] == '4':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line data
      self.data['line_data'] = np.loadtxt(self.file_dir+self.param['linefile'])
      # load sa data
      self.data["sa_data"] = self._load_sa_data()

    elif self.param['flagdim'] == '5':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line2d data
      self.data['line2d_data'] = self._load_line2d_data()
      # load sa data
      self.data["sa_data"] = self._load_sa_data()
    
    else:
      raise Exception("Incorrect FlagDim.")

    print("data dict:", self.data.keys())

    return
  
  def _load_con_rec(self):
    self.results['con_rec'] = np.loadtxt(self.file_dir+"data/con_rec.txt")
    ns = int(self.sample_size)
    nt = int(self.results['con_rec'].shape[0]/ns)
    nc = self.results['con_rec'].shape[1]
    self.results['con_rec'] = np.reshape(self.results['con_rec'], (ns, nt, nc), order='C')
  
  def _load_line_rec(self):
    self.results['line_rec'] = np.loadtxt(self.file_dir+"data/line_rec.txt")
    ns = int(self.sample_size)  
    nt = int(self.results['line_rec'].shape[0]/ns)
    nc = self.results['line_rec'].shape[1]
    self.results['line_rec'] = np.reshape(self.results['line_rec'], (ns, nt, nc), order='C')

  def _load_sa_rec(self):
    # load sa line reconstructions
    self.results["sa_profile_rec"] = np.loadtxt(self.file_dir + 'data/sa_line_rec.txt')
    self.results["sa_phase_rec"] = np.loadtxt(self.file_dir + 'data/sa_phase_rec.txt')
    ns = int(self.sample_size)
    nv = int(self.results['sa_profile_rec'].shape[0]/ns)
    nc = self.results['sa_profile_rec'].shape[1]
    self.results["sa_profile_rec"] = np.reshape(self.results["sa_profile_rec"], (ns, nv, nc), order='C')
    nb = int(self.results['sa_phase_rec'].shape[0]/ns/nv)
    nc = self.results['sa_phase_rec'].shape[1]
    self.results["sa_phase_rec"] = np.reshape(self.results["sa_phase_rec"], (ns, nb, nv, nc), order='C')
  
  def _load_line2d_rec(self):
    # load line2d reconstructions
    self.results['line2d_rec'] = np.loadtxt(self.file_dir+'data/line2d_rec.txt')
    ns = int(self.sample_size)
    nt = int(self.results['line2d_rec'].shape[0]/ns)
    nv = int(self.results['line2d_rec'].shape[1])
    self.results['line2d_rec'] = np.reshape(self.results['line2d_rec'], (ns, nt, nv), order='C')
 
  def _load_tran_rec(self):
    self.results['tran_rec'] = np.loadtxt(self.file_dir+"data/tran_rec.txt")
    ns = int(self.sample_size)  
    nt = int(self.results['tran_rec'].shape[0]/ns)
    nc = self.results['tran_rec'].shape[1]
    self.results['tran_rec'] = np.reshape(self.results['tran_rec'], (ns, nt, nc), order='C')

  def _load_tran2d_rec(self):
    # load tran2d_rec
    self.results['tran2d_rec'] = np.loadtxt(self.file_dir+'data/tran2d_rec.txt')
    ns = int(self.sample_size)
    ntau = int(self.results['tran2d_rec'].shape[0]/ns)
    nv = int(self.results['tran2d_rec'].shape[1])
    self.results['tran2d_rec'] = np.reshape(self.results['tran2d_rec'], (ns, ntau, nv), order='C')
    self.results['tau_rec'] = self.results['tran2d_rec'][:, :, 0]
    self.results['tran2d_rec'] = self.results['tran2d_rec'][:, :, 1:]

  def _load_results(self):
    """
    load results
    """
    self.results={}
    
    if self.param['flagdim'] == '-2':
      self.results['con_sim'] = np.loadtxt(self.file_dir+"data/sim_con.txt")
      self.results['line_sim'] = np.loadtxt(self.file_dir+"data/sim_hb.txt")
      self.results['line2d_sim'] = np.loadtxt(self.file_dir+"data/sim_hb2d.txt")

    elif self.param['flagdim'] == '-1':
      self.results['con_sim'] = np.loadtxt(self.file_dir+"data/sim_con.txt")
      self.results['line_sim'] = np.loadtxt(self.file_dir+"data/sim_hb.txt")
      self.results['line2d_sim'] = np.loadtxt(self.file_dir+"data/sim_hb2d.txt")
    
    elif self.param['flagdim'] == '0':
      self.results['sample'] = np.loadtxt(self.file_dir + "data/posterior_sample_con.txt")
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
    
    elif self.param['flagdim'] == '1':
      self.results['sample'] = np.loadtxt(self.file_dir + "data/posterior_sample_1d.txt")      
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
      self._load_line_rec()
      self._load_tran_rec()
      
    elif self.param['flagdim'] == '2':
      # load posterior samples
      self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_2d.txt')
      # load likelihoods
      self.results['sample_info'] = np.loadtxt(self.file_dir + 'data/posterior_sample_info_2d.txt')
      self.sample_size = self.results['sample'].shape[0]
      
      # load continuum reconstructions
      self._load_con_rec()
      self._load_line2d_rec()
      self._load_tran2d_rec()


    elif self.param['flagdim'] == '3':
      self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa.txt')
      self.sample_size = self.results['sample'].shape[0]
      self._load_sa_rec()

    
    elif self.param['flagdim'] == '4':
      self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa1d.txt')
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
      self._load_line_rec()
      self._load_tran_rec()
      self._load_sa_rec()

    elif self.param['flagdim'] == '5':
      self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa2d.txt')
      self.sample_size = self.results['sample'].shape[0]
      self._load_line2d_rec()
      self._load_tran2d_rec()
      self._load_sa_rec()

    else:
      raise Exception("Incorrect FlagDim.")  
    
    print("results dict:", self.results.keys())
    return
  
  
  def get_con_data(self):
    return self.data['con_data']
  def get_line_data(self):
    return self.data['line_data']
  def get_line2d_data(self):
    return self.data['line2d_data']
  def get_sa_data(self):
    return self.data['sa_data']

  def plot_drw_parameters(self):
    idx = self.num_param_blr_rm+self.num_param_sa
    fig = corner.corner(self.results['sample'][:, idx:idx+3], labels=['syserr', 'ln sigma', 'ln tau'], smooth=True)
    return fig

