# -*- coding: utf-8 -*-

#
# BRAINS
# (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
# Yan-Rong Li, liyanrong@ihep.ac.cn
# Thu, Aug 4, 2016
#
# a plotting backend for BRAINS

import os
import re
import sys
import copy
import corner
import numpy as np 
import configparser as cp
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from matplotlib.backends.backend_pdf import PdfPages
from .postprocess import postprocess as pt

__all__ = ['bplotlib']

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
      self.option_file = file_dir + "/param/OPTIONSCON"
    elif flag_dim == '0':
      self.option_file = file_dir + "/param/OPTIONSCON"
    elif flag_dim == '1':
      self.option_file = file_dir + "/param/OPTIONS1D"
    elif flag_dim == '2':
      self.option_file = file_dir + "/param/OPTIONS2D"
    elif flag_dim == '3':
      self.option_file = file_dir + "/param/OPTIONSLP"
    elif flag_dim == '4':
      self.option_file = file_dir + "/param/OPTIONSSA"
    elif flag_dim == '5':
      self.option_file = file_dir + "/param/OPTIONSSA1D"
    elif flag_dim == '6':
      self.option_file = file_dir + "/param/OPTIONSSA2D"
    elif flag_dim == '7':
      self.option_file = file_dir + "/param/OPTIONSSARM"

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
      self.para_names_file = file_dir + "/data/para_names_lp.txt"
    elif flag_dim == '4':
      self.para_names_file = file_dir + "/data/para_names_sa.txt"
    elif flag_dim == '5':
      self.para_names_file = file_dir + "/data/para_names_sa1d.txt"
    elif flag_dim == '6':
      self.para_names_file = file_dir + "/data/para_names_sa2d.txt"
    elif flag_dim == '7':
      self.para_names_file = file_dir + "/data/para_names_sarm.txt"
    else:
      raise Exception("Incorrect FlagDim.")

    names = np.genfromtxt(self.para_names_file, comments='#', \
      dtype=[int, 'S30', float, float, int, int, float], delimiter=[4, 30, 12, 12, 4, 5, 15])
    
    # the first line is header
    names = np.delete(names, 0)
    
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
    
    self.num_param_blrmodel_rm = 0
    self.num_param_blrmodel_sa = 0
    self.num_param_sa_extra = 0
    for name in self.para_names['name']:
      if re.match("BLR_model", name):
        self.num_param_blrmodel_rm += 1
      
      if re.match("SA_BLR_model", name):
        self.num_param_blrmodel_sa += 1
      
      if re.match("SA_Extra_Par", name):
        self.num_param_sa_extra += 1
    
    self.num_param_sa = self.num_param_blrmodel_sa + self.num_param_sa_extra
    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')
    if len(idx_con[0]) > 0:
      self.num_param_rm_extra = idx_con[0][0] - self.num_param_blrmodel_rm - self.num_param_sa
    else:
      self.num_param_rm_extra = 0

    self.num_param_blr_rm = self.num_param_rm_extra + self.num_param_blrmodel_rm

    if 'time_series' in self.para_names['name']:
      idx = np.nonzero(self.para_names['name'] == 'time_series')
      self.num_param_con = idx[0][0] - self.num_param_sa - self.num_param_blr_rm
    else:
      self.num_param_con = 0

    #print(self.num_param_blr_rm, self.num_param_sa, self.num_param_con)
  
  def print_blrmodel_para_names(self):
    """
    print blr model parameter names
    """
    for i in range(self.num_param_blrmodel_rm):
      print("{:3d} {}".format(i, self.para_names['name'][i]))
  
  def locate_bhmass(self):
    """
    locate the index of black hole mass parameter
    """
    idx = -1
    for i in range(self.num_param_blrmodel_rm):
      if re.match("BLR_model_ln\(Mbh/1e6\)", self.para_names['name'][i]):
        idx = i 
        break
    
    if idx == -1:
      raise ValueError("No black hole mass parameter")
    else:
      return idx


class bplotlib(Param, Options, ParaName):
  """
  backend for ploting BRAINS results
  """
  def __init__(self, fname="./param"):
    Param.__init__(self, fname)
    Options.__init__(self, self.param['filedir'], self.param['flagdim'])
    ParaName.__init__(self, self.param['filedir'], self.param['flagdim'])
    
    self.VelUnit = np.sqrt( 6.672e-8 * 1.0e6 * 1.989e33 / (2.9979e10*8.64e4)) / 1.0e5
    self.C_Unit = 3.0e5/self.VelUnit
    self.con_scale = 1.0
    self.line_scale = 1.0

    self.file_dir = self.param['filedir'] +"/"
    self._load_data()
    self._load_results() 

    self.postfix=["con", "1d","2d", "sa", "sa1d", "sa2d", "sarm"]

    self.redshift = float(self.param['redshift'])

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

    # calculate line scale
    redshift = float(self.param["redshift"])
    linecenter = float(self.param["linecenter"])
    Vline = (line2d_profile[0, :, 0]/(1.0+redshift) - linecenter)/linecenter * (2.9979e5/self.VelUnit)
    dV = Vline[1]-Vline[0]
    flux = np.sum(line2d_profile[:, :, 1], axis=1) - 0.5*(line2d_profile[:, 0, 1] + line2d_profile[:, -1, 1])
    flux *= dV

    self.line_scale = nt/np.sum(flux)
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

      # calculate con_scale and line_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])
    
    elif self.param['flagdim'] == '1':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line data
      self.data['line_data'] = np.loadtxt(self.file_dir+self.param['linefile'])

      # calculate con_scale and line_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])
      self.line_scale = self.data['line_data'].shape[0]/np.sum(self.data['line_data'][:, 1])

    elif self.param['flagdim'] == '2':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line2d data
      self.data['line2d_data'] = self._load_line2d_data()

      # calculate con_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])
    
    elif self.param['flagdim'] == '3':
     # load profile data
     self.data["lp_data"] = np.loadtxt(self.file_dir+self.param['lineprofilefile'])

    elif self.param['flagdim'] == '4':
      self.data["sa_data"] = self._load_sa_data()

    elif self.param['flagdim'] == '5':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line data
      self.data['line_data'] = np.loadtxt(self.file_dir+self.param['linefile'])
      # load sa data
      self.data["sa_data"] = self._load_sa_data()

      # calculate con_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])

    elif self.param['flagdim'] == '6':
      # load continuum data
      self.data['con_data'] = np.loadtxt(self.file_dir+self.param['continuumfile'])
      # load line2d data
      self.data['line2d_data'] = self._load_line2d_data()
      # load sa data
      self.data["sa_data"] = self._load_sa_data()

      # calculate con_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])
    
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
  
  def _load_lp_rec(self):
    self.results['lp_rec'] = np.loadtxt(self.file_dir+"data/lineprofile_rec.txt")
    ns = int(self.sample_size)  
    nv = int(self.results['lp_rec'].shape[0]/ns)
    nc = self.results['lp_rec'].shape[1]
    self.results['lp_rec'] = np.reshape(self.results['lp_rec'], (ns, nv, nc), order='C')

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
      sample = np.loadtxt(self.file_dir + "data/posterior_sample_con.txt")
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
    
    elif self.param['flagdim'] == '1':
      sample = self.results['sample'] = np.loadtxt(self.file_dir + "data/posterior_sample_1d.txt")      
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
      self._load_line_rec()
      self._load_tran_rec()
      
    elif self.param['flagdim'] == '2':
      # load posterior samples
      sample = self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_2d.txt')
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      # load likelihoods
      self.results['sample_info'] = np.loadtxt(self.file_dir + 'data/posterior_sample_info_2d.txt')
      self.sample_size = self.results['sample'].shape[0]
      
      # load continuum reconstructions
      self._load_con_rec()
      self._load_line2d_rec()
      self._load_tran2d_rec()
    
    elif self.param['flagdim'] == '3':
      # load posterior samples
      sample = self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_lp.txt')
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      # load likelihoods
      self.results['sample_info'] = np.loadtxt(self.file_dir + 'data/posterior_sample_info_lp.txt')
      self.sample_size = self.results['sample'].shape[0]

      self._load_lp_rec()

    elif self.param['flagdim'] == '4':
      sample = self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa.txt')
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      self.sample_size = self.results['sample'].shape[0]
      self._load_sa_rec()

    
    elif self.param['flagdim'] == '5':
      sample = self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa1d.txt')
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      self.sample_size = self.results['sample'].shape[0]
      self._load_con_rec()
      self._load_line_rec()
      self._load_tran_rec()
      self._load_sa_rec()

    elif self.param['flagdim'] == '6':
      sample = self.results['sample'] = np.loadtxt(self.file_dir + 'data/posterior_sample_sa2d.txt')
      if sample.ndim == 1:
        self.results['sample'] = np.reshape(sample, (1, sample.shape[0]))
      self.sample_size = self.results['sample'].shape[0]
      self._load_line2d_rec()
      self._load_tran2d_rec()
      self._load_sa_rec()

    else:
      raise Exception("Incorrect FlagDim.")  
    
    print("results dict:", self.results.keys())
    return
  
  
  def get_con_data(self):
    if self.data.__contains__('con_data'):
      return self.data['con_data']
    else:
      raise ValueError("No con data!")
  
  def get_line_data(self):
    if self.data.__contains__('line_data'):
      return self.data['line_data']
    else:
      raise ValueError("No 1d line data!")
    
  def get_line2d_data(self):
    if self.data.__contains__('line2d_data'):
      return self.data['line2d_data']
    else:
      raise ValueError("No 2d line data!")
    
  def get_sa_data(self):
    if self.data.__contains__('sa_data'):
      return self.data['sa_data']
    else:
      raise ValueError("No sa data!")
  
  def get_lp_data(self):
    if self.data__contains__('lp_data'):
      return self.data['lp_data']
    else:
      raise ValueError("No lp data!")

  def plot_drw_parameters(self, doshow=False):
    """
    plot drw paramete distributions
    """
    
    print("# plotting drw parameter histograms.")

    if int(self.param['flagdim']) in [0, 1, 2, 4, 5, 6]:
      idx = self.num_param_blr_rm+self.num_param_sa
      if np.std(self.results['sample'][:, idx]) == 0.0:
        fig = corner.corner(self.results['sample'][:, idx+1:idx+3], labels=['ln sigma', 'ln tau'], smooth=True)
      else:
        fig = corner.corner(self.results['sample'][:, idx:idx+3], labels=['syserr', 'ln sigma', 'ln tau'], smooth=True)
      
      fig.savefig("drw.pdf", bbox_inches='tight')
      if doshow:
        plt.show()
      else:
        plt.close()
    else:
      print("The running does not have continnum reconstruction.")
  
  def get_narrow_line(self, grid_vel, imax, iepoch):
    """
    get narrow line component
    """
    
    linecenter = float(self.param["linecenter"])
    width_na = float(self.param["widthnarrowline"])
    width_na_err = float(self.param["widthnarrowlineerr"])
    shift_na = float(self.param["shiftnarrowline"])
    shift_na_err = float(self.param["shiftnarrowlineerr"])
    
    num_param_na = 0
    if int(self.param["flagnarrowline"]) > 0:
      num_param_na = 3

    # first treat the broadening
    if int(self.param["flaginstres"]) == 0:
      broad = float(self.param["instres"])
    elif int(self.param["flaginstres"]) == 1:
      broad = float(self.param["instres"]) + float(self.param["instreserr"]) \
            * self.results['sample'][imax, self.num_param_blrmodel_rm + num_param_na]
    elif int(self.param["flaginstres"]) == 2:
      broad = float(self.param["instres"]) + float(self.param["instreserr"]) \
            * self.results['sample'][imax, self.num_param_blrmodel_rm + num_param_na + iepoch]
    elif int(self.param["flaginstres"]) == 3:
      broad_all = np.loadtxt(os.path.join(self.param["filedir"], self.param["instresfile"]))
      broad = broad_all[iepoch, 0] + broad_all[iepoch, 1]*self.results['sample'][imax, self.num_param_blrmodel_rm + num_param_na + iepoch]

    if int(self.param["flagnarrowline"]) == 1:
      flux = float(self.param["fluxnarrowline"])
      width = width_na
      shift = shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      # note 1) flux scale. in the code, wavelength is converted into velocity, flux is not.
      # the narrow line flux needs to account for this points. 
      # note 2) the velocity unit. in the code, velocity is scaled by the unit.
      flux /= (linecenter/self.C_Unit)
      prof = flux/(np.sqrt(2.0*np.pi) * width/self.VelUnit) * np.exp( -0.5 *(grid_vel - shift)**2/width_br**2)
      return prof
    elif int(self.param["flagnarrowline"]) == 2:
      flux_na = float(self.param["fluxnarrowline"])
      flux_na_err = float(self.param["fluxnarrowlineerr"])
      flux = self.results['sample'][imax, self.num_param_blrmodel_rm] * flux_na_err + flux_na
      width = self.results['sample'][imax, self.num_param_blrmodel_rm+1] * width_na_err + width_na 
      shift = self.results['sample'][imax, self.num_param_blrmodel_rm+2] * shift_na_err + shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      # note 1) flux scale. in the code, wavelength is converted into velocity, flux is not.
      # the narrow line flux needs to account for this points. 
      # note 2) the velocity unit. in the code, velocity is scaled by the unit.
      flux /= (linecenter/self.C_Unit)
      prof = flux/(np.sqrt(2.0*np.pi) * width/self.VelUnit) * np.exp( -0.5 *(grid_vel - shift)**2/width_br**2)
      return prof
    elif int(self.param["flagnarrowline"]) == 3:
      flux = np.exp(self.results['sample'][imax, self.num_param_blrmodel_rm])
      width = self.results['sample'][imax, self.num_param_blrmodel_rm+1] * width_na_err + width_na 
      shift = self.results['sample'][imax, self.num_param_blrmodel_rm+2] * shift_na_err + shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      # note line scale
      flux /= self.line_scale
      # note 1) flux scale. in the code, wavelength is converted into velocity, flux is not.
      # the narrow line flux needs to account for this points. 
      # note 2) for flagnarrowline==3, the flux parameter is already include this point.
      # note 3) the velocity unit. in the code, velocity is scaled by the unit.
      prof = flux/(np.sqrt(2.0*np.pi) * width/self.VelUnit) * np.exp(-0.5 *(grid_vel - shift)**2/width_br**2)
      return prof

      
  def plot_results_con(self, doshow=False):
    """
    plot continuum results
    """

    if not int(self.param['flagdim']) in [0]:
      print("Flagdim =", self.param['flagdim'], "not continuum analysis.")
      return
    
    print("# plotting results con.")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)
    
    conlc = self.data["con_data"]

    con_mean_err = np.mean(conlc[:, 2])
    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')[0][0]
    sample = self.results['sample']
    syserr_con = (np.exp(np.mean(sample[:, idx_con])) - 1.0) * con_mean_err

    ax.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    ax.plot(con_date, con_mean, color='red', zorder=20) 
    ax.fill_between(con_date, y1=con_mean_low, y2=con_mean_upp, color='grey')
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Flux")
    ax.minorticks_on()
    ax.set_xlim(con_date[0], con_date[-1])

    fig.savefig("results_con.pdf", bbox_inches='tight')
    if doshow:
      plt.show()
    else:
      plt.close()
  
  def plot_results_1d(self, doshow=False, time_range=None, time_start=None):
    """
    plot 1d results
    """
    if not int(self.param['flagdim']) in [1, 5]:
      print("Flagdim =", self.param['flagdim'], "no 1D RM data.")
      return
    
    print("# plotting results 1d.")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    fig = plt.figure(figsize=(8, 4))
    fig.subplots_adjust(hspace=0.05)
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    
    if time_start is None:
      t0 = 0
    else:
      t0 = time_start

    conlc = self.data["con_data"]

    con_mean_err = np.mean(conlc[:, 2])
    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')[0][0]
    sample = self.results['sample']
    syserr_con = (np.exp(np.mean(sample[:, idx_con])) - 1.0) * con_mean_err

    ax1.errorbar(conlc[:, 0]-t0, conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    ax1.plot(con_date-t0, con_mean, color='red', zorder=20) 
    ax1.fill_between(con_date-t0, y1=con_mean_low, y2=con_mean_upp, color='grey')

    # plot different detrend
    idx_diff_trend_all =  np.nonzero(self.para_names['name'] == 'diff_trend')[0]
    num_param_diff_trend = len(idx_diff_trend_all)
    if num_param_diff_trend > 0:
      idx_diff_trend = idx_diff_trend_all[0]
      #note in BRAINS, light curves are converted into rest frame
      t1 = conlc[0, 0]/(1.0+self.redshift)
      t2 = conlc[-1, 0]/(1.0+self.redshift)
      tmed = (t1+t2)/2.0
      tspan = t2-t1

      trend = np.linspace(t1, t2, 200)
      ftrend = np.zeros((sample.shape[0], 200))

      a0 = np.zeros(sample.shape[0])
      for i in range(1, num_param_diff_trend+1):
        a0[:] += -sample[:, idx_diff_trend+i-1]/tspan * ((t2-tmed)**(i+1) - (t1-tmed)**(i+1))/(i+1)
      
      ftrend += a0[:, np.newaxis]
      for i in range(1, num_param_diff_trend+1):
        ftrend +=  sample[:, idx_diff_trend+i-1, np.newaxis] * (trend[np.newaxis, :] - tmed)**(i)
      
      fcon_mean = np.median(conlc[:, 1])
      ftrend_med, ftrend_low, ftrend_upp = np.quantile(ftrend, axis=0, q=(0.5, (1.0-0.683)/2, 1.0-(1.0-0.683)/2))/self.con_scale + fcon_mean
      ax1.plot(trend*(1.0+self.redshift) - t0, ftrend_med, ls='--', color='grey', lw=1)
      ax1.fill_between(trend*(1.0+self.redshift) - t0, y1=ftrend_low, y2=ftrend_upp, color='grey', alpha=0.5)

    #ax1.set_xlabel("Time")
    ax1.set_ylabel("Flux")
    ax1.set_xticklabels([])
    ax1.minorticks_on()
    fmax = np.max(conlc[:, 1])
    fmin = np.min(conlc[:, 1])
    ax1.set_ylim(fmin-0.2*(fmax-fmin), fmax+0.2*(fmax-fmin))
    ax1.set_xlim(con_date[0], con_date[-1])
    xlim = ax1.get_xlim()

    linelc = self.data["line_data"]
    line_mean_err = np.mean(linelc[:, 2])
    idx_line =  np.nonzero(self.para_names['name'] == 'sys_err_line')[0][0]
    syserr_line = (np.exp(np.mean(sample[:, idx_line])) - 1.0) * line_mean_err

    ax2.errorbar(linelc[:, 0]-t0, linelc[:, 1], yerr=np.sqrt(linelc[:, 2]*linelc[:, 2] + syserr_line*syserr_line), \
                marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    line_rec = self.results['line_rec']
    line_date = line_rec[0, :, 0]
    line_mean = np.quantile(line_rec[:, :, 1], axis=0, q=0.5)
    line_mean_upp = np.quantile(line_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    line_mean_low = np.quantile(line_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    ax2.plot(line_date-t0, line_mean, color='red', zorder=20) 
    ax2.fill_between(line_date-t0, y1=line_mean_low, y2=line_mean_upp, color='grey')
    
    if time_start is None:
      ax2.set_xlabel("Time")
    else:
      ax2.set_xlabel("Time - %d"%t0)
      
    ax2.set_ylabel("Flux")
    ax2.minorticks_on()
    ax2.set_xlim(xlim[0], xlim[1])
    
    fig.align_ylabels()

    if time_range is not None:
      ax1.set_xlim(time_range[0], time_range[1])
      ax2.set_xlim(time_range[0], time_range[1])

    fig.savefig("results_1d.pdf", bbox_inches='tight')
    if doshow:
      plt.show()
    else:
      plt.close()

  def plot_results_2d_style2018(self, doshow=False, time_range=None):
    """
    plot 2d results in the style of 2018 ApJ paper
    """

    if not int(self.param['flagdim']) in [2, 5, 6]:
      print("Flagdim =", self.param['flagdim'], "no 2D RM data.")
      return
    
    print("# plotting results 2d, with the 2018 style.")

    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    wave0 = float(self.param['linecenter'])

    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')[0][0]
    idx_line =  np.nonzero(self.para_names['name'] == 'sys_err_line')[0][0]
    
    date_line =  self.data['line2d_data']['time']
    prof = self.data['line2d_data']['profile'][:,:, 1]
    prof_err = self.data['line2d_data']['profile'][:, :, 2]

    nt, nv = prof.shape
    grid_wav = np.zeros(nv)
    grid_wav[:] = self.data['line2d_data']['profile'][0, :, 0]
    grid_vel = (grid_wav[:]/(1.0+float(self.param['redshift'])) - wave0)/wave0 *3.0e5
    dV = (grid_vel[1]-grid_vel[0]) / self.VelUnit
    
    conlc = self.data["con_data"]
    hblc=np.zeros((nt, 3))
    hblc[:, 1]=np.sum(prof, axis=1) * dV
    hblc[:, 2]=np.sqrt(np.sum(prof_err**2, axis=1)) * dV
    
    date0 = conlc[0, 0]
    conlc[:, 0]=conlc[:, 0]-date0
    date_line = date_line - date0 
    hblc[:, 0] = date_line
    
    con_mean_err = np.mean(conlc[:, 2])
    line_mean_err = np.mean(prof_err)
    sample = self.results['sample']
    syserr_con = (np.exp(np.mean(sample[:, idx_con])) - 1.0) * con_mean_err
    syserr_line = (np.exp(np.mean(sample[:, idx_line])) - 1.0) * line_mean_err
    # systematic error to line fluxes
    hblc_syserr = np.sqrt(syserr_line**2 * nv) * dV  

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    #========================================================================
    fig=plt.figure(1, figsize=(9, 8))
    cmap=plt.get_cmap('jet')
    
    
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    #========================================================================
    # subfig 5
    ax5=fig.add_axes([0.1, 0.2, 0.52, 0.15])

    prof_rec = np.zeros((nt, nv))
    prof_rec_max = np.zeros((nt, nv))
    line_rec = np.zeros(nt)
    line_rec_max = np.zeros(nt)
    
    chi = np.zeros(sample.shape[0])
    for i in range(self.results['line2d_rec'].shape[0]):
      prof_rec = self.results['line2d_rec'][i, :, :]
      line_rec = np.sum(prof_rec, axis=1) * dV
      #chi[i] = np.sum( (line_rec - hblc[:, 1]) * (line_rec - hblc[:, 1]))
      chi[i] = np.sum( ((prof_rec - prof)/prof_err)**2 )
    
    
    imax = np.argmin(chi)
    prof_rec_max = self.results['line2d_rec'][imax]
    line_rec_max = np.sum(prof_rec_max, axis=1) * dV

    plt.errorbar(hblc[:, 0], hblc[:, 1] * wave0*self.VelUnit/3e5, \
                 yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + hblc_syserr*hblc_syserr)* wave0*self.VelUnit/3e5, \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=0.7, markeredgewidth=0.5, zorder=32)
    #plt.plot(date_line, line_rec_max* wave0*self.VelUnit/3e5, color='red')
    
    # 68.3% quantile of reconstruction
    line_rec = np.sum(self.results['line2d_rec'], axis=2) * dV
    line_mean = np.quantile(line_rec, axis=0, q=0.5)
    line_mean_upp = np.quantile(line_rec, axis=0, q=(1.0-0.683)/2.0)
    line_mean_low = np.quantile(line_rec, axis=0, q=1.0 - (1.0-0.683)/2.0)
    ax5.plot(date_line, line_mean* wave0*self.VelUnit/3e5, color='red', zorder=20)
    ax5.fill_between(date_line, y1 = line_mean_upp* wave0*self.VelUnit/3e5, y2 = line_mean_low * wave0*self.VelUnit/3e5, color='grey')
    
    ax5.set_xlabel(r'$\rm Time\ (day)$')
    ax5.set_ylabel(r'$F_{\rm H\beta}$')
    
    ymax = np.max(hblc[:, 1]* wave0*self.VelUnit/3e5)
    ymin = np.min(hblc[:, 1]* wave0*self.VelUnit/3e5)
    ax5.set_ylim(ymin - 0.2*(ymax-ymin), ymax + 0.2*(ymax- ymin))
    
    xmax = hblc[-1, 0]
    xmin = hblc[0, 0]
    ax5.set_xlim(xmin - 0.2*(xmax-xmin), xmax + 0.2*(xmax- xmin))
    
    #========================================================================
     # subfig 4
    ax4=fig.add_axes([0.1, 0.35, 0.52, 0.15])

    # set a proper range
    xmin = min(conlc[0, 0], date_line[0])
    xmax = max(conlc[-1, 0], date_line[-1])
    dx = xmax-xmin
    xmin -= 0.05*dx 
    xmax += 0.05*dx 

    plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    idx = np.where((con_date-date0>=xmin)&(con_date-date0<=xmax))[0]
    ax4.plot(con_date[idx]-date0, con_mean[idx], color='red', zorder=20) 
    ax4.fill_between(con_date[idx]-date0, y1=con_mean_low[idx], y2=con_mean_upp[idx], color='grey')

    ax4.set_ylabel(r'$F_{\rm 5100}$')

    [i.set_visible(False) for i in ax4.get_xticklabels()]

    
    ax4.set_xlim(xmin, xmax)
    ax5.set_xlim(xmin, xmax)
    ax4.minorticks_on()
    ax5.minorticks_on()

    if time_range is not None:
      ax4.set_xlim(time_range[0], time_range[1])
      ax5.set_xlim(time_range[0], time_range[1])
    
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['xtick.top'] = False 
    plt.rcParams['ytick.right'] = False
    #========================================================================
    # subfig 1
    ax1 = fig.add_axes([0.1, 0.6, 0.25, 0.3])
    
    plt.imshow(prof, cmap=cmap, interpolation='gaussian', aspect='auto', origin='lower', \
               extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax = np.amax(prof), vmin=np.amin(prof))
    ax1.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax1.set_ylabel(r'$\rm Epoch~Number$')
    plt.text(0.08, 0.9, r'$\rm Data$', color='white', transform=ax1.transAxes)
    
    #========================================================================
    # subfig 2
    ax2=fig.add_axes([0.37, 0.6, 0.25, 0.3])
    
    plt.imshow(prof_rec_max, cmap=cmap, interpolation='gaussian',  aspect='auto', origin='lower', \
               extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax = np.amax(prof), vmin=np.amin(prof))
    ax2.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax2.text(0.08, 0.9, r'$\rm Model$', color='white', transform=ax2.transAxes)
    
    [i.set_visible(False) for i in ax2.get_yticklabels()]

    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = False 
    plt.rcParams['ytick.right'] = True
    #========================================================================
    # subfig 3, plot the best two epochs
    ax3=fig.add_axes([0.7, 0.6, 0.25, 0.3])
    chifit = np.zeros(nt)
    for i in range(nt):
      chifit[i] = np.sum( (prof_rec_max[i, :] - prof[i, :])**2/(prof_err[i, :]**2 + syserr_line**2) )
    chifit_sort = np.sort(chifit)

    offset = np.max(prof.flatten()) * 0.25
    
    idx = np.where(chifit == chifit_sort[0])
    i = idx[0][0]
    j = 0
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='grey', capsize=1, markeredgewidth=1, zorder=0)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[0], lw=2, zorder=1)

    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[0], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', label='Narrow Line', zorder=1)
    
    idx = np.where(chifit == chifit_sort[1])
    i = idx[0][0]
    j = 1
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='grey', capsize=1, markeredgewidth=1, zorder=0)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[1], lw=2, zorder=1)
    
    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[1], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', zorder=1)
    
    ax3.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax3.set_ylabel(r'$\rm Flux + offset$')
    ax3.set_xlim([grid_vel[0]/1.0e3, grid_vel[-1]/1.0e3])
    ax3.text(0.08, 0.9, r'$\rm Profile$', color='k', transform=ax3.transAxes)

    if int(self.param["flagnarrowline"]) > 0:
      ax3.legend(fontsize=8, handlelength=1.0, handletextpad=0.2, loc=(0.02, 0.8), frameon=False)
    

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    #========================================================================
    # subfig 6
    ax6=fig.add_axes([0.7, 0.2, 0.25, 0.3])
    
    idx = self.locate_bhmass()
    nhist, bhist = np.histogram(sample[:, idx]/np.log(10.0)+6.0, bins=25, density=True)
    nhist = gaussian_filter(nhist, 1.0)
    x0 = np.array(list(zip(bhist[:-1], bhist[1:]))).flatten()
    y0 = np.array(list(zip(nhist, nhist))).flatten()
    
    plt.stairs(nhist, bhist, fill=True)
    
    ax6.set_xlabel(r'$\log(M_\bullet/M_\odot)$')
    ax6.set_ylabel(r'$\rm Hist$')
    
    fig.savefig("results_2d.pdf", bbox_inches='tight')
    if doshow:
      plt.show()
    else:
      plt.close()
    
    # back to default
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
  
  def plot_results_2d_style2022(self, doshow=False, time_range=None):
    """
    plot 2d results in the style of 2022 ApJ paper
    """

    if not int(self.param['flagdim']) in [2, 5, 6]:
      print("Flagdim =", self.param['flagdim'], "no 2D RM data.")
      return
    
    print("# plotting results 2d, with the 2022 style.")

    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    wave0 = float(self.param['linecenter'])

    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')[0][0]
    idx_line =  np.nonzero(self.para_names['name'] == 'sys_err_line')[0][0]
    
    date_line =  self.data['line2d_data']['time']
    prof = self.data['line2d_data']['profile'][:,:, 1]
    prof_err = self.data['line2d_data']['profile'][:, :, 2]

    nt, nv = prof.shape
    grid_wav = np.zeros(nv)
    grid_wav[:] = self.data['line2d_data']['profile'][0, :, 0]
    grid_vel = (grid_wav[:]/(1.0+float(self.param['redshift'])) - wave0)/wave0 *3.0e5
    dV = (grid_vel[1]-grid_vel[0]) / self.VelUnit
    
    conlc = self.data["con_data"]
    hblc=np.zeros((nt, 3))
    hblc[:, 1]=np.sum(prof, axis=1) * dV
    hblc[:, 2]=np.sqrt(np.sum(prof_err**2, axis=1)) * dV
    
    date0 = conlc[0, 0]
    conlc[:, 0]=conlc[:, 0]-date0
    date_line = date_line - date0 
    hblc[:, 0] = date_line
    
    con_mean_err = np.mean(conlc[:, 2])
    line_mean_err = np.mean(prof_err)
    sample = self.results['sample']
    syserr_con = (np.exp(np.mean(sample[:, idx_con])) - 1.0) * con_mean_err
    syserr_line = (np.exp(np.mean(sample[:, idx_line])) - 1.0) * line_mean_err
    # systematic error to line fluxes
    hblc_syserr = np.sqrt(syserr_line**2 * nv) * dV  

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    #========================================================================
    fig=plt.figure(1, figsize=(9, 8))
    cmap=plt.get_cmap('jet')
    
    
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    #========================================================================
    # subfig 6
    ax6=fig.add_axes([0.37, 0.2, 0.55, 0.15])

    prof_rec = np.zeros((nt, nv))
    prof_rec_max = np.zeros((nt, nv))
    line_rec = np.zeros(nt)
    line_rec_max = np.zeros(nt)
    
    chi = np.zeros(sample.shape[0])
    for i in range(self.results['line2d_rec'].shape[0]):
      prof_rec = self.results['line2d_rec'][i, :, :]
      line_rec = np.sum(prof_rec, axis=1) * dV
      #chi[i] = np.sum( (line_rec - hblc[:, 1]) * (line_rec - hblc[:, 1]))
      chi[i] = np.sum( ((prof_rec - prof)/prof_err)**2 )
    
    imax = np.argmin(chi)
    prof_rec_max = self.results['line2d_rec'][imax]
    line_rec_max = np.sum(prof_rec_max, axis=1) * dV

    plt.errorbar(hblc[:, 0], hblc[:, 1] * wave0*self.VelUnit/3e5, \
                 yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + hblc_syserr**2)* wave0*self.VelUnit/3e5, \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=0.7, markeredgewidth=0.5, zorder=32)
    #plt.plot(date_line, line_rec_max* wave0*self.VelUnit/3e5, color='red')
    
    # 68.3% quantile of reconstruction
    line_rec = np.sum(self.results['line2d_rec'], axis=2) * dV
    line_mean = np.quantile(line_rec, axis=0, q=0.5)
    line_mean_upp = np.quantile(line_rec, axis=0, q=(1.0-0.683)/2.0)
    line_mean_low = np.quantile(line_rec, axis=0, q=1.0 - (1.0-0.683)/2.0)
    ax6.plot(date_line, line_mean* wave0*self.VelUnit/3e5, color='red', zorder=20)
    ax6.fill_between(date_line, y1 = line_mean_upp* wave0*self.VelUnit/3e5, y2 = line_mean_low * wave0*self.VelUnit/3e5, color='grey')
    
    ax6.set_xlabel(r'$\rm Time\ (day)$')
    ax6.set_ylabel(r'$F_{\rm H\beta}$')
    
    ymax = np.max(hblc[:, 1]* wave0*self.VelUnit/3e5)
    ymin = np.min(hblc[:, 1]* wave0*self.VelUnit/3e5)
    ax6.set_ylim(ymin - 0.2*(ymax-ymin), ymax + 0.2*(ymax- ymin))
    
    xmax = hblc[-1, 0]
    xmin = hblc[0, 0]
    ax6.set_xlim(xmin - 0.2*(xmax-xmin), xmax + 0.2*(xmax- xmin))
    ax6.tick_params(right=True, labelright=True, left=True, labelleft=False)
    ax6.yaxis.set_label_position('right')
    #========================================================================
     # subfig 5
    ax5=fig.add_axes([0.37, 0.35, 0.55, 0.15])

    # set a proper range
    xmin = min(conlc[0, 0], date_line[0])
    xmax = max(conlc[-1, 0], date_line[-1])
    dx = xmax-xmin
    xmin -= 0.05*dx 
    xmax += 0.05*dx 

    plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    idx = np.where((con_date-date0>=xmin)&(con_date-date0<=xmax))[0]
    ax5.plot(con_date[idx]-date0, con_mean[idx], color='red', zorder=20) 
    ax5.fill_between(con_date[idx]-date0, y1=con_mean_low[idx], y2=con_mean_upp[idx], color='grey')

    ax5.set_ylabel(r'$F_{\rm 5100}$')

    [i.set_visible(False) for i in ax5.get_xticklabels()]

    
    ax5.set_xlim(xmin, xmax)
    ax5.tick_params(labelright=True, labelleft=False)
    ax5.yaxis.set_label_position('right')
    
    xlim = ax5.get_xlim()
    ax6.set_xlim(xlim[0], xlim[1])
    ax5.minorticks_on()
    ax6.minorticks_on()

    if time_range is not None:
      ax5.set_xlim(time_range[0], time_range[1])
      ax6.set_xlim(time_range[0], time_range[1])

    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['xtick.top'] = False
    plt.rcParams['ytick.right'] = False
    #========================================================================
    # subfig 1
    ax1 = fig.add_axes([0.1, 0.6, 0.25, 0.3])
    
    plt.imshow(prof, cmap=cmap, interpolation='gaussian', aspect='auto', origin='lower',\
               extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax = np.amax(prof), vmin=np.amin(prof))
    ax1.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax1.set_ylabel(r'$\rm Epoch~Number$')
    plt.text(0.08, 0.9, r'\bf Data', color='white', transform=ax1.transAxes)
    
    #========================================================================
    # subfig 2
    ax2=fig.add_axes([0.37, 0.6, 0.25, 0.3])
    
    plt.imshow(prof_rec_max, cmap=cmap, interpolation='gaussian',  aspect='auto', origin='lower', \
               extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax = np.amax(prof), vmin=np.amin(prof))
    ax2.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax2.text(0.08, 0.9, r'\bf Model', color='white', transform=ax2.transAxes)
    
    [i.set_visible(False) for i in ax2.get_yticklabels()]

    #========================================================================
    # subfig 3
    prof_diff = prof - prof_rec_max
    ax3=fig.add_axes([0.64, 0.6, 0.31, 0.3])
    img = prof_diff/np.sqrt(prof_err**2 + syserr_line**2)
    vmax = np.max((abs(np.max(img)), abs(np.min(img)), 5.5))
    vmin = -vmax
    cmap=ax3.imshow(img,  aspect='auto', origin='lower', cmap='jet',
                    extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax = vmax, vmin = vmin)
    
    plt.colorbar(cmap, ticks=[-5, -2.5, 0.0, 2.5, 5.0])
    ax3.text(0.08, 0.9, r'\bf Residuals', color='white', transform=ax3.transAxes)
    ax3.set_xlabel(r'$\rm Velocity\ (10^3\ km\ s^{-1})$')
    [i.set_visible(False) for i in ax3.get_yticklabels()]

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
    #========================================================================
    # subfig 4, plot the best two epochs
    ax4=fig.add_axes([0.1, 0.2, 0.25, 0.3])
    chifit = np.zeros(nt)
    for i in range(nt):
      chifit[i] = np.sum( (prof_rec_max[i, :] - prof[i, :])**2 /(prof_err[i, :]**2 + syserr_line**2) )
    chifit_sort = np.sort(chifit)

    offset = np.max(prof.flatten()) * 0.2
    
    idx = np.where(chifit == chifit_sort[0])
    i = idx[0][0]
    j = 0
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='grey', capsize=1, markeredgewidth=1)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[0], lw=2)
    
    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[0], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', zorder=1, label='Narrow Line')

    idx = np.where(chifit == chifit_sort[1])
    i = idx[0][0]
    j = 1
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='grey', capsize=1, markeredgewidth=1)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[1], lw=2)
    
    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[1], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', zorder=1)

    ax4.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
    ax4.set_ylabel(r'$\rm Flux + offset$')
    ax4.set_xlim([grid_vel[0]/1.0e3, grid_vel[-1]/1.0e3])
    ax4.text(0.08, 0.9, r'$\rm Profile$', color='k', transform=ax4.transAxes)
    
    if int(self.param["flagnarrowline"]) > 0:
      ax4.legend(fontsize=8, handlelength=1.0, handletextpad=0.2, loc=(0.02, 0.8), frameon=False)
    
    fig.savefig("results_2d.pdf", bbox_inches='tight')
    if doshow:
      plt.show()
    else:
      plt.close()
    
    # back to default
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True
  
  def plot_results_lp(self, doshow=False):
    """
    plot lp results
    """
    if not int(self.param['flagdim']) in [3]:
      print("Flagdim =", self.param['flagdim'], "no lineprofile data.")
      return
    
    print("# plotting results lineprofile.")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    profile = self.data['lp_data']
    profile_rec = self.results['lp_rec']
   
    wave = profile_rec[0, :, 0]
    profile_mean     = np.quantile(profile_rec[:, :, 1], axis=0, q=0.5)
    profile_mean_upp = np.quantile(profile_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    profile_mean_low = np.quantile(profile_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)

    # get the mean and 

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    ax.errorbar(profile[:, 0], profile[:, 1], yerr=profile[:, 2], ls='none')
    ax.plot(wave, profile_mean)
    ax.fill_between(wave, y1=profile_mean_low, y2=profile_mean_upp, color='grey', zorder=0)
    
    ax.minorticks_on()
    ax.set_xlabel("Wavelength")
    ax.set_ylabel("Flux")
    
    fig.savefig("results_lp.pdf", bbox_inches='tight')
    
    if doshow:
      plt.show()
    else:
      plt.close()

  
  def plot_results_sa(self, doshow=False, show_offset=False, subtract_offset=False, phase_limit=None, column_first=True,
                      average_baseline=None):
    if not int(self.param['flagdim']) in [4, 5, 6]:
      print("Flagdim =", self.param['flagdim'], "no  SA data.")
      return

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    #============================================
    # figures 
    fig_avg = plt.figure(figsize=(6, 6))
    ax_avg = fig_avg.add_axes((0.1, 0.1, 0.8, 0.55))
    ax2_avg = fig_avg.add_subplot((0.1, 0.65, 0.8, 0.25))
    #============================================
    
    ns = self.data["sa_data"]["profile"].shape[0] # number of spectra
    nv = self.data["sa_data"]["profile"].shape[1] # number of wavelength
    nb = self.data["sa_data"]["phase"].shape[0]   # number of baseline
    nb_per_ns = int(nb/ns)
    print(ns, nv, nb)

    profile_rec = self.results["sa_profile_rec"]
    phase_rec = self.results["sa_phase_rec"]
    baseline = self.data["sa_data"]["baseline"]
    wave = profile_rec[0, :, 0] # reconstruction and data have the same wavelength bin

    # calculate continuum offset
    if subtract_offset == True or show_offset == True:
      Factor = 1.0/360.0 * 3.08567758e+18/(3e10*24*3600.0)
      phase_offset_rec = np.zeros((phase_rec.shape[0], nb, nv))
      phase_offset_rec_mean = np.zeros((nb, nv))
      idx_sa_extpar = np.nonzero(self.para_names['name'] == 'SA_Ext_Par_DA')[0][0]
      sample = self.results['sample']      
      
      # offset = fline/(1+fline) * xc/DA
      # phase_offset = -360 * offset * B/lambda [deg]
      # note the negative sign. In the code, phase_obs = phase_blr - phase_offset
      alphac = profile_rec[:, :, 1]/(1.0+profile_rec[:, :, 1])*sample[:, idx_sa_extpar+4, np.newaxis]/np.exp(sample[:, idx_sa_extpar,  np.newaxis])  # ld/Mpc
      betac  = profile_rec[:, :, 1]/(1.0+profile_rec[:, :, 1])*sample[:, idx_sa_extpar+5, np.newaxis]/np.exp(sample[:, idx_sa_extpar,  np.newaxis])  # ld/Mpc
      for loop_b in range(nb):
        phase_offset_rec[:, loop_b, :] = -(alphac*baseline[loop_b, 0] + betac*baseline[loop_b, 1])/wave/Factor
        phase_offset_rec_mean[loop_b, :] = np.median(phase_offset_rec[:, loop_b, :], axis=0)

    #========================================================================
    ncol = ns
    nrow = int(nb/ns)+1
    if nb > 6:
      ncol = int(np.ceil(nb/6))
      nrow = 6+1

    fig, axs = plt.subplots(nrow, ncol)
    fig.set_size_inches(4*ncol, 8)
    axs = np.reshape(axs, (nrow, ncol))
    
    fig.subplots_adjust(wspace=0, hspace=0)

    # line profile
    ax_spec = axs[0, :]
    for loop_e in range(ns, ncol):
      ax = ax_spec[loop_e]
      ax.set_axis_off()
    
    flux_rec_mean = np.mean(profile_rec, axis=0)
    flux_rec_upp = np.quantile(profile_rec, axis=0, q=1.0-(1.0-0.683)/2.0)
    flux_rec_low = np.quantile(profile_rec, axis=0, q=(1.0-0.683)/2.0)
    for loop_e in range(ns):
      ax = ax_spec[loop_e]
      ax.plot(flux_rec_mean[:, 0], flux_rec_mean[:, 1], color="C0", lw=2)
      ax.fill_between(flux_rec_mean[:, 0], y1=flux_rec_upp[:, 1], y2=flux_rec_low[:, 1], color='C0', alpha=0.5)

      ax2_avg.plot(flux_rec_mean[:, 0], flux_rec_mean[:, 1], color="C1", lw=2)
      ax2_avg.fill_between(flux_rec_mean[:, 0], y1=flux_rec_upp[:, 1], y2=flux_rec_low[:, 1], color='C1', alpha=0.5)


    profile = self.data["sa_data"]["profile"]
    for loop_e in range(ns):
      ax = ax_spec[loop_e]
      ax.errorbar(profile[loop_e, :, 0], profile[loop_e, :, 1], yerr=profile[loop_e, :, 2], ls='none', color="C0", 
                  marker='o', markersize=2, zorder=0, elinewidth=1)
      
      ax2_avg.errorbar(profile[loop_e, :, 0], profile[loop_e, :, 1], yerr=profile[loop_e, :, 2], ls='none', color="C0", 
                       marker='o', zorder=0, elinewidth=1)

      if loop_e > 0:
        ax.set_yticklabels([])
      else:
        ax.set_ylabel("Flux")
      
      ax.minorticks_on()
    
    # then phase 
    ax_vphi = axs[1:, :] 
    phase_rec_mean = np.mean(phase_rec, axis=0)
    phase_rec_low = np.quantile(phase_rec, axis=0, q=(1.0-0.683)/2)
    phase_rec_upp = np.quantile(phase_rec, axis=0, q=1.0-(1.0-0.683)/2)

    if column_first == True:
      print("Column first changes in baselines!")
    else:
      print("Row first changes in baselines!")

    for loop_b in range(nb):
      if column_first == True:  # column first changes
        loop_e = loop_b%(ncol)
        loop_b_p = loop_b//(ncol)
      else:  # row first changes 
        loop_b_p = loop_b%(nrow-1)
        loop_e = loop_b//(nrow-1)
      ax = ax_vphi[loop_b_p, loop_e]
      phi = phase_rec_mean[loop_b, :, :]
      phi_upp = phase_rec_upp[loop_b, :, 1]  
      phi_low = phase_rec_low[loop_b, :, 1]                
      ax.plot(phi[:, 0], phi[:, 1], color="C{0}".format(loop_b_p), lw=2)
      ax.fill_between(phi[:, 0], y1=phi_upp, y2=phi_low, color="C{0}".format(loop_b_p), alpha=0.5)

      if show_offset == True:
        #note: continuum offset already included in reconstruction
        #      in the code,  phase_obs =  phase_blr - phase_offset
        phi_offset = phase_offset_rec_mean[loop_b, :]
        ax.plot(phi[:, 0], -phi_offset, ls='dotted', color="C{0}".format(loop_b_p), zorder=5)
        ax.plot(phi[:, 0], phi[:, 1]+phi_offset, ls='--', color="C{0}".format(loop_b_p), zorder=5)

    phase = self.data["sa_data"]["phase"]
    ymin=np.finfo(np.float64).max
    ymax=np.finfo(np.float64).min
    for loop_b in range(nb):
      if column_first == True:
        loop_e = loop_b%(ncol)
        loop_b_p = loop_b//(ncol)
      else:
        loop_b_p = loop_b%(nrow-1)
        loop_e = loop_b//(nrow-1)
      ax = ax_vphi[loop_b_p, loop_e]
      x = phase[loop_b, :, 0]
      y = phase[loop_b, :, 1]
      e = phase[loop_b, :, 2]
      ecopy = copy.copy(e)
      
      # cope with masked points 
      idx_mask = np.where(e<=0.0)[0]
      if len(idx_mask) > 0:
        ecopy[idx_mask] = 0.0       
        # plot masked points with black color
        ax.errorbar(x[idx_mask], y[idx_mask], yerr=ecopy(idx_mask), ls="none", marker="o", markersize=2, color="k", elinewidth=1, zorder=10)
      
      ax.errorbar(x, y, yerr=ecopy, ls="none", marker="o", markersize=2, color="C{0}".format(loop_b_p), elinewidth=1)

      if loop_b_p < nrow-2:
        ax.set_xticklabels([])
      else:
        ax.set_xlabel(r"Wavelength")

      if loop_e > 0:
        ax.set_yticklabels([])
      else:
        ax.set_ylabel(r"$\phi$")
      
      ax.minorticks_on()
      ylim = ax.get_ylim()
      ymin = np.min((ylim[0], ymin))
      ymax = np.max((ylim[1], ymax))
    
    for loop_b in range(nb):
      if column_first == True:
        loop_e = loop_b%(ncol)
        loop_b_p = loop_b//(ncol)
      else:
        loop_b_p = loop_b%(nrow-1)
        loop_e = loop_b//(nrow-1)
        
      ax = ax_vphi[loop_b_p, loop_e]

      if phase_limit is not None:
        ax.set_ylim(phase_limit[0], phase_limit[1])
      else:
        ax.set_ylim(ymin, ymax)

    fig.align_labels()

    fig.savefig("results_sa.pdf", bbox_inches='tight')
    
    # averaged phase
    # do averaging on selected baselines
    if average_baseline is not None:
      if isinstance(average_baseline, int):
        phase_rec_mean_max = np.max(np.abs(phase_rec_mean[:, :, 1]+phase_offset_rec_mean), axis=1)
        # note in ascending order 
        tmp = np.sort(phase_rec_mean_max)
        if average_baseline <= nb:
          tmp_sel = tmp[nb-average_baseline]
        else:
          tmp_sel = tmp[0]
        idx_baseline= np.where(phase_rec_mean_max>=tmp_sel)[0]

      elif isinstance(average_baseline, list):
        idx_baseline = average_baseline
        
    else:
      idx_baseline = np.array(np.arange(nb))
    
    print("averaging baselines:", idx_baseline)

    phase_avg = np.zeros((nv, 3))
    norm = np.zeros(nv)
    for loop_b in idx_baseline:
      y = phase[loop_b, :, 1]
      e = phase[loop_b, :, 2]

      if subtract_offset == True:
        phi_offset = phase_offset_rec_mean[loop_b, :]
        y = y + phi_offset 

      weight = 1.0/e**2
      norm[:] += weight
      phase_avg[:, 1] += y*weight
      phase_avg[:, 2] += e**2*weight**2

    phase_avg[:, 0] = phase[0, :, 0]
    phase_avg[:, 1] /= norm 
    phase_avg[:, 2] = np.sqrt(phase_avg[:, 2])/norm
    
    ax_avg.errorbar(phase_avg[:, 0], phase_avg[:, 1], yerr=phase_avg[:, 2], ls='--', marker='o', color="C0")

    # averaged phase reconstruction 
    if subtract_offset == False:
      phase_rec_avg = np.mean(phase_rec_mean[idx_baseline, :, 1], axis=0)
      phase_rec_upp_avg = np.mean(phase_rec_upp[idx_baseline, :, 1], axis=0)
      phase_rec_low_avg = np.mean(phase_rec_low[idx_baseline, :, 1], axis=0)
    else:
      phase_rec_avg = np.mean(phase_rec_mean[idx_baseline, :, 1]    + phase_offset_rec_mean[idx_baseline, :], axis=0)
      phase_rec_upp_avg = np.mean(phase_rec_upp[idx_baseline, :, 1] + phase_offset_rec_mean[idx_baseline, :], axis=0)
      phase_rec_low_avg = np.mean(phase_rec_low[idx_baseline, :, 1] + phase_offset_rec_mean[idx_baseline, :], axis=0)

    ax_avg.plot(wave, phase_rec_avg, color='C1')
    ax_avg.fill_between(wave, y1=phase_rec_upp_avg, y2=phase_rec_low_avg, alpha=0.5, color="C1")
    ax_avg.axhline(y=0.0, ls='--', color='grey', lw=1)

    ax_avg.set_xlabel(r"Wavelength")
    ax_avg.set_ylabel(r"$\phi$")
    ax_avg.minorticks_on()
    ax2_avg.set_xticklabels([])
    ax2_avg.set_ylabel("Flux")

    fig_avg.savefig("results_sa_avg.pdf", bbox_inches='tight')

    if doshow==True:
      plt.show()
    else:
      plt.close()

  
  def find_max_prob(self):
    """
    find the parameter set with the maximum prob
    """
    idx = np.argmax(self.results['sample_info'])
    return idx, self.results['sample'][idx, :]
  
  def check_cloud_dist_p14(self):
    """
    check Appendix A in Pancoast et al. 2014a
    """
    Ec = -0.5
    Lc = 1.0
    lam = 0.1
    n = 10000
    chi = np.random.randn(n) * lam
    E =(chi+ 1.0) * Ec
    #plt.hist(E)
    Lmax = np.sqrt(2.0 * (E + 1.0))
    L = Lmax * lam * np.log( (np.exp(1.0/lam) - 1.0) * np.random.random(n)  + 1.0)
    
    vr2 = 2.0*(E + 1.0) - L*L
    u = np.random.rand(n)
    sig = [1 if ut > 0.5 else -1 for ut in u]
    
    idx = np.where(vr2 < 0)
    vr2[idx[0]] = 0
    vr = np.sqrt(vr2) * sig
    
    u = np.random.rand(n)
    sig = [1 if ut > 0.5 else -1 for ut in u]
    vphi = L * sig
    
    
    plt.plot(vr, vphi, ls='none', marker='.', markersize=1)
    
    r= 1.0
    theta = np.linspace(0, 2.0*np.pi, 100)
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    plt.plot(x, y)
    plt.show()
  
  def plot_limits(self):
    """
    plot CDNest limits of parameters with levels
    """
    if int(self.param['flagdim']) < 0:
      print("Flagdim =", self.param["flagdim"], "no limits")
      return 
    
    print("# plotting limits.")

    str_dim = self.postfix[int(self.param['flagdim'])] 

    pdf = PdfPages("limits_"+str_dim+".pdf")
  
    limits = np.loadtxt("../data/limits_"+str_dim+".txt", comments='#')
  
    for i in range(self.num_param_blrmodel_rm):
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.plot(limits[:, 0], limits[:, (i+1)*2-1], marker='o')
      ax.plot(limits[:, 0], limits[:, (i+1)*2], marker='o')
      ax.set_xlabel("Level")
      ax.set_ylabel("Limits " + self.para_names['name'][i])
      pdf.savefig()
      plt.close()
  
    pdf.close()
  
  def plot_tran2d(self, tau_range=None, vel_range=None, doshow=False):
    """
    plot transfer function with the maximum prob
    """
    print("# plotting 2d transfer function.")

    if int(self.param['flagdim']) in [2, 5, 6]:
      idx, par = self.find_max_prob()
      
      extent = [self.data['line2d_data']['profile'][0,0,0],self.data['line2d_data']['profile'][0,-1,0], \
                        self.results['tau_rec'][idx][0], self.results['tau_rec'][idx][-1]]

      plt.rcParams['xtick.direction'] = 'out'
      plt.rcParams['ytick.direction'] = 'out'
      fig = plt.figure()
      ax = fig.add_subplot(111)
      ax.imshow(self.results['tran2d_rec'][idx], aspect='auto', origin='lower', \
                extent=extent,\
                interpolation='gaussian')
      
      if tau_range is not None:
        ax.set_ylim(tau_range[0], tau_range[1])
      if vel_range is not None:
        ax.set_xlim(vel_range[0], vel_range[1])

      ax.set_xlabel(r'Wavelength')
      ax.set_ylabel(r'Time Lag')
      ax.set_title("Transfer function with the maximum prob")
      plt.rcParams['xtick.direction'] = 'in'
      plt.rcParams['ytick.direction'] = 'in'
      fig.savefig("tran2d.pdf", bbox_inches='tight')
      if doshow:
        plt.show()
      else:
        plt.close()
    else:
      print("Flagdim =", self.param["flagdim"], "no trans2d")
      return
  
  def plot_tran1d(self, tau_range=None, doshow=False):
    """
    plot transfer function 1d
    """

    print("# plotting 1d transfer function.")

    if int(self.param['flagdim']) in [1, 5]:
      tran_rec = self.results['tran_rec']
      fig = plt.figure()
      ax = fig.add_subplot(111)
      for i in range(0, tran_rec.shape[0], int(tran_rec.shape[0]/100+1.0)):
        ax.plot(tran_rec[i, :, 0], tran_rec[i, :, 1], lw=0.5)
      
      if tau_range is not None:
        ax.set_xlim(tau_range[0], tau_range[1])
      
      fig.savefig("tran1d.pdf", bbox_inches='tight')
      if doshow:
        plt.show()
      else:
        plt.close()

    elif int(self.param['flagdim']) in [2, 6, 7]:
      tau_rec = self.results['tau_rec']
      tran2d_rec = self.results['tran2d_rec']
      tran1d_rec = np.sum(tran2d_rec, axis=2)
      
      fig = plt.figure()
      ax = fig.add_subplot(111)
      for i in range(0, tran1d_rec.shape[0], int(tran1d_rec.shape[0]/100+1.0)):
        ax.plot(tau_rec[i, :], tran1d_rec[i, :], lw=0.5)
      
      if tau_range is not None:
        ax.set_xlim(tau_range[0], tau_range[1])

      ax.set_xlabel(r"Time Lag")
      ax.set_ylabel(r"Transfer Function")
      fig.savefig("tran1d.pdf", bbox_inches='tight')
      if doshow:
        plt.show()
      else:
        plt.close()
    else:
      print("Flagdim =", self.param["flagdim"], "no trans1d")
      return
  
  def plot_blrmodel_para_hist(self, para_indx=None, doshow=False):
    """
    plot histograms of BLR model parameters
    """
    print("# plotting blrmodel parameter histograms.")

    prange = np.zeros((self.results['sample'].shape[1], 2))
    prange[:, 0] = np.min(self.results['sample'], axis=0)
    prange[:, 1] = np.max(self.results['sample'], axis=0)
    dprange = prange[:, 1]-prange[:, 0]
    prange[:, 0] -= 0.1*dprange 
    prange[:, 1] += 0.1*dprange
    
    print("BLR para num:", self.num_param_blrmodel_rm)
    if int(self.param['flagdim']) >= 1:
      if para_indx is None:
        label = []
        for i in range(self.num_param_blrmodel_rm):
          label += [self.para_names['name'][i][10:]]

        fig = corner.corner(self.results['sample'][:, :self.num_param_blrmodel_rm], \
                            range = prange[:self.num_param_blrmodel_rm, :], \
                            smooth=True, smooth1d=True, labels=label)
      else:
        label = []
        for i in para_indx:
          label += [self.para_names['name'][i][10:]]

        fig = corner.corner(self.results['sample'][:, para_indx], range = prange[para_indx, :], \
                            smooth=True, smooth1d=True, labels=label)
      
      if doshow:
        plt.show()
      else:
        plt.close()
  
      fig.savefig("hist.pdf", bbox_inches='tight')
    else:
      print("Flagdim =", self.param["flagdim"], "no BLR parameters")
      return
  
  def plot_clouds(self, fname=None, cmap='bwr', range=[-10, 10], objname=None, velocity=True, format="jpg", doshow=False):
    """
    plot clouds' distribution
    """
    
    if fname is None:
      raise ValueError("please input a file that restores the data for clouds' distribution.")
    
    if int(self.param["flagdim"]) == 1:
      velocity = False 
    
    print("# plotting clouds.")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(4, 8))
    ax1 = fig.add_axes((0.2, 0.405, 0.6, 0.3))
    ax2 = fig.add_axes((0.2, 0.1, 0.6, 0.3))

    if velocity == True:
      ax3 = fig.add_axes((0.2, 0.715, 0.6, 0.02))
    
    clouds  = np.loadtxt(fname)
    if clouds.shape[0] > 20000:
      clouds_rm = clouds[::10, :] 
    else:
      clouds_rm = clouds

    if velocity==True:
      vlos_max = np.max(np.abs(clouds_rm[:, 3]))
      
      norm = colors.Normalize(vmin=-vlos_max, vmax=vlos_max)
      clouds_rm[:, 3] = -clouds_rm[:, 3]  # positive is recending
    
    if velocity == True:
      cax=ax1.scatter(clouds_rm[:, 0], clouds_rm[:, 1], c=clouds_rm[:, 3], cmap=cmap, norm=norm, alpha=0.7, s=10, edgecolors='k', linewidths=0.5)
      ax2.scatter(clouds_rm[:, 0], clouds_rm[:, 2], c=clouds_rm[:, 3], cmap=cmap, norm=norm, alpha=0.7, s=10, edgecolors='k', linewidths=0.5)
    else:
      ax1.scatter(clouds_rm[:, 0], clouds_rm[:, 1], alpha=0.7, s=10, linewidths=0.5)
      ax2.scatter(clouds_rm[:, 0], clouds_rm[:, 2], alpha=0.7, s=10, linewidths=0.5)
    
    # check the range 
    clouds_max = np.max(clouds_rm[:, 0:3])
    if np.max(np.abs(range)) < clouds_max/10.0:
      print("cloud locations significant exceed input/default range=[%.2f, %.2f], better to set a larger range."%(range[0], range[1]))

    ax1.set_xlim(range[0], range[1])
    ax1.set_ylim(range[0], range[1])
    ax2.set_xlim(range[0], range[1])
    ax2.set_ylim(range[0], range[1])
    ax1.minorticks_on()
    ax2.minorticks_on()

    xt = ax1.get_xticks()
    ax1.set_yticks(xt)
    xt = ax2.get_xticks()
    ax2.set_yticks(xt)

    if objname is not None:
      ax1.text(0.06, 0.95, r"\bf "+ objname, transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top', bbox=dict(boxstyle="round", facecolor='w', alpha=0.7))
    
    [xt.set_visible(False) for xt in ax1.get_xticklabels()]
    
    ax1.set_ylabel(r"$y$ (lt-day)")
    ax2.set_xlabel(r"$x$ (lt-day)")
    ax2.set_ylabel(r"$z$ (lt-day)")
    
    if velocity == True:
      cbar = plt.colorbar(cax, cax=ax3, location='top', label=r"LOS Velocity (km~s$^{-1}$)")

    if format=="jpg":    
      fig.savefig("fig_clouds.jpg", bbox_inches='tight', dpi=500)
    else:
      fig.savefig("fig_clouds.pdf", bbox_inches='tight')
    
    if doshow:
      plt.show()
    else:
      plt.close()
  
  def plot_clouds_los(self, fname=None, cmap='bwr', range=[-10, 10], objname=None, velocity=True, format="jpg", doshow=False):
    """
    plot clouds' distribution viewed from line of sight
    """
    if fname is None:
      raise ValueError("please input a file that restores the data for clouds' distribution.")
    
    if int(self.param["flagdim"]) == 1:
      velocity = False
    
    print("# plotting clouds los.")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    fig = plt.figure(figsize=(4, 4))
    ax1 = fig.add_axes((0.2, 0.2, 0.6, 0.6))

    if velocity == True:
      ax2 = fig.add_axes((0.2, 0.815, 0.6, 0.02))
    
    clouds  = np.loadtxt(fname)
    if clouds.shape[0] > 20000:
      clouds_rm = clouds[::10, :] 
    else:
      clouds_rm = clouds

    if velocity == True:
      vlos_max = np.max(np.abs(clouds_rm[:, 3]))
      
      norm = colors.Normalize(vmin=-vlos_max, vmax=vlos_max)
      clouds_rm[:, 3] = -clouds_rm[:, 3]  # positive is recending
    
    if velocity == True:
      cax=ax1.scatter(clouds_rm[:, 1], clouds_rm[:, 2], c=clouds_rm[:, 3], cmap=cmap, norm=norm, alpha=0.7, s=10, edgecolors='k', linewidths=0.5)
    else:
      ax1.scatter(clouds_rm[:, 1], clouds_rm[:, 2], alpha=0.7, s=10, linewidths=0.5)
    
    # check the range 
    clouds_max = np.max(clouds_rm[:, 0:3])
    if np.max(np.abs(range)) < clouds_max/10.0:
      print("cloud locations significant exceed input/default range=[%.2f, %.2f], better to set a larger range."%(range[0], range[1]))

    ax1.set_xlim(range[0], range[1])
    ax1.set_ylim(range[0], range[1])
    ax1.minorticks_on()

    #set the same ticks
    xt = ax1.get_xticks()
    ax1.set_yticks(xt)

    if objname is not None:
      ax1.text(0.06, 0.95, r"\bf "+ objname, transform=ax1.transAxes, horizontalalignment='left', verticalalignment='top', bbox=dict(boxstyle="round", facecolor='w', alpha=0.7))
        
    ax1.set_xlabel(r"$y$ (lt-day)")
    ax1.set_ylabel(r"$z$ (lt-day)")
    ax1.invert_xaxis()

    # plot N and E directions
    ax1.arrow(0.95, 0.05,  0.0, 0.1, transform=ax1.transAxes, lw=2)
    ax1.arrow(0.95, 0.05, -0.1, 0.0, transform=ax1.transAxes, lw=2)
    ax1.text(0.95, 0.2, r"\bf N", va='center', ha='center', transform=ax1.transAxes)
    ax1.text(0.8, 0.05, r"\bf E", va='center', ha='center', transform=ax1.transAxes)
    
    if velocity == True:
      cbar = plt.colorbar(cax, cax=ax2, location='top', label=r"LOS Velocity (km~s$^{-1}$)")

    if format=="jpg":    
      fig.savefig("fig_clouds_los.jpg", bbox_inches='tight', dpi=500)
    else:
      fig.savefig("fig_clouds_los.pdf", bbox_inches='tight')
    
    if doshow:
      plt.show()
    else:
      plt.close()


  def postprocess(self, temperature=1, doshow=False):
    """
    do posterior process, output results into "cdnest.pdf"
    """
    print("# doing posterior process, plotting cdnest.pdf.")
    pt(int(self.param["flagdim"]), temp=temperature, fcut=0.0, fdir=self.param["filedir"], doshow=doshow)
    return
