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
import corner
import numpy as np 
import configparser as cp
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
      self.option_file = file_dir + "/param/OPTIONSSA"
    elif flag_dim == '4':
      self.option_file = file_dir + "/param/OPTIONSSA1D"
    elif flag_dim == '5':
      self.option_file = file_dir + "/param/OPTIONSSA2D"
    elif flag_dim == '6':
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
      self.para_names_file = file_dir + "/data/para_names_sa.txt"
    elif flag_dim == '4':
      self.para_names_file = file_dir + "/data/para_names_sa1d.txt"
    elif flag_dim == '5':
      self.para_names_file = file_dir + "/data/para_names_sa2d.txt"
    elif flag_dim == '6':
      self.para_names_file = file_dir + "/data/para_names_sarm.txt"
    else:
      raise Exception("Incorrect FlagDim.")

    names = np.genfromtxt(self.para_names_file, skip_header=7, \
      dtype=[int, 'S30', float, float, int, int, float], delimiter=[4, 30, 12, 12, 4, 5, 15])
    
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
      if re.match("BLR model", name):
        self.num_param_blrmodel_rm += 1
      
      if re.match("SA BLR model", name):
        self.num_param_blrmodel_sa += 1
      
      if re.match("SA Extra Par", name):
        self.num_param_sa_extra += 1
    
    self.num_param_sa = self.num_param_blrmodel_sa + self.num_param_sa_extra
    idx_con =  np.nonzero(self.para_names['name'] == 'sys_err_con')
    self.num_param_rm_extra = idx_con[0][0] - self.num_param_blrmodel_rm - self.num_param_sa
    self.num_param_blr_rm = self.num_param_rm_extra + self.num_param_blrmodel_rm

    if 'time series' in self.para_names['name']:
      idx = np.nonzero(self.para_names['name'] == 'time series')
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
      if re.match("BLR model ln\(Mbh\)", self.para_names['name'][i]):
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
    self.con_scale = 1.0
    self.line_scale = 1.0

    self.file_dir = self.param['filedir'] +"/"
    self._load_data()
    self._load_results() 

    self.postfix=["con", "1d","2d", "sa", "sa1d", "sa2d", "sarm"]

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

      # calculate con_scale and line_scale
      self.con_scale = self.data['con_data'].shape[0]/np.sum(self.data['con_data'][:, 1])
    
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

  def plot_drw_parameters(self):
    """
    plot drw paramete distributions
    """

    if int(self.param['flagdim']) in [0, 1, 2, 4, 5, 6]:
      idx = self.num_param_blr_rm+self.num_param_sa
      if np.std(self.results['sample'][:, idx]) == 0.0:
        fig = corner.corner(self.results['sample'][:, idx+1:idx+3], labels=['ln sigma', 'ln tau'], smooth=True)
      else:
        fig = corner.corner(self.results['sample'][:, idx:idx+3], labels=['syserr', 'ln sigma', 'ln tau'], smooth=True)
      plt.show()
      return fig
    else:
      print("The running does not have continnum reconstruction.")
  
  def get_narrow_line(self, grid_vel, imax, iepoch):
    """
    get narrow line component
    """
    
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

    if int(self.param["flagnarrowline"]) == 1:
      flux = float(self.param["fluxnarrowline"])
      width = width_na
      shift = shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      prof = flux * np.exp( -0.5 *(grid_vel - shift)**2/width_br**2)
      return prof
    elif int(self.param["flagnarrowline"]) == 2:
      flux_na = float(self.param["fluxnarrowline"])
      flux_na_err = float(self.param["fluxnarrowlineerr"])
      flux = self.results['sample'][imax, self.num_param_blrmodel_rm] * flux_na_err + flux_na
      width = self.results['sample'][imax, self.num_param_blrmodel_rm+1] * width_na_err + width_na 
      shift = self.results['sample'][imax, self.num_param_blrmodel_rm+2] * shift_na_err + shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      prof = flux * np.exp( -0.5 *(grid_vel - shift)**2/width_br**2)
      return prof
    elif int(self.param["flagnarrowline"]) == 3:
      flux = np.exp(self.results['sample'][imax, self.num_param_blrmodel_rm])
      width = self.results['sample'][imax, self.num_param_blrmodel_rm+1] * width_na_err + width_na 
      shift = self.results['sample'][imax, self.num_param_blrmodel_rm+2] * shift_na_err + shift_na
      
      width_br = np.sqrt(width**2 + broad**2)
      flux *= width/width_br
      # note line scale
      flux /= self.line_scale
      prof = flux * np.exp(-0.5 *(grid_vel - shift)**2/width_br**2)
      return prof

      
  def plot_results_con(self):
    """
    plot continuum results
    """

    if not int(self.param['flagdim']) in [0]:
      print("Flagdim =", self.param['flagdim'], "not continuum analysis.")
      return
    
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

    fig.savefig("results_con.pdf", bbox_inches='tight')
    plt.close()

  def plot_results_2d_style2018(self):
    """
    plot 2d results in the style of 2018 ApJ paper
    """

    if not int(self.param['flagdim']) in [2, 5, 6]:
      print("Flagdim =", self.param['flagdim'], "no 2D RM data.")
      return
    
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

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=15)

    #========================================================================
    fig=plt.figure(1, figsize=(9, 8))
    cmap=plt.get_cmap('jet')
    
    
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
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
                 yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + syserr_line * syserr_line * dV*dV)* wave0*self.VelUnit/3e5, \
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
    plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    ax4.plot(con_date-date0, con_mean, color='red', zorder=20) 
    ax4.fill_between(con_date-date0, y1=con_mean_low, y2=con_mean_upp, color='grey')

    ax4.set_ylabel(r'$F_{\rm 5100}$')

    [i.set_visible(False) for i in ax4.get_xticklabels()]

    # set a proper range
    xmin = min(conlc[0, 0], date_line[0])
    xmax = max(conlc[-1, 0], date_line[-1])
    dx = xmax-xmin
    xmin -= 0.05*dx 
    xmax += 0.05*dx 
    ax4.set_xlim(xmin, xmax)
    ax5.set_xlim(xmin, xmax)
    ax4.minorticks_on()
    ax5.minorticks_on()
    
    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
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

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
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
                 ls='none', ecolor='k', capsize=1, markeredgewidth=1, zorder=10)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[0], lw=2, zorder=1)

    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[0], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', label='Narrow Line', zorder=1)
    
    idx = np.where(chifit == chifit_sort[1])
    i = idx[0][0]
    j = 1
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='k', capsize=1, markeredgewidth=1, zorder=10)
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
    plt.show()
  
  def plot_results_2d_style2022(self):
    """
    plot 2d results in the style of 2022 ApJ paper
    """

    if not int(self.param['flagdim']) in [2, 5, 6]:
      print("Flagdim =", self.param['flagdim'], "no 2D RM data.")
      return
    
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
                 yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + syserr_line * syserr_line * dV*dV)* wave0*self.VelUnit/3e5, \
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
    plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), \
                 marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
    
    con_rec = self.results['con_rec']
    con_date = con_rec[0, :, 0]
    con_mean = np.quantile(con_rec[:, :, 1], axis=0, q=0.5)
    con_mean_upp = np.quantile(con_rec[:, :, 1], axis=0, q=(1.0-0.683)/2.0)
    con_mean_low = np.quantile(con_rec[:, :, 1], axis=0, q=1.0 - (1.0-0.683)/2.0)
    
    ax5.plot(con_date-date0, con_mean, color='red', zorder=20) 
    ax5.fill_between(con_date-date0, y1=con_mean_low, y2=con_mean_upp, color='grey')

    ax5.set_ylabel(r'$F_{\rm 5100}$')

    [i.set_visible(False) for i in ax5.get_xticklabels()]

    # set a proper range
    xmin = min(conlc[0, 0], date_line[0])
    xmax = max(conlc[-1, 0], date_line[-1])
    dx = xmax-xmin
    xmin -= 0.05*dx 
    xmax += 0.05*dx 
    ax5.set_xlim(xmin, xmax)
    ax5.tick_params(labelright=True, labelleft=False)
    ax5.yaxis.set_label_position('right')
    
    xlim = ax5.get_xlim()
    ax6.set_xlim(xlim[0], xlim[1])
    ax5.minorticks_on()
    ax6.minorticks_on()

    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
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
    cmap=ax3.imshow(img,  aspect='auto', origin='lower', cmap='jet',
                    extent=[grid_vel[0]/1.0e3, grid_vel[nv-1]/1.0e3, 1, prof.shape[0]], vmax=7.0, vmin = -7.0)
    
    plt.colorbar(cmap, ticks=[-5, -2.5, 0.0, 2.5, 5.0])
    ax3.text(0.08, 0.9, r'\bf Residuals', color='white', transform=ax3.transAxes)
    ax3.set_xlabel(r'$\rm Velocity\ (10^3\ km\ s^{-1})$')
    [i.set_visible(False) for i in ax3.get_yticklabels()]

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

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
                 ls='none', ecolor='k', capsize=1, markeredgewidth=1)
    plt.plot(grid_vel/1.0e3, prof_rec_max[i, :]+j*offset, color=cycle[0], lw=2)
    
    if int(self.param["flagnarrowline"]) > 0:
      prof_na = self.get_narrow_line(grid_vel, imax, i)
      plt.plot(grid_vel/1.0e3, (prof_rec_max[i, :]-prof_na)+j*offset, lw=1, color=cycle[0], ls='--', zorder=1)
      plt.plot(grid_vel/1.0e3, prof_na+j*offset, lw=1, color='k', ls='--', zorder=1, label='Narrow Line')

    idx = np.where(chifit == chifit_sort[1])
    i = idx[0][0]
    j = 1
    plt.errorbar(grid_vel/1.0e3, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr_line*syserr_line), \
                 ls='none', ecolor='k', capsize=1, markeredgewidth=1)
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
    plt.show()
  
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
  
  def plot_tran2d(self, tau_range=None, vel_range=None):
    """
    plot transfer function with the maximum prob
    """

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
      plt.show()
      fig.savefig("tran2d.pdf", bbox_inches='tight')
    else:
      print("Flagdim =", self.param["flagdim"], "no trans2d")
      return
  
  def plot_tran1d(self, tau_range=None):
    """
    plot transfer function 1d
    """
    if int(self.param['flagdim']) in [1, 4]:
      tran_rec = self.results['tran_rec']
      fig = plt.figure()
      ax = fig.add_subplot(111)
      for i in range(0, tran_rec.shape[0], int(tran_rec.shape[0]/100+1.0)):
        ax.plot(tran_rec[i, :, 0], tran_rec[i, :, 1], lw=0.5)
      
      if tau_range is not None:
        ax.set_xlim(tau_range[0], tau_range[1])

      plt.show()
    elif int(self.param['flagdim']) in [2, 5, 6]:
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
      plt.show()
    else:
      print("Flagdim =", self.param["flagdim"], "no trans1d")
      return
  
  def plot_blrmodel_para_hist(self, para_indx=None):
    """
    plot histograms of BLR model parameters
    """
    
    prange = np.zeros((self.results['sample'].shape[1], 2))
    prange[:, 0] = np.min(self.results['sample'], axis=0)
    prange[:, 1] = np.max(self.results['sample'], axis=0)
    dprange = prange[:, 1]-prange[:, 0]
    prange[:, 0] -= 0.1*dprange 
    prange[:, 1] += 0.1*dprange

    if int(self.param['flagdim']) >= 1:
      if para_indx is None:
        fig = corner.corner(self.results['sample'][:, :self.num_param_blrmodel_rm], \
                            range = prange[:self.num_param_blrmodel_rm, :], \
                            smooth=True, smooth1d=True)
      else:
        fig = corner.corner(self.results['sample'][:, para_indx], range = prange[para_indx, :], \
                            smooth=True, smooth1d=True)
      fig.savefig("hist.pdf", bbox_inches='tight')
      plt.show()
    else:
      print("Flagdim =", self.param["flagdim"], "no BLR parameters")
      return
  
  def postprocess(self, temperature=1):
    """
    do posterior process, output results into "cdnest.pdf"
    """
    pt(int(self.param["flagdim"]), temp=temperature, fcut=0.0, fdir=self.param["filedir"])
    return
