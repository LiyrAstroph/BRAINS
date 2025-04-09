import sys
import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

__all__ = ["postprocess"]

def logsumexp(values):     # log SUM( exp(values) )
  biggest = np.max(values)
  x = values - biggest
  result = np.log(np.sum(np.exp(x))) + biggest
  return result

def logdiffexp(x1, x2):    # log( exp(x1) - exp(x2))
  biggest = x1
  xx1 = x1 - biggest
  xx2 = x2 - biggest
  result = np.log(np.exp(xx1) - np.exp(xx2)) + biggest
  return result

def postprocess(ndim, temp=1.0, fcut=0.0, fdir="../", doshow=False):

  pdf = PdfPages('cdnest.pdf')

  cut = fcut
  numResampleLogX=1
  compression_bias_min=1.
  compression_scatter=0.
  temperature = temp
  zoom_in = True
  moreSamples = 1
 
  if ndim == 0:
    str_dim = '_con'
  elif ndim == 1:
    str_dim = '_1d'
  elif ndim == 2:
    str_dim = '_2d'
  elif ndim == 3:
    str_dim = '_lp'
  elif ndim == 4:
    str_dim = '_sa'
  elif ndim == 5:
    str_dim = '_sa1d'
  elif ndim == 6:
    str_dim = '_sa2d'
  elif ndim == 7:
    str_dim = '_sarm'
  else:
    raise Exception('incorrect dimension.')
  
  print(fdir)
  
  try:
    levels_orig = np.loadtxt(fdir+"/data/levels"+str_dim+".txt", comments='#')
  except:
    levels_orig = np.genfromtxt(fdir+"/data/levels"+str_dim+".txt", comments='#', skip_footer=1)

  try:  
    sample_info = np.loadtxt(fdir+"/data/sample_info"+str_dim+".txt", comments='#')
  except:
    sample_info = np.genfromtxt(fdir+"/data/sample_info"+str_dim+".txt", comments='#', skip_footer=1)

  try:
    sample = np.atleast_2d(np.loadtxt(fdir+"/data/sample"+str_dim+".txt"))
  except:
    sample = np.atleast_2d(np.genfromtxt(fdir+"/data/sample"+str_dim+".txt", skip_footer=1))
  
  sample = sample[int(cut*sample.shape[0]):, :]
  sample_info = sample_info[int(cut*sample_info.shape[0]):, :]
  
  if sample.shape[0] != sample_info.shape[0]:
    print('# Size mismatch. Truncating...')
  lowest = np.min([sample.shape[0], sample_info.shape[0]])
  sample = sample[0:lowest, :]
  sample_info = sample_info[0:lowest, :]

  # reset level assignment
  idx = (sample_info[:, 0] > levels_orig.shape[0] - 1)
  sample_info[idx, 0] = levels_orig.shape[0] - 1
  
  fig = plt.figure(figsize=(12, 8))
  fig.subplots_adjust(wspace=0.2, hspace=0.3)
  ax1 = fig.add_subplot(3, 1, 1)
  ax1.plot(sample_info[:,0], "k", lw=1, zorder=0)
  ax1.set_xlabel("Iteration")
  ax1.set_ylabel("Level")

  # resterization
  ax1.set_rasterization_zorder(0)
  
  ax2=fig.add_subplot(3,2,3)
  ax2.plot(np.diff(levels_orig[:,0]), "k")
  ax2.set_ylabel("Compression")
  ax2.set_xlabel("Level")
  xlim = plt.gca().get_xlim()
  ax2.axhline(-1., color='g')
  ax2.axhline(-np.log(10.), color='g', linestyle="--")
  ax2.set_ylim(ymax=0.05)
  
  ax3=fig.add_subplot(3,2,5)
  good = np.nonzero(levels_orig[:,4] > 0)[0]
  ax3.plot(levels_orig[good,3]/levels_orig[good,4], marker='o')
  ax3.set_xlim(xlim)
  ax3.set_ylim([0., 1.])
  ax3.set_xlabel("Level")
  ax3.set_ylabel("MH Acceptance")
  
  logl_levels = [(levels_orig[i,1], levels_orig[i, 2]) for i in range(0, levels_orig.shape[0])]
  logl_samples = [(sample_info[i, 1], sample_info[i, 2], i) for i in range(0, sample.shape[0])]
  logx_samples = np.zeros((sample_info.shape[0], numResampleLogX))
  logp_samples = np.zeros((sample_info.shape[0], numResampleLogX))
  logP_samples = np.zeros((sample_info.shape[0], numResampleLogX))
  P_samples = np.zeros((sample_info.shape[0], numResampleLogX))
  logz_estimates = np.zeros((numResampleLogX, 1))
  H_estimates = np.zeros((numResampleLogX, 1))
  
  sandwich = sample_info[:,0].copy().astype('int')
  for i in range(0, sample.shape[0]):
    while sandwich[i] < levels_orig.shape[0] -1 and logl_samples[i] > logl_levels[sandwich[i] + 1]:
      sandwich[i] += 1
  
  for z in range(0, numResampleLogX):
    levels = levels_orig.copy()
    compressions = -np.diff(levels[:, 0])
    compressions *= compression_bias_min + ( 1.0 - compression_bias_min) * np.random.rand()
    compressions *= np.exp(compression_scatter * np.random.randn(compressions.size))
    levels[1:, 0] = -compressions
    levels[:, 0] = np.cumsum(levels[:, 0])
  
    for i in range(0, levels.shape[0]):
      which = np.nonzero(sandwich == i)[0]
      logl_samples_thisLevel = []
      for j in range(0, len(which)):
        logl_samples_thisLevel.append(copy.deepcopy(logl_samples[which[j]]))
      logl_samples_thisLevel = sorted(logl_samples_thisLevel)
  
      N = len(logl_samples_thisLevel)
  
      logx_max = levels[i, 0]
      if i == levels.shape[0] - 1:
        logx_min = -1E300   # the minmial X -> 0
      else:
        logx_min = levels[i+1, 0]
  
      Umin = np.exp(logx_min - logx_max)
  
      if N ==0 or numResampleLogX > 1:
        U = Umin + (1.0 - Umin) * np.random.rand(len(which))
      else:
        U = Umin + (1.0 - Umin) *np.linspace(1.0/(N + 1), 1.0 - 1.0/(N+1), N)
  
      logx_samples_thisLevel = np.sort(logx_max + np.log(U))[::-1]
  
      for j in range(0, which.size):
        logx_samples[logl_samples_thisLevel[j][2]][z] = logx_samples_thisLevel[j]
  
        if j != which.size -1:
          left = logx_samples_thisLevel[j+1]
        elif i == levels.shape[0] -1:
          left = -1E300
        else:
          left = levels[i+1][0]
  
        if j != 0:
          right = logx_samples_thisLevel[j-1]
        else:
          right = levels[i][0]
  
        logp_samples[logl_samples_thisLevel[j][2]][z] = np.log(0.5) + logdiffexp(right, left)
  
    logl = sample_info[:,1]/temperature
    
    logp_samples[:,z] = logp_samples[:,z] - logsumexp(logp_samples[:,z]) # make sure that Sum(p) = 1.0
    logP_samples[:,z] = logp_samples[:,z] + logl
    logz_estimates[z] = logsumexp(logP_samples[:,z])
    logP_samples[:,z] -= logz_estimates[z]
    P_samples[:,z] = np.exp(logP_samples[:,z])
    H_estimates[z] = -logz_estimates[z] + np.sum(P_samples[:,z]*logl)
  
    ax4=fig.add_subplot(3,2,4)
    ax4.plot(logx_samples[:,z], sample_info[:,1], 'b.', label='Samples', zorder=0)
    ax4.plot(levels[1:,0], levels[1:,1], 'r.', label='Levels')
    ax4.legend(numpoints=1, loc='lower left')
    ax4.set_ylabel('log(L)')
    ax4.set_xlabel('log(X)')

    # resterization
    ax4.set_rasterization_zorder(0)

    ax1.set_title(str(z+1) + "/" + str(numResampleLogX) + ", log(Z) = " + str(logz_estimates[z][0]))
    
    # Use all plotted logl values to set ylim
    combined_logl = np.hstack([sample_info[:,1], levels[1:, 1]])
    combined_logl = np.sort(combined_logl)
    lower = combined_logl[int(0.1*combined_logl.size)]
    upper = combined_logl[-1]
    diff = upper - lower
    lower -= 0.05*diff
    upper += 0.05*diff
    if zoom_in:
      ax4.set_ylim([lower, upper])
  
    xlim = ax4.get_xlim()
  
    ax5=fig.add_subplot(3,2,6)
    ax5.plot(logx_samples[:,z], P_samples[:,z], 'b.', zorder=0)
    ax5.set_ylabel('Posterior Weights')
    ax5.set_xlabel('log(X)')
    ax5.set_xlim(xlim)
    
    xlim = ax5.get_xlim()
    ylim = ax5.get_ylim()
    ax5.text(xlim[0]+0.1*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), 'T='+str(temperature))
    
    # resterization
    ax5.set_rasterization_zorder(0)

    fig.align_ylabels()

    pdf.savefig()
    
    if z < numResampleLogX-1:
      ax4.clear()
      ax4.remove()
      ax5.clear()
      ax5.remove()
  
  if doshow:
    plt.show()
  else:
    plt.close()

  pdf.close()
  P_samples = np.mean(P_samples, 1)
  P_samples = P_samples/np.sum(P_samples)
  logz_estimate = np.mean(logz_estimates)
  logz_error = np.std(logz_estimates)
  H_estimate = np.mean(H_estimates)
  H_error = np.std(H_estimates)
  ESS = np.exp(-np.sum(P_samples*np.log(P_samples+1E-300)))
  
  print("log(Z) = " + str(logz_estimate) + " +- " + str(logz_error))
  print("Information = " + str(H_estimate) + " +- " + str(H_error) + " nats.")
  print("Effective sample size = " + str(ESS))
  
  N = int(moreSamples*ESS)
  posterior_sample = np.zeros((N, sample.shape[1]))
  w = P_samples
  w = w/np.max(w)
  
  #np.savetxt('../data/weigths'+str_dim+'.txt', w)
  
  for i in range(0, N):
    while True:
      which = np.random.randint(sample.shape[0])
      if np.random.rand() < w[which]:
        break
    posterior_sample[i, :] = sample[which, :]
  
  #np.savetxt('../data/posterior_sample'+str_dim+'.txt', posterior_sample)


if __name__ == "__main__":
  temp = 1.0
  fcut = 0.0
  assert len(sys.argv) >= 2, "No dimension specified."
  
  if len(sys.argv) >= 3:
    temp = float(sys.argv[2])

  if len(sys.argv) >= 4:
    fcut = float(sys.argv[3])

  postprocess(int(sys.argv[1]), temp=temp, fcut=fcut)
