import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter

def set_cov_Pmat(sigma, tau, alpha, Tcon):
  Pmat = np.zeros((Tcon.shape[0], Tcon.shape[0]))
  xv, yv = np.meshgrid(Tcon, Tcon)
  Pmat = sigma*sigma*np.exp( - pow( np.fabs( xv - yv )/tau, alpha) )
  return Pmat

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=15)


fp=open("1h0323_hb2d.txt", "r")
line=fp.readline()
text=line.split()
nt=int(text[1])
nv=int(text[2])

print(nv, nt)

# read data
date_hb=np.zeros(nt)
grid_vel=np.zeros(nv)
prof = np.zeros((nt, nv))
prof_err = np.zeros((nt, nv))

for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel[i], date_hb[j], prof[j, i], prof_err[j, i] = line.split()
  
  line=fp.readline()

#grid_vel = grid_vel +160.0

  
fp.close()
grid_wave = grid_vel / 3e5 * 4861.0 + 4861.0


# read sim
fp=open("pline2d_data.txt", "r")

date_hb_sim=np.zeros(nt)
grid_vel_sim=np.zeros(nv)
prof_sim = np.zeros((nt, nv))
prof_err_sim = np.zeros((nt, nv))

#fp.readline()
for j in range(0, nt):
  for i in range(0, nv):
    line=fp.readline()
    grid_vel_sim[i], date_hb_sim[j], prof_sim[j, i] = line.split()
  
  line=fp.readline()



# read light curves
conlc=np.loadtxt("1h0323_con.txt")
conlc_sim=np.loadtxt("pcon.txt")
hblc=np.zeros((nt, 3))
hblc[:, 1]=np.sum(prof, axis=1) * (grid_vel[1]-grid_vel[0]) * 4861/3e5
hblc[:, 2]=np.sqrt(np.sum(prof_err**2, axis=1)) * (grid_vel[1]-grid_vel[0]) * 4861/3e5
hblc_sim=np.sum(prof_sim, axis=1)*(grid_vel_sim[1]-grid_vel_sim[0]) * 4861/3e5

date0 = conlc[0, 0]
conlc[:, 0]=conlc[:, 0]-date0
conlc_sim[:, 0]=conlc_sim[:, 0]-date0
date_hb_sim = date_hb_sim-date0
date_hb = date_hb-date0 
hblc[:, 0] = date_hb

grid_vel /= 1.0e3
grid_vel_sim /= 1.0e3

# read pt
con_scale = np.mean(conlc[:, 1])
line_scale = np.mean(hblc[:, 0])
syserr_con = np.exp(-25.049281) * con_scale
syserr = np.exp(-6.421502) * line_scale
hd=np.loadtxt("sample2d.txt", skiprows=1)
phd = np.loadtxt("posterior_sample2d.txt")
hd_info = np.loadtxt("sample_info2d.txt", skiprows=1)
level = np.loadtxt("levels2d.txt", skiprows=1)
idx = np.where(hd_info[:, 0]>level.shape[0] - 50)
hd_sort=np.sort(hd[idx[0], 8]/np.log(10.0)+6.0)
mbh1=hd_sort[len(hd_sort)*0.1585]
mbh2=hd_sort[len(hd_sort)*(1.0-0.1585)]
print(mbh1, mbh2)


fig=plt.figure(1, figsize=(9, 8))
cmap=plt.get_cmap('jet')

ax1 = fig.add_axes([0.1, 0.6, 0.25, 0.3])

plt.imshow(prof, cmap=cmap, interpolation=None, aspect='auto', extent=[grid_vel[0], grid_vel[nv-1], date_hb[-1], date_hb[0]], vmax = np.amax(prof), vmin=np.amin(prof))
ax1.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
ax1.set_ylabel(r'$\rm Time\ (+2\ 450\ 500)$')

xlim=ax1.get_xlim()
ylim=ax1.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Data$', color='white')

#ax1.xaxis.set_major_locator(MultipleLocator(2000))
#ax1.xaxis.set_minor_locator(MultipleLocator(400))


ax2=fig.add_axes([0.37, 0.6, 0.25, 0.3])

plt.imshow(prof_sim, cmap=cmap, interpolation=None,  aspect='auto', extent=[grid_vel[0], grid_vel[nv-1], date_hb[-1], date_hb[0]], vmax = np.amax(prof_sim), vmin=np.amin(prof_sim))
ax2.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
#ax2.set_ylabel('Time (+2 450 000)')
xlim=ax2.get_xlim()
ylim=ax2.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Model$', color='white')


[i.set_visible(False) for i in ax2.get_yticklabels()]
#ax2.xaxis.set_major_locator(MultipleLocator(2000))
#ax2.xaxis.set_minor_locator(MultipleLocator(400))


ax3=fig.add_axes([0.7, 0.6, 0.25, 0.3])

offset = np.max(prof.flatten()) * 0.3
for j in range(0, 2):
  i = np.random.randint(nt)
  plt.errorbar(grid_vel, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr*syserr/10000), ls='none', ecolor='k', capsize=1, markeredgewidth=1)
  plt.plot(grid_vel_sim, prof_sim[i, :]+j*offset, color='b', lw=2)

ax3.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
ax3.set_ylabel(r'$\rm Flux$')
ax3.set_xlim([grid_vel[0], grid_vel[-1]])
xlim=ax3.get_xlim()
ylim=ax3.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Profile$', color='black')

#ax3.xaxis.set_major_locator(MultipleLocator(2000))
#ax3.xaxis.set_minor_locator(MultipleLocator(400))


ax4=fig.add_axes([0.1, 0.35, 0.52, 0.15])

offset = 12
con_scale = np.mean(conlc[:, 1])
con = np.zeros(conlc_sim.shape[0])
for i in range(0, phd.shape[0], int(phd.shape[0]/10.0+1.0)):
  if(phd.ndim == 1):
    Pmat = set_cov_Pmat(np.exp(phd[1+offset]), np.exp(phd[2+offset]), 1.0, conlc_sim[:, 0])
    break
  else:
    Pmat = set_cov_Pmat(np.exp(phd[i, 1+offset]), np.exp(phd[i, 2+offset]), 1.0, conlc_sim[:, 0])
  Mmat = np.linalg.cholesky(Pmat)
  #Mmat = Mmat.T

  con = np.matmul(Mmat, phd[i, 4+offset:]) + phd[i, 3+offset]
  ax4.plot(conlc_sim[:, 0], con*con_scale, color='grey', linewidth=0.1)

plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), marker='o', markersize=3, ls='none', lw=1, capsize=1, markeredgewidth=0.5)
plt.plot(conlc_sim[:, 0], conlc_sim[:, 1], color='red')
#plt.plot(conlc_sim[:, 0], conlc_sim[:, 1]+conlc_sim[:, 2], ls='dotted', color='red')
#plt.plot(conlc_sim[:, 0], conlc_sim[:, 1]-conlc_sim[:, 2], ls='dotted', color='red')
ax4.set_xlim(conlc_sim[0, 0], conlc_sim[-1, 0])
ymax = np.max(conlc[:, 1])
ymin = np.min(conlc[:, 1])
ax4.set_ylim(ymin - 0.2*(ymax-ymin), ymax + 0.2*(ymax- ymin))
ax4.set_ylabel(r'$F_{\rm 5100}$')
#ax4.set_ylim([30, 95])

#ax4.xaxis.set_minor_locator(MultipleLocator(4.0))
#ax4.yaxis.set_major_locator(MultipleLocator(40.0))
#ax4.yaxis.set_minor_locator(MultipleLocator(8.0))


[i.set_visible(False) for i in ax4.get_xticklabels()]

ax5=fig.add_axes([0.1, 0.2, 0.52, 0.15])

plt.errorbar(hblc[:, 0], hblc[:, 1], yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + syserr*syserr/10000), marker='o', markersize=3, ls='none', lw=1, capsize=1, markeredgewidth=0.5)
plt.plot(date_hb_sim, hblc_sim, color='red')
#plt.plot(hblc_sim[:, 0], hblc_sim[:, 1]+hblc_sim[:, 2], ls='dotted', color='red')
#plt.plot(hblc_sim[:, 0], hblc_sim[:, 1]-hblc_sim[:, 2], ls='dotted', color='red')

ax5.set_xlabel(r'$\rm HJD\ (+2\ 450\ 500)$')
ax5.set_ylabel(r'$F_{\rm H\beta}$')
xlim=ax4.get_xlim()
ax5.set_xlim([xlim[0], xlim[1]])
#ax4.set_xlim([0, 110])

#ax5.xaxis.set_minor_locator(MultipleLocator(4.0))
#ax5.yaxis.set_major_locator(MultipleLocator(0.5))
#ax5.yaxis.set_minor_locator(MultipleLocator(0.1))
  

ax6=fig.add_axes([0.7, 0.2, 0.25, 0.3])

hp, bins, patches=plt.hist(hd[idx[0], 8]/np.log(10.0)+6.0, bins=10, normed=1)
ax6.set_xlabel(r'$\log(M_\bullet/M_\odot)$')
ax6.set_ylabel(r'$\rm Hist$')
ylim=ax6.get_ylim()

ax6.set_ylim(ylim)
#ax6.set_xlim((6.0, 8.5))
ax6.xaxis.set_major_locator(MultipleLocator(1.0))


plt.savefig('lc.eps', format='eps', bbox_inches='tight')
#plt.show()
