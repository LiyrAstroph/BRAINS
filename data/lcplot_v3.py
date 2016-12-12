import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import copy

VelUnit = np.sqrt( 6.672e-8 * 1.0e6 * 1.989e33 / (2.9979e10*8.64e4)) / 1.0e5

def set_cov_Pmat(sigma, tau, alpha, Tcon):
  Pmat = np.zeros((Tcon.shape[0], Tcon.shape[0]))
  xv, yv = np.meshgrid(Tcon, Tcon)
  Pmat = sigma*sigma*np.exp( - pow( np.fabs( xv - yv )/tau, alpha) )
  return Pmat

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=15)

obj = "sim"

fp=open(obj + "_hb2d.txt", "r")
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

fp.close()
grid_wave = grid_vel / 3e5 * 4861.0 + 4861.0
dV = (grid_vel[1]-grid_vel[0]) / VelUnit

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
fp.close()


# read light curves
conlc=np.loadtxt(obj + "_con.txt")
conlc_sim=np.loadtxt("pcon.txt")
hblc=np.zeros((nt, 3))
hblc[:, 1]=np.sum(prof, axis=1) * dV
hblc[:, 2]=np.sqrt(np.sum(prof_err**2, axis=1)) * dV
hblc_sim=np.sum(prof_sim, axis=1)*dV

date0 = conlc[0, 0]
conlc[:, 0]=conlc[:, 0]-date0
conlc_sim[:, 0]=conlc_sim[:, 0]-date0
date_hb_sim = date_hb_sim-date0
date_hb = date_hb-date0 
hblc[:, 0] = date_hb

grid_vel /= 1.0e3
grid_vel_sim /= 1.0e3

# read sample
con_scale = 1.0/np.mean(conlc[:, 1])
line_scale = 1.0/np.mean(hblc[:, 1])
hd=np.loadtxt("sample2d.txt", skiprows=1)
phd = np.loadtxt("posterior_sample2d.txt")
hd_info = np.loadtxt("sample_info2d.txt", skiprows=1)
level = np.loadtxt("levels2d.txt", skiprows=1)
idx_mbh = np.where(hd_info[:, 0]>level.shape[0] - 40)
hd_sort=np.sort(hd[idx_mbh[0], 8]/np.log(10.0)+6.0)
mbh1=hd_sort[len(hd_sort)*0.1585]
mbh2=hd_sort[len(hd_sort)*(1.0-0.1585)]
print(mbh1, mbh2)

syserr_con = np.exp(np.mean(phd[:, 12]))/con_scale
syserr = np.exp(np.mean(phd[:, 11]))/line_scale

print("scale:", con_scale, line_scale)
print("syserr:", syserr_con, syserr)

logP = np.loadtxt("posterior_sample_info2d.txt", skiprows=1)
nmax = np.argmax(logP)
print(nmax, logP[nmax])

#========================================================================
fig=plt.figure(1, figsize=(9, 8))
cmap=plt.get_cmap('jet')


#========================================================================
# subfig 5
ax5=fig.add_axes([0.1, 0.2, 0.52, 0.15])

prof_rec = np.zeros((nt, nv))
prof_rec_max = np.zeros((nt, nv))
line_rec = np.zeros(nt)
line_rec_max = np.zeros(nt)
fp = open("line2d_rec.txt", "r")
chi = np.zeros(len(logP))
 
for i in range(len(logP)):
  for j in range(nt):
    line = fp.readline()
    prof_rec[j, :]=line.split()
    
  line_rec = np.sum(prof_rec, axis=1) * dV
  chi[i] = np.sum( (line_rec - hblc[:, 1]) * (line_rec - hblc[:, 1]))
  #chi[i] = np.sum( (prof_rec - prof) * (prof_rec - prof))
  fp.readline()
  
fp.close()

fp = open("line2d_rec.txt", "r")
ichimin = np.argmin(chi)
ichimin = nmax
print(ichimin, chi[ichimin])
print(ichimin, logP[ichimin])
print(nmax, phd[nmax, 11])
print(ichimin, phd[ichimin, 11])
print(nmax, phd[nmax, 12])
print(ichimin, phd[ichimin, 12])
print(ichimin, phd[ichimin, 0])

plotnum = 0
for i in range(len(logP)):
  for j in range(nt):
    line = fp.readline()
    prof_rec[j, :]=line.split()
    
  line_rec = np.sum(prof_rec, axis=1) * dV
  if i == ichimin:
    prof_rec_max = copy.copy(prof_rec)
    line_rec_max = copy.copy(line_rec)

  fp.readline()
  #if (chi[i] < chi[ichimin] * (1.0 + 0.5)) and (plotnum < 100):
  ax5.plot(date_hb, line_rec* 4861*VelUnit/3e5, color='grey', lw=0.1, alpha=0.1)
  plotnum += 1
    
fp.close()

plt.errorbar(hblc[:, 0], hblc[:, 1] * 4861*VelUnit/3e5, yerr=np.sqrt(hblc[:, 2]*hblc[:, 2] + syserr * syserr * dV*dV)* 4861*VelUnit/3e5, marker='None', markersize=3, ls='none', lw=1.0, capsize=0.7, markeredgewidth=0.5, zorder=32)
plt.plot(date_hb, line_rec_max* 4861*VelUnit/3e5, color='red')

ax5.set_xlabel(r'$\rm HJD\ (+2\ 450\ 500)$')
ax5.set_ylabel(r'$F_{\rm H\beta}$')

ymax = np.max(hblc[:, 1]* 4861*VelUnit/3e5)
ymin = np.min(hblc[:, 1]* 4861*VelUnit/3e5)
ax5.set_ylim(ymin - 0.2*(ymax-ymin), ymax + 0.2*(ymax- ymin))

xmax = hblc[-1, 0]
xmin = hblc[0, 0]
ax5.set_xlim(xmin - 0.2*(xmax-xmin), xmax + 0.1*(xmax- xmin))

#ax5.xaxis.set_minor_locator(MultipleLocator(4.0))
#ax5.yaxis.set_major_locator(MultipleLocator(0.5))
#ax5.yaxis.set_minor_locator(MultipleLocator(0.1))

#========================================================================
# subfig 4
ax4=fig.add_axes([0.1, 0.35, 0.52, 0.15])

con_scale = np.mean(conlc[:, 1])
con = np.zeros(conlc_sim.shape[0])
date = np.zeros(conlc_sim.shape[0])
chi = np.zeros(len(logP))

fp = open("con_rec.txt", "r")
for i in range(len(logP)):
  for j in range(len(con)):
    line = fp.readline()
    date[j], con[j] = line.split()
  
  confit = np.interp(conlc[:, 0], date, con)
  chi[i] = np.sum( (confit - conlc[:, 1]) * (confit - conlc[:, 1]) / conlc[:, 2]/conlc[:, 2])
  fp.readline()

fp.close()
print(chi[ichimin])

plotnum = 0
fp = open("con_rec.txt", "r")
for i in range(len(logP)):
  for j in range(len(con)):
    line = fp.readline()
    date[j], con[j] = line.split()
  if(i==ichimin):
    con_max = copy.copy(con)  
  fp.readline()
  #if (chi[i] < chi[ichimin] * (1.0 + 0.01)) and (plotnum < 100):
  ax4.plot(date-date0, con, color='grey', lw=0.1, alpha=0.1)
  plotnum += 1

fp.close()

con = copy.copy(con_max)
plt.errorbar(conlc[:, 0], conlc[:, 1], yerr=np.sqrt(conlc[:, 2]*conlc[:, 2] + syserr_con*syserr_con), marker='None', markersize=3, ls='none', lw=1.0, capsize=1, markeredgewidth=0.5, zorder=32)
plt.plot(date-date0, con, color='red')
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

xlim=ax5.get_xlim()
ax4.set_xlim([xlim[0], xlim[1]])

#========================================================================
# subfig 1
ax1 = fig.add_axes([0.1, 0.6, 0.25, 0.3])

plt.imshow(prof, cmap=cmap, interpolation=None, aspect='auto', extent=[grid_vel[0], grid_vel[nv-1], date_hb[-1], date_hb[0]], vmax = np.amax(prof), vmin=np.amin(prof))
ax1.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
ax1.set_ylabel(r'$\rm Time\ (+2\ 450\ 500)$')

xlim=ax1.get_xlim()
ylim=ax1.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Data$', color='white')

#ax1.xaxis.set_major_locator(MultipleLocator(2000))
#ax1.xaxis.set_minor_locator(MultipleLocator(400))

#========================================================================
# subfig 2
ax2=fig.add_axes([0.37, 0.6, 0.25, 0.3])

plt.imshow(prof_rec_max, cmap=cmap, interpolation=None,  aspect='auto', extent=[grid_vel[0], grid_vel[nv-1], date_hb[-1], date_hb[0]], vmax = np.amax(prof_rec_max), vmin=np.amin(prof_rec_max))
ax2.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
#ax2.set_ylabel('Time (+2 450 000)')
xlim=ax2.get_xlim()
ylim=ax2.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Model$', color='white')


[i.set_visible(False) for i in ax2.get_yticklabels()]
#ax2.xaxis.set_major_locator(MultipleLocator(2000))
#ax2.xaxis.set_minor_locator(MultipleLocator(400))

#========================================================================
# subfig 3
ax3=fig.add_axes([0.7, 0.6, 0.25, 0.3])

chifit = np.zeros(nt)
for i in range(nt):
  chifit[i] = np.sum( (prof_rec_max[i, :] - prof[i, :]) * (prof_rec_max[i, :] - prof[i, :]) )

chifit_sort = np.sort(chifit)

offset = np.max(prof.flatten()) * 0.3

idx = np.where(chifit == chifit_sort[0])
i = idx[0][0]
j = 0
print(i)
plt.errorbar(grid_vel, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr*syserr), ls='none', ecolor='k', capsize=1, markeredgewidth=1, zorder=32)
plt.plot(grid_vel, prof_rec_max[i, :]+j*offset, color='b', lw=2)

idx = np.where(chifit == chifit_sort[1])
i = idx[0][0]
j = 1
print(i)
plt.errorbar(grid_vel, prof[i, :]+j*offset, yerr=np.sqrt(prof_err[i, :]*prof_err[i, :] + syserr*syserr), ls='none', ecolor='k', capsize=1, markeredgewidth=1)
plt.plot(grid_vel, prof_rec_max[i, :]+j*offset, color='b', lw=2)

ax3.set_xlabel(r'$\rm Velocity\ (10^3km\ s^{-1})$')
ax3.set_ylabel(r'$\rm Flux$')
ax3.set_xlim([grid_vel[0], grid_vel[-1]])
xlim=ax3.get_xlim()
ylim=ax3.get_ylim()
plt.text(xlim[0]+0.08*(xlim[1]-xlim[0]), ylim[1]-0.1*(ylim[1]-ylim[0]), r'$\rm Profile$', color='black')

#ax3.xaxis.set_major_locator(MultipleLocator(2000))
#ax3.xaxis.set_minor_locator(MultipleLocator(400))

#========================================================================
# subfig 6
ax6=fig.add_axes([0.7, 0.2, 0.25, 0.3])

hp, bins, patches=plt.hist(hd[idx_mbh[0], 8]/np.log(10.0)+6.0, bins=10, normed=1)
ax6.set_xlabel(r'$\log(M_\bullet/M_\odot)$')
ax6.set_ylabel(r'$\rm Hist$')
ylim=ax6.get_ylim()

ax6.set_ylim(ylim)
#ax6.set_xlim((6.0, 8.5))
#ax6.xaxis.set_major_locator(MultipleLocator(1.0))


plt.savefig('lc.pdf', bbox_inches='tight')
#plt.show()
