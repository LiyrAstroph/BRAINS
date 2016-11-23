##
# \file plotcheck.py
# \brief check if the recovered light curves are correct.
#

import numpy as np 
import matplotlib.pyplot as plt 
import copy 

ncheck = 600

logP = np.loadtxt("../data/posterior_sample_info2d.txt", skiprows=1)

ncon = 100
nline = 25
nv = 25
ntau = 100

con = np.zeros(ncon)
date_con = np.zeros(ncon)
con_check = np.zeros(ncon)

fp = open("con_rec.txt", "r")
for i in range(len(logP)):
  for j in range(len(con)):
    line = fp.readline()
    date_con[j], con[j] = line.split()
  
  fp.readline()

  if i == ncheck:
    con_check = copy.copy(con)

fp.close()


prof_rec = np.zeros((nline, nv))
prof_rec_max = np.zeros((nline, nv))

prof_check = np.zeros((nline, nv)) 
fp = open("line2d_rec.txt", "r")
 
for i in range(len(logP)):
  for j in range(nline):
    line = fp.readline()
    prof_rec[j, :]=line.split()
    
  fp.readline()
  
  if i == ncheck:
    prof_check = copy.copy(prof_rec)

fp.close()


tran = np.zeros((ntau, nv))
fp = open("../data/tran2d_rec.txt", "r")
for i in range(len(logP)):
  for j in range(0, ntau):
    line=fp.readline()
    tran[j, :] = line.split()
  fp.readline()

  if i == ncheck:
    tran_check = copy.copy(tran)

fp.close()

tran1d = np.sum(tran_check, axis = 1)

date_line = np.loadtxt("../data/mrk493_hb.txt", usecols=(0,))

line1d_check = np.sum(prof_check, axis = 1)

fig = plt.figure(0)

plt.plot(date_line, line1d_check)

line1d = np.zeros(nline)


grid_tau = np.linspace(0, 50.0, ntau)
print(grid_tau)
for i in range(nline):
  tl = date_line[i]
  line1d[i] = 0.0
  for j in range(ntau):
    tau = grid_tau[j]
    tc = tl-tau
    if(tc > date_con[0]):
      fcon = np.interp(tc, date_con, con_check)
      line1d[i] = line1d[i] +  pow(fcon, 1.0-0.2) * tran1d[j] 

line1d = line1d * (grid_tau[1]-grid_tau[0]) 
plt.plot(date_line, line1d)

date_line_rec = np.linspace(date_line[0], date_line[-1], 100)
line1d_rec = np.zeros(100)
for i in range(100):
  tl = date_line_rec[i]
  line1d_rec[i] = 0.0
  for j in range(ntau):
    tau = grid_tau[j]
    tc = tl-tau
    if(tc > date_con[0]):
      fcon = np.interp(tc, date_con, con_check)
      line1d_rec[i] = line1d_rec[i] +  pow(fcon, 1.0-0.2) * tran1d[j] 
      
line1d_rec = line1d_rec * (grid_tau[1]-grid_tau[0]) 
plt.plot(date_line_rec, line1d_rec)

fig = plt.figure(1)
plt.plot(grid_tau, tran1d)

fig = plt.figure(2)
plt.imshow(tran_check, origin="low", aspect='auto', extent=[grid_vel[0], grid_vel[-1], grid_tau[0], grid_tau[-1]])

plt.show()
