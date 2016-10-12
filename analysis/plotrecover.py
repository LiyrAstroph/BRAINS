##
# \file plotrecover.py
# \brief plot the recovered light curves.
#

import sys
import numpy as np
import matplotlib.pyplot as plt

import plotline2d

def set_cov_Pmat(sigma, tau, alpha, Tcon):
	Pmat = np.zeros((Tcon.shape[0], Tcon.shape[0]))
	xv, yv = np.meshgrid(Tcon, Tcon)
	Pmat = sigma*sigma*np.exp( - pow( np.fabs( xv - yv )/tau, alpha) )
	return Pmat


def reconstruct_con():
	data = np.loadtxt("../data/mrk486_con.txt")
	scale =  np.mean(data[:, 1])
	
	sample  = np.loadtxt("../data/posterior_sample.txt")
	
	pcon = np.loadtxt("../data/pcon.txt")
	
	fig = plt.figure(figsize=(6, 4))
	con = np.zeros(pcon.shape[0])
	for i in range(sample.shape[0]):
		Pmat = set_cov_Pmat(np.exp(sample[i, 1]), np.exp(sample[i, 2]), 1.0, pcon[:, 0])
		Mmat = np.linalg.cholesky(Pmat)
		#Mmat = Mmat.T

		con = np.matmul(Mmat, sample[i, 4:]) + sample[i, 3]
		plt.plot(pcon[:, 0], con * scale, color='grey', lw=0.05)

	plt.plot(pcon[:, 0], pcon[:, 1], color='r', lw=2)
	plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], ls='none', marker='o', color='b')
	plt.xlabel("HJD")
	plt.ylabel("Fcon")
	plt.savefig("rec_con.pdf", bbox_inches='tight')
	
	plt.close()


def reconstruct_line1d():
	con_data = np.loadtxt("../data/mrk486_con.txt")
	con_scale =  np.mean(con_data[:, 1])
	
	line_data = np.loadtxt("../data/mrk486_hb.txt")
	#line_data[:, 1:] /= np.mean(line_data[:, 1])

	
	sample  = np.loadtxt("../data/posterior_sample1d.txt")
	
	pcon = np.loadtxt("../data/pcon.txt")
	pline = np.loadtxt("../data/pline.txt")

	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)

	offset = 9
	con = np.zeros(pcon.shape[0])
	for i in range(sample.shape[0]):
		Pmat = set_cov_Pmat(np.exp(sample[i, 1+offset]), np.exp(sample[i, 2+offset]), 1.0, pcon[:, 0])
		Mmat = np.linalg.cholesky(Pmat)
		#Mmat = Mmat.T

		con = np.matmul(Mmat, sample[i, 4+offset:]) + sample[i, 3+offset]
		ax1.plot(pcon[:, 0], con * con_scale, color='grey')

	ax1.errorbar(con_data[:, 0], con_data[:, 1], yerr=con_data[:, 2], ls='none', marker='o', color='r')
	ax1.plot(pcon[:, 0], pcon[:, 1], lw = 2, color='b')
	ax1.set_ylabel("Fcon")

	ax2.errorbar(line_data[:, 0], line_data[:, 1], yerr = line_data[:, 2], ls='none', marker='o', color='r')
	ax2.plot(pline[:, 0], pline[:, 1]+0.1)
	ax2.set_xlim(ax1.get_xlim())
	ax2.set_xlabel("HJD")
	ax2.set_ylabel("Fline")
	
	fig.savefig("rec_line1d.pdf", bbox_inches='tight')
   
def reconstruct_line2d():
	plotline2d.plotline2d()
  
if __name__ == "__main__":
	if(len(sys.argv) < 2):
		print("No dimension specified.")
		exit(0);

	if int(sys.argv[1]) == 0:
		reconstruct_con()
	elif int(sys.argv[1]) == 1:
		reconstruct_line1d()
	elif int(sys.argv[1]) == 2:
		reconstruct_line2d()
	else:
		print("Incorrect dimension.")




