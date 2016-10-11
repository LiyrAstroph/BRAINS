import sys
import numpy as np
import matplotlib.pyplot as plt

def set_cov_Pmat(sigma, tau, alpha, Tcon):
	Pmat = np.zeros((Tcon.shape[0], Tcon.shape[0]))
	xv, yv = np.meshgrid(Tcon, Tcon)
	Pmat = sigma*sigma*np.exp( - pow( np.fabs( xv - yv )/tau, alpha) )
	return Pmat


def reconstruct_con():
	data = np.loadtxt("../data/mcg06_con.txt")
	scale =  np.mean(data[:, 1])
	data[:, 1:] /= scale
	
	sample  = np.loadtxt("../data/posterior_sample.txt")
	
	pcon = np.loadtxt("../data/pcon.txt")
	
	fig = plt.figure(figsize=(6, 4))
	con = np.zeros(pcon.shape[0])
	for i in range(sample.shape[0]):
		Pmat = set_cov_Pmat(np.exp(sample[i, 1]), np.exp(sample[i, 2]), 1.0, pcon[:, 0])
		Mmat = np.linalg.cholesky(Pmat)
		#Mmat = Mmat.T

		con = np.matmul(Mmat, sample[i, 4:]) + sample[i, 3]
		plt.plot(pcon[:, 0], con, color='grey', lw=0.05)

	plt.plot(pcon[:, 0], pcon[:, 1]/scale, color='r', lw=2)
	plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], ls='none', marker='o', color='b')
	
	plt.savefig("rec_con.pdf", bbox_inches='tight')
	
	plt.close()


def reconstruct_line1d():
	con_data = np.loadtxt("../data/mcg06_con.txt")
	con_scale = np.mean(con_data[:, 1])
	con_data[:, 1:] /= con_scale
	
	line_data = np.loadtxt("../data/mcg06_hb.txt")
	#line_data[:, 1:] /= np.mean(line_data[:, 1])

	
	sample  = np.loadtxt("../data/posterior_sample1d.txt")
	
	pcon = np.loadtxt("../data/pcon.txt")
	pline = np.loadtxt("../data/pline.txt")

	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)

	offset = 8
	con = np.zeros(pcon.shape[0])
	for i in range(sample.shape[0]):
		Pmat = set_cov_Pmat(np.exp(sample[i, 0+9]), np.exp(sample[i, 1+9]), 1.0, pcon[:, 0])
		Mmat = np.linalg.cholesky(Pmat)
		#Mmat = Mmat.T

		con = np.matmul(Mmat, sample[i, 3+9:]) + sample[i, 2+9]
		ax1.plot(pcon[:, 0], con, color='grey')

	ax1.errorbar(con_data[:, 0], con_data[:, 1], yerr=con_data[:, 2], ls='none', marker='o', color='r')
	ax1.plot(pcon[:, 0], pcon[:, 1]/con_scale, lw = 2, color='b')
	
	ax2.errorbar(line_data[:, 0], line_data[:, 1], yerr = line_data[:, 2], ls='none', marker='o', color='r')
	ax2.plot(pline[:, 0], pline[:, 1]+0.1)
	ax2.set_xlim(ax1.get_xlim())
	
	fig.savefig("rec_line1d.pdf", bbox_inches='tight')
   
def reconstruct_line2d():
 	con_data = np.loadtxt("../data/mrk493_con.txt")
	con_scale = np.mean(con_data[:, 1])

	fp = open("../data/posterior_sample2d.txt", "r")
	line = fp.readline()
	st = line.split()
	num_ps = int(st[1])
	fp.close()
	
	pcon = np.loadtxt("../data/pcon.txt")
	ncon = pcon.shape[0]
	
	con_rec = np.zeros((ncon, num_ps))
	fp = open("../data/con_rec.txt", "r")
	for i in range(num_ps):
	  for j in range(ncon):
	    line = fp.readline()
	    print(line)
	    tmp, con_rec[j , i] = line.split()
	  
	  fp.readline()
	
	fp=open("mrk493_hb2d.txt", "r")
	line=fp.readline()
        text=line.split()
	nt=int(text[1])
	nv=int(text[2])

	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax2 = fig.add_subplot(212)

	ax1.plot(pcon[:, 0], con_rec)
	ax1.errorbar(con_data[:, 0], con_data[:, 1], yerr=con_data[:, 2], ls='none', marker='o')

	
	fig.savefig("rec_line2d.pdf", bbox_inches='tight')

  
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




