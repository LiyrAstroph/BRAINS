import numpy as np
import matplotlib.pyplot as plt

def set_cov_Pmat(sigma, tau, alpha, Tcon):
	Pmat = np.zeros((Tcon.shape[0], Tcon.shape[0]))
	xv, yv = np.meshgrid(Tcon, Tcon)
	Pmat = sigma*sigma*np.exp( - pow( np.fabs( xv - yv )/tau, alpha) )
	return Pmat


def reconstruct_con():
	data = np.loadtxt("../data/arp151_con.txt")
	data[:, 1:] /= np.mean(data[:, 1])
	
	sample  = np.loadtxt("../data/posterior_sample.txt")
	
	pcon = np.loadtxt("../data/pcon.txt")
	
	con = np.zeros(pcon.shape[0])
	for i in range(sample.shape[0]):
		Pmat = set_cov_Pmat(np.exp(sample[i, 0]), np.exp(sample[i, 1]), 1.0, pcon[:, 0])
		Mmat = np.linalg.cholesky(Pmat)
		#Mmat = Mmat.T

		con = np.matmul(Mmat, sample[i, 3:]) + sample[i, 2]
		plt.plot(pcon[:, 0], con, color='grey')

	plt.errorbar(data[:, 0], data[:, 1], yerr=data[:, 2], ls='none', marker='o', color='b')
	
	plt.show()
	
if __name__ == "__main__":
	reconstruct_con()




