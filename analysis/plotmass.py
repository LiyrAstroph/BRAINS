import numpy as np
import matplotlib.pyplot as plt
import sys

def read_params():
  """
  read parameter file
  """
  fp = open("./src/param")

  tags = ["FileDir", "FlagBLRModel", "FlagDim"]
  params={}

  for line in fp.readlines():
    if line[0] == "#" or line.strip() == "":
      continue

    sp = line.split()
    for tag in tags:
      if tag == sp[0]:
        params[tag] = sp[1].lower()
    
  fp.close()
  
  return params

def plotmass():
  params = read_params()
  
  if int(params["FlagDim"]) < 2:
  	print "not 2d RM analysis. no mass parameter."
  	exit(0)

  sample = np.loadtxt(params["FileDir"]+"/data/posterior_sample_2d.txt", skiprows=1)
  sample_info = np.loadtxt(params["FileDir"]+"/data/posterior_sample_info_2d.txt", skiprows=1)
    
  idx_mass = [8, 8, 8, 8, 11, 10, 15]
  
  idx = idx_mass[int(params["FlagBLRModel"])-1]

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.hist(sample[:, idx]/np.log(10.0)+6.0, bins=20, density=True)
  ax.set_xlabel(r"$\log (M_\bullet/M_\odot)$")
  plt.show()

if __name__ == "__main__":

	if len(sys.argv) < 2:
		print "Please specify a parameter file."
		exit(0)

	plotmass()
