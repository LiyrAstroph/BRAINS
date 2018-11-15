# BRAINS
Bayesian Reverberation-mapping Analysis Integrated with Nested Sampling

A package for dynamically modeling broad-line regions and analyzing reverberation mapping data, and measuring the black hole mass.

reference: Li, Y.-R., Songshen, Y.-Y., Qiu, J., et al. 2018, ApJ in press

# Compiling
To compile ``RECON``, the following third-party packages are required:
- **MPICH**, an MPI implementation (version 1.0 or higher), available at http://www-unix.mcs.anl.gov/mpi/mpich/
- **FFTW**, a fast Fourier transform library (version 3.0 or higher), available at http://www.fftw.org/
- **GSL**, the GNU Scientific Library (version 2.2.1 or higher), available at http://www.gnu.org/software/gsl
- **DNest**, the diffusive nested sampling library, developed by Yan-Rong Li, which is a C version of the DNest code by Brendon Brewer (https://github.com/eggplantbren/DNest3), available at https://github.com/LiyrAstroph/DNest_C

After installing the above packages, edit the corresponding library paths in ``Makefile``, and then use ``make`` to create ``brains``.

# Running

```bash
mpiexec -n np ./brains src/param
```

where ``np`` is the number of cores and ``param`` is the parameter file, which specifies configurations for ``brains``.

An exemplary reveberation mapping dataset is provided in the subdirectory ``data/``, containing four files:
```bash
sim_con.txt                 # continuum light curve
sim_hb.txt                  # 1d broad-line flux light curve   
sim_hb2d.txt                # 2d broad-line time series
sim_broadening.txt          # spectral broadening data
```

One can try to run the above command to test ``brains`` with the provided dataset.

# MCMC Samping
The output Markov chain is stored in ``data/posterior_sample.txt`` for continuum reconstuction, ``data/posterior_sample1d.txt`` for 1d reverberation mapping analysis, and `data/posterior_sample2d.txt`` for 2d reverberation mapping analysis.

One need to tune the corresponding option files ``OPTIONSCON``, ``OPTIONS1d``, and ``OPTIONS2D`` accordingly, which specify configurations for nested sampling.

# An Exemplary Test
Application to a mock reverberation mapping dataset, see Li, Y.-R., Songshen, Y.-Y., Qiu, J., et al. 2018, ApJ in press.
![Application to a mock reverberation mapping dataset](https://github.com/liyropt/MyGithubPic/blob/master/fig_sim_brains.jpg)

**Contact the author Li, Yan-Rong (liyanrong@mail.ihep.ac.cn) for any problem.**

**A more detailed usage guideline is coming.**

