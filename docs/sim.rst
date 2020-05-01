**************************
Generating Simulation Data
**************************

The code can generate simulation data according the ``FlagDim`` option in parameter file.

* ``FlagDim = -1``
  
  Generate simulation data on time and wavelength grids same as the input data.


* ``FlagDim = -2``
  
  Generate fully random simulation data.

The output simulation data files are ``sim_con.txt``, ``sim_line.txt``, ``sim_line2d.txt``, ``sim_broadening.txt``,
and ``sim_sa.txt`` in the subdirectory ``data/``.