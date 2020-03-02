**************************
Generating Simulation Data
**************************

The code can generate simulation data according the ``FlagDim`` option in parameter file.

* ``FlagDim = -1``
  
  Generate simulation data on time and wavelength grids same as the input data.


* ``FlagDim = -2``
  
  Generate fully random simulation data.

The output simulation data files are ``sim_con.txt``, ``sim_hb.txt``, and ``sim_hb2d.txt``
in the subdirectory ``data/``.