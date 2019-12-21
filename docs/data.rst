****************
Data
****************

Data Input
==========

The input data in ``brains`` depend on the option ``FlagDim`` specified in the parameter file. 

* ``FlagDim = -2``

  create fully random simulation data; no need to input data.


* ``FlagDim = -1``

  create simulation data according to input continuum and 2D line data.
  ``ContinuumFile`` and ``Line2DFile`` need to be specified.


* ``FlagDim = 0``
  
  only perform continuum reconstruction. 
  ``ContinuumFile`` needs to be specified. 


* ``FlagDim = 1``
  
  perform 1D line analysis.
  ``ContinuumFile`` and ``Line1DFile`` need to be specified.


* ``FlagDim = 2``

  perform 2D line analysis.
  ``ContinuumFile`` and ``Line2DFile`` need to be specified.


Data Format
===========

* Light curve data of continuum and 1D line should 
  contain three columns, which represent time, flux, and 
  flux uncertainty, respectively.


* 2D line data are formated as follows:
  
  The first line is ``# ne nb``, which specify the number ``ne`` of epochs
  and the number ``nb`` of velocity bins.

  Then there are ``ne`` blocks. Each blocks looks like::

    velocity bin 0, time, flux, error
    velocity bin 1, time, flux, error
    ...
    velocity bin nb, time, flux, error
  
  where time is the epoch of the current block. Blocks are separated by a blank line.

  .. note::
    Velocity is in a unit of km/s and velocity bins should be equally spaced.

Data Mask
=========

For 2D line data, ``brains`` can mask those velocity bins with negative 
errors.