****************
Data
****************

Data Input
==========

The input data in ``brains`` depend on the option ``FlagDim`` specified in the parameter file (see :ref:`Parameter File`). 

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
    
    # time
    wavelength bin 0, flux, error
    wavelength bin 1, flux, error
    ...
    wavelength bin nb, flux, error
  
  where time is the epoch of the current block. Blocks are separated by a blank line.

  .. note::

    * Time should be given in **observed frame**. After reading in the data, 
      the code automatically converts the time into 
      **rest frame** by dividing with a factor :math:`(1+z)`, to account for redshift effects.

    * Wavelength bins should be equally spaced and the unit of wavelength should be the same
      as that of ``LineCenter`` option in :ref:`Parameter File`.

    * Wavelength is converted into velocity as 

      .. math::

        V = \frac{\lambda/(1+z) - \lambda_0}{\lambda_0}, 
      
      where :math:`\lambda_0` is the rest-frame line center and :math:`z` is the redshift.

Data Mask
=========

For 2D line data, ``brains`` can mask those velocity bins with negative 
errors.