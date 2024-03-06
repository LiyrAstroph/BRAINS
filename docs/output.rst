.. _outputs:

********
Outputs
********
BRAINS outputs following files.

- **para_names_*.txt**: parameter names, prior types, and prior ranges.
 
  .. note::
    
    * Some of parameters are in natural logarithmic space.
    * The black hole mass paramter is in a unit of :math:`10^6M_\odot`.
    * The inclination parameter is in a form of :math:`\cos i`, where :math:`i` is the angle 
      between the line of sight and the rotating axis of the BLR.
  

- **posteroir_sample_*.txt**:
  posterior sample of parameters from MCMC sampling. 
  
  The file contains a number of rows and columns. The number of rows represents the effective 
  sample size. The number of columns equal to the number of model parameters, that is to say, 
  different columns correspond to different parameters. The arranging order is the same as 
  that in para_names_*.txt.

- **con_rec.txt**: reconstruction of continuum light curve.

- **line_rec.txt**: reconstruction of 1D line light curve.
  
- **line2d_rec.txt**: reconstruction of 2D line light curve.
  
- **tran_rec.txt**: 1D transfer function.
  
- **tran2d_rec.txt**: 2D tranfer function.

See :ref:`plot_label` for plotting the results.