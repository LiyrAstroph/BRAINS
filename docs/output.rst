.. _outputs:

********
Outputs
********
BRAINS outputs following files.

- **para_names_*.txt**: parameter names, prior types, and prior ranges.
 
  .. note::
    
    * Some of parameters are in natural logarithmic space.
    * The black hole mass paramter is in a unit of** :math:`10^6M_\odot`.
  

- **posteroir_sample_*.txt**:
  posterior sample of parameters from MCMC sampling. 
  
  In each row, the parameters are aranged in an order same as they appear in para_names_*.txt.

- **con_rec.txt**: reconstruction of continuum light curve.

- **line_rec.txt**: reconstruction of 1D line light curve.
  
- **line2d_rec.txt**: reconstruction of 2D line light curve.
  
- **tran_rec.txt**: 1D transfer function.
  
- **tran2d_rec.txt**: 2D tranfer function.