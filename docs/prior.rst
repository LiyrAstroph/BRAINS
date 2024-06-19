.. _prior_label:

***********************
Adjusting Prior Ranges
***********************

BRAINS sets default prior ranges for parameters in MCMC sampling. In some cases, these default prior 
ranges might be inapproprite.  BRAINS allows to adjust prior ranges as follows.

First, generate the default prior ranges adopted in BRAINS with the command 

.. code-block:: bash

    ./brains param -n 

This will output the prior ranges to a file named `data/para_names_xx.txt`. The content looks like::

   # Par                                 Min        Max Prior  Fix             Val       Mean(Gau)        Std(Gau)
   0 BLR_model_ln(Rblr)            -2.302585   2.504504    2    0  -1.797693e+308    0.000000e+00    0.000000e+00
   1 BLR_model_beta                 0.001000   2.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00
   2 BLR_model_F                    0.000000   1.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00
   3 BLR_model_Inc                  0.000000   1.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00
   4 BLR_model_Opn                  0.000000  90.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00

Here, `Min` and `Max` columns represent the lower and upper limits. For `Prior` column, `1` means Gaussian prior 
and `2` means uniform prior. For `Fix` column, `0` means not fixed and `1` means fixed. When the parameter is fixed,
the `Val` columns represent the fixed value. When the prior is Gaussian, the `Mean(Gau)` and `Std(Gau)` represent
the mean and standard deviation of the Gaussian, respectively. The lines starting with "#" will be neglected.

**Note that do not change the format the file, otherwise, there will be an error when reading in it.**

If the prior range of some parameters needs to change, just adjust the corresponding `Min` and `Max` columns. 
Then save the edited file to a new file, e.g., say, `data/new_prior.txt`. 
Note that for the sake of brevity, one can only keep those lines for the parameters to be adjusted. The rest lines can be removed. 
However, it is still fine to keep all lines. For example, one can edit the 1st and 3rd parameters (counting from 0) as::

   1 BLR_model_beta                 0.001000   2.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00
   3 BLR_model_Inc                  0.000000   1.000000    2    0  -1.797693e+308    0.000000e+00    0.000000e+00

Afterwards, pass this new prior file to BRAINS as 

.. code-block:: bash

    mpiexec -n 6 ./brains param -l data/new_prior.txt

BRAINS will read in the prior ranges and use them for MCMC sampling.