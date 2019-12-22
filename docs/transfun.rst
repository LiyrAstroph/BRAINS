****************************
Analytical Transfer Function
****************************
``brains`` also allows users to directly define analytical transfer functions (without BLR dynamical modeling).


To define analytical transfer functions, edit the source code file ``user_transfun.c`` as follows:

* assign the variables ``num_params_MyTransfun1d`` and ``num_params_MyTransfun2d``.

* define the function ``set_blr_range_mytransfun``, which sets the prior range of model parameters for the transfer function.

* define the function ``transfun_1d_cal_mytransfun``, which assigns the transfer function arrays ``TransTan``
  and ``Trans1D``. This function is used for 1D modeling.

* define the function ``transfun_2d_cal_mytransfun``, which assigns the transfer function arrays
  ``TransTau`` and transv. This function is used for 2D modeling.

* define the function ``set_par_value_mytransfun_sim``, which inputs parameter values for generating mock data.

Then re-compile the code and set ``FlagBLRModel`` to ``-1``
in the parameter file.