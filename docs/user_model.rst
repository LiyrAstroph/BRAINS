******************
User's BLR Model
******************

``brains`` allows user to define their own BLR model. 
At present, there is one way to do so, but we plan to develop a 
more flexible interface in near future.

To define new BLR model, edit the file ``user_blr_model.c`` as follows:

* define the struct ``MyBLRmodel``, and provide the variables 
  ``num_params_MyBLRmodel1d`` and ``num_params_MyBLRmodel2d``.

* define the function ``set_blr_range_mymodel``, which sets the prior range of model parameters.

* define the function ``transfun_1d_cloud_sample_mymodel``, which generates clouds' time lag and weight, 
  stored in arrays "clouds_tau" and "clouds_weight", respectively. This function is used for 1D modeling.

* define the function ``transfun_2d_cloud_sample_mymodel``, which generates clouds' time lag, velocity, and weight, 
  stored in arrays "clouds_tau", "clouds_vel", and "clouds_weight", respectively. This function is used for 2D modeling.

* define the function ``set_par_value_mymodel_sim``, which inputs parameter values for generating mock data.

Then re-compile the code and set ``FlagBLRModel`` to ``0``
in the parameter file.
