******************
User's BLR Model
******************

``brains`` allows user to define their own BLR model. 
At present, there is one way to do so, but we plan to develop a 
more flexible interface in near future.

To define new BLR model, edit the file ``user_blr_model.c`` as follows:

* define the struct ``MyBLRmodel``, and provide the variables 
  ``num_params_MyBLRmodel1d`` and ``num_params_MyBLRmodel2d``.

* define the function ``set_blr_range_mymodel``.

* define the function ``transfun_1d_cloud_sample_mymodel``.

* define the function ``transfun_2d_cloud_sample_mymodel``.

* define the function ``set_par_value_mymodel_sim``.

Then re-compile the code and set ``FlagBLRModel`` to ``0``
in the parameter file.
