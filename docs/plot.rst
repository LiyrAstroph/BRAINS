.. _plot_label:

**********************
Plotting Interface 
**********************

BRAINS provides a Python plotting interface (named `bbackend`) to visualize the results. The source code is placed in the subdirectory `analysis/`. 
The interface is independent from the BRAINS and can be installed separately.

Change the path into `analysis` and install the plotting interface as follows.

.. code-block:: bash

    python setup.py install --user


Then open up a Python environment and import the interface as 

.. code-block:: python

    from bbackend import bplotlib 

Load the parameter file (e.g., `../param/param`) as 

.. code-block:: python

    pb = bplotlib("../param/param")

Then plot the results as follows. An examplary script `example.py` is also provided to guide the useage.

.. code-block:: python

    flagdim = int(pb.param["flagdim"])

    #===============================================
    # continuum results, flagdim == 0
    if flagdim == 0:
    # plot continuum DRW parameters
    pb.plot_drw_parameters(doshow=False)
    
    # plot CDNest diagnoistics
    temperature = 1
    pb.postprocess(temperature, doshow=False)

    pb.plot_results_con()

    #===============================================
    # RM 1D results, flagdim == 1
    if flagdim == 1:
    # plot continuum DRW parameters
    pb.plot_drw_parameters(doshow=False)
    
    # plot CDNest diagnoistics
    temperature = 1
    pb.postprocess(temperature, doshow=False)

    pb.plot_results_1d(doshow=False)

    # plot 1d transfer function
    # pb.plot_tran1d(tau_range=[a, b])
    pb.plot_tran1d(doshow=False)

    # plot histograms of BLR model parameters
    pb.plot_blrmodel_para_hist(doshow=False)

    # plot clouds' distribution
    pb.plot_clouds("../data/clouds.txt", doshow=False)
    #pb.plot_clouds("../data/clouds.txt", range=[-10, 10], objname="target", format="jpg", velocity=False, doshow=False)

    # plot clouds' distribution viewd from line of sight
    pb.plot_clouds_los("../data/clouds.txt", doshow=False)
    #pb.plot_clouds_los("../data/clouds.txt", range=[-10, 10], objname="target", format="jpg", velocity=False, doshow=False)


    #===============================================
    # RM 2D results, flagdim == 2
    # 
    if flagdim == 2:
    # print parameter names
    pb.print_blrmodel_para_names()
    
    # plot continuum DRW parameters
    pb.plot_drw_parameters(doshow=False)
    
    # plot CDNest diagnoistics
    temperature = 1
    pb.postprocess(temperature, doshow=False)
    
    # get continuum data 
    con_data = pb.get_con_data()
    
    # get line 2D data 
    line2d_data = pb.get_line2d_data()
    
    # pb.plot_results_2d_style2018(doshow=False)
    pb.plot_results_2d_style2022(doshow=False)
    
    # plot 2d tranfer function with the maximum prob
    # pb.plot_tran2d(tau_range=[a, b])
    # pb.plot_tran2d(vel_range=[a, b])
    # pb.plot_tran2d(tau_range=[a, b], vel_range=[c, d], doshow=False)
    pb.plot_tran2d(doshow=False)
    
    # plot 1d transfer function
    # pb.plot_tran1d(tau_range=[a, b], doshow=False)
    pb.plot_tran1d(doshow=False)
    
    # plot histograms of BLR model parameters
    pb.plot_blrmodel_para_hist(doshow=False)

    # plot clouds' distribution
    pb.plot_clouds("../data/clouds.txt", doshow=False)
    #pb.plot_clouds("../data/clouds.txt", range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)
    
    
    # plot clouds' distribution viewd from line of sight
    pb.plot_clouds_los("../data/clouds.txt", doshow=False)
    #pb.plot_clouds_los("../data/clouds.txt", range=[-10, 10], objname="target", format="jpg", velocity=True, doshow=False)
