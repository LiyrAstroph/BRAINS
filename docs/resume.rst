.. _resume_label:

************************
Resuming from Last Run
************************

``brains`` supports resuming from last run. The adopted DNest sampling code ``CDNest`` 
outputs a number of sampling snapshots with a name as ``restart_dnest.txt_xxx`` for continuum reconstruction,
``restart1d_dnest.txt_xxx`` for 1D and ``restart2d_dnest.txt_xxx`` for 2D RM analysis. Here, 
`xxx` mean the steps when the snapshot is saved. Such snapshots are saved every one fifth of ``MaxNumberSaves``.

For 2D RM analysis, rename the above files as 

.. code:: bash

  cp restart2d_dnest.txt_xxx restart2d_dnest.txt 
  #namely, just remove "_xxx" in the file name

Then change the step numbers in the option files (e.g., param/OPTIONS2D) and resume from last run by adding 
an option "-r" as

.. code:: bash 
  
  mpiexec -n np ./brains para/param -r 

where "np" is the number of cores. It must be the same as the number of cores used in last run.