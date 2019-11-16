************************
Broad-line region models
************************

``brains`` already encapsulates several BLR dynamical models detailed below. To specify those models in 
running, edit the option "FlagBLRModel" in the parameter file.

Broad-line regions are generally assumed to be composed of a larg number of 
point-like clouds. These clouds respond to the central ionizing continuum
and emit broad emission lines.

.. figure:: _static/fig_blr_disk.jpg
  :align: center

  Schematic of a disk-like broad-line region.

BLR model 1
===========
Clouds' distribution has a disk-like shape (see the figure) and is axis-symetric.

* Radial distribution: :math:`\Gamma`-distribution

  .. math::
  
    \Gamma(r|\alpha, \theta) = 
    \frac{1}{\Gamma(\alpha)\theta^{\alpha-1}}r^{\alpha-1}\exp\left(-\frac{r}{\theta}\right)

* Dynamics: clouds' orbital angular momentum and energy are randomly assigned following distributions

  .. math::

    E = \left(\frac{1}{1+\exp(-\chi)}\right)E_{\rm min},~~~
    E_{\rm min}=-\frac{GM_\bullet}{r}, ~~~\chi\sim N(0, \lambda^2)\\
    p(L)\sim \exp\left(-\frac{|L|}{\lambda L_{\rm max}}\right),~~~
    |L| < L_{\rm max} = \sqrt{2r^2\left(E+\frac{GM_\bullet}{r}\right)}

* Emissivity: anisotropic prescription

  .. math::
  
    w(\phi) = \frac{1}{2} + \kappa \cos\phi

  where :math:`\phi` is the angle between the observer's line of sight to the central ionizing 
  source and the cloud's line of sight to the central source.

BLR model 2
===========

BLR model 3
===========

BLR model 4
===========

BLR model 5
===========

BLR model 6
===========
This is compatible with Pancoast et al. (2014)'s model.

BLR model 7
===========
This is the shadowed model in Li et al. (2018).

.. figure:: _static/fig_blr_twozone.jpg
  :align: center 

  Schematic of a disk-like broad-line region with two zones.

