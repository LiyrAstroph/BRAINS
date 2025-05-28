.. _narrowline_label:

************************
Narrow Line Modeling
************************

``brains`` can also include the narrow line to the line profile modeling. 
The narrow line is modeled as a single Gaussian component.
In the parameter file (see :ref:`Parameter File`),  ``FlagNarrowLine`` 
and other related options can be specified to control the behavior of the narrow line modeling.
The options include::

  %========================================================================
  % RM narrow-line component
  % use a gaussian to model the narrow-line component
  % width is the standard deviation of the Gaussian
  % note that here width is the instrinsic width (after correcting for spectral broandening)
  % 
  FlagNarrowLine              0                           % 0, no narrow line; 
                                                          % 1, add fixed narrow line; 
                                                          % 2, add Gaussian priors of the flux for narrow line; 
                                                          % 3, add logarithmic prior of the flux for narrow line

  FluxNarrowLine              1.5                         % flux of narrow line
  FluxNarrowLineErr           1.50                        % flux error of narrow line
  WidthNarrowLine             93.0                        % width km/s
  WidthNarrowLineErr          10.0                         % width error
  ShiftNarrowLine             0.0                        % shift, km/s, with respect to broad line center.  
  ShiftNarrowLineErr          0.0                         % shift error

  FluxNarrowLineLimitLow      0.01                        % when FlagNarrowLine=3
  FluxNarrowLineLimitUpp      1.00                        % when FlagNarrowLine=3

Here, ``FlagNarrowLine`` controls whether the narrow line is included in the model.

- If ``FlagNarrowLine=0``, no narrow line is included.
- If ``FlagNarrowLine=1``, a fixed narrow line is included with the flux and width specified 
  by ``FluxNarrowLine`` and ``WidthNarrowLine``.

- If ``FlagNarrowLine=2``, a Gaussian prior is applied to the flux, width, and shift of the narrow line, 
  The flux prior is specified by ``FluxNarrowLine`` and ``FluxNarrowLineErr``, the width prior 
  is specified by ``WidthNarrowLine`` and ``WidthNarrowLineErr``, and the shift prior is 
  specified by ``ShiftNarrowLine`` and ``ShiftNarrowLineErr``.

- If ``FlagNarrowLine=3``, a logarithmic prior is applied to the flux of the narrow line, 
  with the flux and width specified by ``FluxNarrowLine`` and ``WidthNarrowLine``. 
  The limits for the flux are specified by ``FluxNarrowLineLimitLow`` and ``FluxNarrowLineLimitUpp``.
  The width and shift priors are the same as in case 2, specified by 
  ``WidthNarrowLine``, ``WidthNarrowLineErr``, ``ShiftNarrowLine``, and ``ShiftNarrowLineErr``.

The narrow-line flux refers to the total flux of the narrow line, with the unit consistent with the input line 
profile data. One usually encountered error is that the flux unit is not properly accounted for or the 
flux is not properly specified, leading to a failure in the narrow-line modeling.

For model parameters with Gaussian priors, the parameterized quantities are the deviations from the prior 
means in unit of the standard deviations.

Note that in the code, the flux of the narrow line uses a unit different from the normal unit for the 
convenience of calculations. The code will automatically convert the flux unit and the user does not 
need to worry about it. However, the output Markov chain still uses the internal unit in the code, so
the user should be aware of this when directly inspecting the posterior samples. Specifically,
the conversion relation between the narrow-line flux in the code and the real data is

.. math::

  f_{\rm real} = f_{\rm code}\times \frac{\lambda_0}{C_{\rm unit}} \times \frac{1}{f_{\rm scale}},

where :math:`\lambda_0` is the central wavelength of the line, :math:`f_{\rm scale}` 
is the scaling factor for the flux, which will be output on the screen when running the code, 
and

.. math::
  
  C_{\rm unit} = \frac{c}{\rm VelUnit},

  {\rm VelUnit} = \sqrt{10^6 G M_\odot /{\rm lightday}} / 10^5,

where :math:`c` is the speed of light, :math:`G` is the gravitational constant,
:math:`M_\odot` is the solar mass, and :math:`{\rm lightday}` is the light travel 
distance in one day. Here, VelUit is equivalent to the Keplerian velocity (in unit of km/s) at a distance of 
1 light-day to a :math:`10^6 M_\odot` black hole. 