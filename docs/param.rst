**************
Parameter File
**************
  
``brains`` reads a parameter file to configure its running. A typical parameter file looks like::
  
  % parameter file
  % lines beginning with '%' are regarded as comments and are neglected
  % 
  
  FileDir                     /home/liyropt/Projects/GIT/BRAINS      % directory for where all the output files are stored
  
  FlagDim                     2                           % -2, generate fully random lc; -1, mock lc based on input data; 
                                                          % 0, only continuum modeling; 1, 1d line RM; 2, 2d line RM
  
  Redshift                    0.0449                         % redshift of the object

  FlagBLRModel                6                           % 1, 2, 3, 4, 5, 6, 7, 8; BLR model type
                                                          % 5 is double power-law model; 6 is Pancoast's model; 7 is two-zone model
                                                          % 0 is user defined
                                                          % -1 is user defined, analytical transfer function
  LineCenter                  4861.0                      %  rest frame, A
  FlagLineCenter              0                           % -1, epoch-dependent offset; 0, not included; 1, uniform offset
  LineCenterErr               50.0                        % km/s
  %========================================================================
  % data file
  
  ContinuumFile               data/mrk142_con.txt          % file for continuum data
  LineFile                    data/mrk142_hb.txt           % file for line data
  Line2DFile                  data/mrk142_hb2d.txt         % file for line 2d data
  
  %========================================================================
  % reconstruction
  NConRecon                   200                         % number of points for continuum reconstruction
  FlagTrend                   0                           % 0, mean; 1, first-order trend
  FlagTrendDiff               0                           % 0, no difference; 1 or larger, add difference in the long-term trends between continuum and line
                                                          % the different trend is modelled by a polynomical with the order set by the value of FlagTrendDiff.
  
  ConConstructFileOut         data/pcon.txt               % output filename for continuum reconstruction
  
  FlagFixVar                  0                           % 0, not fixed; 1, fix the parameters of variation from continuum data.
  NLineRecon                  100                         % number of points for line reconstruction along time axis
  LineConstructFileOut        data/pline.txt              % output filename for 1d line reconstruction
  TranFileOut                 data/tran.txt               % output filename for 1d transfer function
  
  NVelRecon                   42                          % number of points for line reconstruction along velocity axis
  Line2DConstructFileOut      data/pline2d.txt            % output filename for 2d line reconstrction
  Line2DDataConstructFileOut  data/pline2d_data.txt       % output filename for 2d transfer function at points same with data
  Tran2DFileOut               data/tran2d.txt             % output filename for 2d transfer function
  Tran2DDataFileOut           data/tran2d_data.txt        % output filename for 2d transfer function at velocity points same with data
  
  NCloudPerCore               5000                        % number of clouds per task
  NVPerCloud                  5                           % number of velocities per cloud
  
  NTau                        200                         % number of time-lag points calculated in transfer function
  RCloudMax                   -1                          % outer edge of broad-line region in a unit of light-day; -1, set automatically 
  TimeBack                    -1                          % time prior to the start of continuum in reconstruction; -1, set automatically
  
  FlagCloudsOut               1                           % 1, save clouds at the last run; 0, do not save
  CloudsFileOut               data/clouds.txt             % output filename for clouds 
  
  %========================================================================
  %
  FlagCloudsForceUpdate       0                           % default 0; 
                                                          % 0, only update when BLR parameters are updated 
                                                          % 1, update every MCMC perturb (continuum + BLR)
  
  FlagConSysErr               0                           % 0, not include systematic error of contininuum; 1, include
  FlagLineSysErr              0                           % 0, not include systematic error of line; 1, include
  
  FlagNonlinear               1                           % 0, linear response; 1, non-linear response
  %========================================================================
  % spectral broadening
  
  FlagInstRes                 1                           % 0, uniform prior, epoch-independent
                                                          % 1, uniform prior, but epoch-dependent parameterization
                                                          % 2, epoch-dependent parameterization, prior stored in "InstResFile"
  
  InstRes                     220                           % instrument broadening (modeled by a Gaussian), in km/s, for FlagInstRes=0, or 1
                                                            
  InstResErr                  50.0                         % instrument broadening error, in km/s, for FlagInstRes=0, or 1
  
  InstResFile                 data/sim_broaden.txt       % file for storing epoch-dependent instrument broadening
                                                            % two columns: broadening width and error (km/s), in the order of time as the 2d line data
  
  %========================================================================
  % narrow-line component
  % use a gaussian to model the narrow-line component
  
  FlagNarrowLine              0                           % 0, no narrow line; 
                                                          % 1, add fixed narrow line; 
                                                          % 2, add Gaussian priors of the flux for narrow line; 
                                                          % 3, add logarithmic prior of the flux for narrow line
  
  FluxNarrowLine              1.5                         % flux of narrow line
  FluxNarrowLineErr           0.50                        % flux error of narrow line
  WidthNarrowLine             93.0                        % width km/s
  WidthNarrowLineErr          10.0                         % width error
  ShiftNarrowLine             0.0                        % shift, km/s, with respect to broad line center.  
  ShiftNarrowLineErr          0.0                         % shift error  
  
  %========================================================================
  % set fixed BLR parameters and their fixed values
  % do not put sapce in the strings
  % 1: fixed; 0: not fixed;
  % values are separated by ":"
  
  BLRParFix                   0000000000
  BLRParFixVal                2.0:1.0

.. note::
  In the subdirectory ``example/``, some examples of parameter file are provided. Users can choose appropriate 
  parameter files with their purposes.