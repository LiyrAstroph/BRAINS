/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/*! \file main.c
 *  \brief start of the program
 */
 
/*!
 *  This function initializes the MPI communication packages.
 */ 
int main(int argc, char **argv)
{
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  
  if(thistask == roottask)
  {
    printf("===============BRAINS==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }

  if(argc<2)
  {
    if(thistask == roottask)
    {
      fprintf(stderr, "# Error: No parameter file specified!\n");
      fprintf(stdout, "Ends incorrectly.\n");
    }
    MPI_Finalize();
    return 0;
  }

  strcpy(parset.param_file, argv[1]); /* copy input parameter file */

  begin_run();    /* implementation run */

  end_run();      /* end run */
  
  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
    printf("Ends successfully.\n");
    printf("===============BRAINS==================\n");
  }
  return 0;
}

/* ----------------------------------------------------------------------
   The rest of this file contains documentation for compiling and
   running the code, in a format appropriate for 'doxygen'.
   ----------------------------------------------------------------------
 */

/*! \mainpage Reference documentation of BRAINS

\author Yan-Rong Li\n
        liyanrong@ihep.ac.cn\n
        Institute of High Energy Physics

\section introduction Introduction

\b BRAINS is an abbreviation for  <b>B</b>LR <b>R</b>everberation-mapping 
   <b>A</b>nalysis <b>I</b>ntegrated with <b>N</b>ested <b>S</b>ampling.

\section install Compilation 

BRAINS needs the following non-standard libraries for compilation:

- \b MPI - the Message Passing Interface (version 1.0 or higher).\n
  The widely used MPI implementation MPICH is recommended, which can 
  be downloaded at http://www-unix.mcs.anl.gov/mpi/mpich/. Note that 
  any vendor supplied versions exist.

- \b GSL - the GNU Scientific Library (version 2.2.1 or higher).\n
  This open-source package can be downloaded at http://www.gnu.org/software/gsl.
  BRAINS uses this library to generate random numbers, perform FFT
  tranformation, and access to CBLAS library. 

- \b LAPACKE - the C-interface of LAPACK (version 3.6.1 or higher).\n
  This open-source package can be downloaded at http://www.netlib.org/lapack/.
  BRAINS uses this library to perform linear algebra computation. 
  Note that to LAPACKE is contained in LAPACK source. One needs to switch on
  LAPACKE when compiling LAPACK.

Note that in Linux system, there are package managers that can install the above 
libraries convienently. If so, use them. In this case, the libraries usually are 
installed in standard environment path. Otherwise, any of the above libraries is 
not installed in standard locations on your system, the \ref BRAINS-Makefile 
"Makefile" provided with the code may need slight adjustments.

\section howtorun Running the code

To run the code in a terminal, type:\n
\code{.sh} 
mpiexec -n np ./brains param
\endcode
where np is the number of cores, param is the parameter file.
In order to start the code, a \ref parameterfile "parameterfile" 
and a \ref optionsfile "optionfile" needs to be specified.
 */


/*! \page parameterfile Parameterfile of BRAINS

FileDir                     ./    \n

FlagDim                     2                           % 0, only continuum; 1, 1d line; 2, 2d line \n
%=============================== \n
% data file
ContinnumFile               data/mrk382_con.txt         % file for contnnum data \n
LineFile                    data/mrk382_hb.txt          % file for line data     \n
Line2DFile                  data/mrk382_hb2d.txt        % file for line 2d data  \n
%=============================== \n
% reconstruction
NConRecon                   100                         % number of points for continuum reconstruction \n
ConConstructFileOut         data/pcon.txt               % output filename  \n

NLineRecon                  200                         % number of points for line reconstruction along time axis \n
LineConstructFileOut        data/pline.txt              % output filename \n
TranFileOut                 data/tran.txt               % output filename for 1d transfer function \n

NVelRecon                   100                         % number of points for line reconstruction along velocity axis \n      
Line2DConstructFileOut      data/pline2d.txt            % output filename   \n
Line2DDataConstructFileOut  data/pline2d_data.txt       % output filename for 2d transfer function at points same with data \n
Tran2DFileOut               data/tran2d.txt             % output filename for 2d transfer function \n

NCloudPerCore               20000
NVPerCloud                  10                          % number of velocities per cloud \n

NTau                        100                         % number of time-lag points calculated in transfer function \n
TauMinSet                   0.0                         % minimum time lag calculated \n
TauMaxSet                   50.0                        % maximum time lag calculated \n

FlagCloudsOut               0                           % 1, save clouds at the last run; 0, do not save \n
CloudsFileOut               data/clouds.txt             % output filename for clouds  \n

InstRes                     280                         % instrument resolution, in km/s \n

%=============================== \n
% set fixed BLR parameters and their fixed values \n
% do not put sapce in the strings \n
% 1: fixed; 0: not fixed; \n
% values are separated by ":" \n

BLRParFix                   0000000000 \n
BLRParFixVal                1.38:22.0:34.5   \n         

 */

/*! \page optionsfile Options file of BRAINS
  
  The following options constrols the implementation of diffusive nested 
  sampling (dnest). The meanings of these options are as follows, 

  \b 2     # Number of particles \n
  \b 5000  # new level interval \n
  \b 5000  # save interval \n
  \b 200   # threadSteps - how many steps each thread should do independently before communication \n
  \b 120   # maximum number of levels \n
  \b 10    # Backtracking scale length (lambda in the paper) \n
  \b 100   # Strength of effect to force histogram to equal push. 0-10 is best. (beta in the paper) \n
  \b 1000  # Maximum number of saves (0 = infinite) \n
  \b data/sample1d.txt            # sample file \n
  \b data/sample_info1d.txt       # sample_info file \n
  \b data/levels1d.txt            # \n
  \b data/sampler_state1d.txt     #  \n
  \b data/posterior_sample1d.txt  # \n
  \b data/limits1d.txt            #
 */

/*! \page BRAINS-Makefile  Makefile of BRAINS
 
 */