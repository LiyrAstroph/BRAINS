/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
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
  double t0, t1, dt;
  int opt;
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  
  if(thistask == roottask)
  {
    t0 = second();
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

  /* cope with command options. */
  if(thistask == roottask)
  {
    opterr = 0;
    parset.flag_postprc = 0; /* default value */
    parset.temperature = 1.0; /* default value */
    while( (opt = getopt(argc, argv, "pt:")) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocess */
          parset.flag_postprc = 1;
          printf("# MCMC samples available, only do post-process.\n");
          break;
        case 't': /* temperature for postprocess */
          parset.temperature = atof(optarg);
          printf("# Set a temperature %f.\n", parset.temperature);
          if(parset.temperature == 0.0)
          {
            printf("# Incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(parset.temperature < 1.0)
          {
            printf("# Temperature should >= 1.0\n");
            exit(0);
          }
          break;
        case '?':
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;
        default:
          break;
      }
    }
    strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
  }

  begin_run();    /* implementation run */

  end_run();      /* end run */
  
  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
    int ht, mt;
    double st;
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
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
  be downloaded at http://www-unix.mcs.anl.gov/mpi/mpich/. \n Note that 
  any vendor supplied versions exist.

- \b GSL - the GNU Scientific Library (version 2.2.1 or higher).\n
  This open-source package can be downloaded at http://www.gnu.org/software/gsl.\n
  BRAINS uses this library to generate random numbers, perform FFT
  tranformation, and access to CBLAS library. 

- \b LAPACKE - the C-interface of LAPACK (version 3.6.1 or higher).\n
  This open-source package can be downloaded at http://www.netlib.org/lapack/.\n
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

% lines beginning with "%" are regarded as comments and are neglected

<table>
<caption id="parameterfile">Parameter File</caption>
<tr><th>Parameter                   <th>Value                       <th>Note
<tr><td>FileDir                     <td> ./                         <td>
<tr><td>FlagDim                     <td>2                           <td> % 0, only continuum; 1, 1d line; 2, 2d line 
<tr> <td colspan="3">%=============================== \n
% data file
<tr><td>ContinnumFile               <td>data/mrk382_con.txt         <td>% file for contnnum data \n
<tr><td>LineFile                    <td>data/mrk382_hb.txt          <td>% file for line data     \n
<tr><td>Line2DFile                  <td>data/mrk382_hb2d.txt        <td>% file for line 2d data  \n
<tr><td colspan="3">%=============================== \n
% reconstruction
<tr><td>NConRecon                   <td>100                         <td>% number of points for continuum reconstruction \n
<tr><td>ConConstructFileOut         <td>data/pcon.txt               <td>% output filename  \n

<tr><td>NLineRecon                  <td>200                         <td>% number of points for line reconstruction along time axis \n
<tr><td>LineConstructFileOut        <td>data/pline.txt              <td>% output filename \n
<tr><td>TranFileOut                 <td>data/tran.txt               <td>% output filename for 1d transfer function \n

<tr><td>NVelRecon                   <td>100                         <td>% number of points for line reconstruction along velocity axis \n      
<tr><td>Line2DConstructFileOut      <td>data/pline2d.txt            <td>% output filename   \n
<tr><td>Line2DDataConstructFileOut  <td>data/pline2d_data.txt       <td>% output filename for 2d transfer function at points same with data \n
<tr><td>Tran2DFileOut               <td>data/tran2d.txt             <td>% output filename for 2d transfer function \n
<tr><td>Tran2DDataFileOut           <td>data/tran2d_data.txt        <td>% output filename for 2d transfer function \n

<tr><td>NCloudPerCore               <td>20000                       <td>
<tr><td>NVPerCloud                  <td>10                          <td>% number of velocities per cloud \n

<tr><td>NTau                        <td>100                         <td>% number of time-lag points calculated in transfer function \n
<tr><td>TauMinSet                   <td>0.0                         <td>% minimum time lag calculated \n
<tr><td>TauMaxSet                   <td>50.0                        <td>% maximum time lag calculated \n

<tr><td>FlagCloudsOut               <td>0                           <td>% 1, save clouds at the last run; 0, do not save \n
<tr><td>CloudsFileOut               <td>data/clouds.txt             <td>% output filename for clouds  \n

<tr><td>InstRes                     <td>280                         <td>% instrument resolution, in km/s \n

<tr><td colspan="3">%=============================== \n
% set fixed BLR parameters and their fixed values \n
% do not put sapce in the strings \n
% 1: fixed; 0: not fixed; \n
% values are separated by ":" \n

<tr><td>BLRParFix                   <td>0000000000                   <td>
<tr><td>BLRParFixVal                <td>1.38:22.0:34.5               <td>       

</table>
 */

/*! \page optionsfile Options file of BRAINS
  
  The following options constrols the implementation of diffusive nested 
  sampling (dnest). The meanings of these options are as follows. \n

  - Do not change the order of lines. \n
  - Lines beginning with '#' are regarded as comments. \n

 <table>
 <tr><th>Value                       <th>Note
 <tr><td>\b 2     <td># Number of particles \n
 <tr><td>\b 5000  <td># new level interval \n
 <tr><td>\b 5000  <td># save interval \n
 <tr><td>\b 200   <td># threadSteps - how many steps each thread should do independently before communication \n
 <tr><td>\b 120   <td># maximum number of levels \n
 <tr><td>\b 10    <td># Backtracking scale length (lambda in the paper) \n
 <tr><td>\b 100   <td># Strength of effect to force histogram to equal push. 0-10 is best. (beta in the paper) \n
 <tr><td>\b 1000  <td># Maximum number of saves (0 = infinite) \n
 <tr><td>\b data/sample1d.txt            <td># sample file \n
 <tr><td>\b data/sample_info1d.txt       <td># sample_info file \n
 <tr><td>\b data/levels1d.txt            <td># \n
 <tr><td>\b data/sampler_state1d.txt     <td>#  \n
 <tr><td>\b data/posterior_sample1d.txt  <td># \n
 <tr><td>\b data/posterior_sample_info1d.txt  <td># \n
 <tr><td>\b data/limits1d.txt            <td> #
  </table>
 */

/*! \page BRAINS-Makefile  Makefile of BRAINS
 Makefile
 --------

 <table>
 <tr><th> Command                     <th> Note 
 <tr><td>SHELL=/bin/bash              <td> Bash shell
 <tr><td>CC       = mpicc             <td> MPI C Compiler
 <tr><td>OPTIMIZE = -O2 -Wall         <td> Compile option
 <tr><td>#OPTIMIZE += -DDebug         <td> Debug with fixed random seed
 <tr><td>
     #SYSTEM="Darwin"\n
      SYSTEM="Linux" \n
     #SYSTEM="Cluster" \n                 <td> Target system
 </table>

 */