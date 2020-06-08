/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file command_line.c
 *  \brief cope with command-line options.
 */

#include <getopt.h>

#include "brains.h"

int command_line_options(int argc, char** argv)
{
  int opt, opt_idx, info=EXIT_SUCCESS;

  /* cope with command options. */
  if(thistask == roottask)
  {
    static struct option long_options[] = 
    {
      {"para_name", no_argument, 0, 'n'},
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"restart", no_argument, 0, 'r'},
      {"post_proc", no_argument, 0, 'p'},
      {"exam_prior", no_argument, 0, 'e'},
      {"recalc_info", no_argument, 0, 'c'},
      {"temperature", required_argument, 0, 't'},
      {"temp", required_argument, 0, 't'},
      {"seed", required_argument, 0, 's'},
      {"gravity", no_argument, 0, 'a'},
      {0, 0, 0, 0}
    };

    opterr = 0; /* reset getopt. */
    optind = 0; /* reset getopt. */
    parset.flag_postprc = 0; /* default value, 0 means postprocessing after runing MCMC sampling. */
    parset.temperature = 1.0; /* default value */
    parset.flag_restart = 0;
    parset.flag_sample_info = 0;
    parset.flag_temp = 0;
    parset.flag_exam_prior= 0;
    parset.flag_rng_seed = 0;
    parset.flag_help = 0;
    parset.flag_end = 0;
    parset.flag_para_name = 0;
    parset.flag_force_run = 0;
    parset.flag_gravity = 0;

    while( (opt = getopt_long(argc, argv, "pt:rcs:ehvnfa", long_options, &opt_idx)) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocessing */
          parset.flag_postprc = 1;
          parset.temperature = 1.0;
          printf("# MCMC samples available, only do post-processing.\n");
          break;
        case 't': /* temperature for postprocessing */
          parset.flag_temp = 1;
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

        case 'r':   /* restart from restored file */
          parset.flag_restart = 1;
          printf("# Restart run.\n");
          break;

        case 'c':  /* calculate sample information */
          printf("# Recalculate the sample info.\n");
          parset.flag_sample_info = 1;
          break;

        case 's':  /* set random number generator seed */
          parset.flag_rng_seed = 1;
          parset.rng_seed = atoi(optarg);
          printf("# Set random seed %d.\n", parset.rng_seed);
          break;

        case 'e':  /* examine the priors assigned */
          printf("# Examine priors.\n");
          parset.flag_exam_prior = 1;
          break;
        
        case 'f':
          parset.flag_force_run = 1;
          break;

        case 'h':  /* print help */
          parset.flag_help = 1;
          print_help();
          break;

        case 'v': /* print version */
          parset.flag_help = 1;
          print_version();
          break;

        case 'n': /* print parameter names */
          printf("# Print parameter name.\n");
          parset.flag_para_name = 1;
          break;

        case 'a': /* use gravity's baseline */
          printf("# Use GRAVITY's 3C 273 baselines.\n");
          parset.flag_gravity = 1;
          break;
          
        case '?':
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          info = EXIT_FAILURE;
          break;

        default:
          break;
      }
    }
    
    if(info != EXIT_FAILURE)
    {
      if(parset.flag_postprc == 1 || parset.flag_sample_info == 1)
        parset.flag_restart = 0;
      
      if(parset.flag_help == 0) /* not only print help */
      {
        if(argv[optind] != NULL) /* parameter file is specified */
          strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
        else
        {
          parset.flag_end = 1;
          fprintf(stderr, "# Error: No parameter file specified!\n");
        }
      }
    }
  }
  
  MPI_Bcast(&info, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  if(info != EXIT_FAILURE)
  {
    /* broadcast parset */
    MPI_Bcast(&parset, sizeof(PARSET), MPI_BYTE, roottask, MPI_COMM_WORLD);
  
    if(parset.flag_end == 1 && parset.flag_help ==0 )
    {
      if(thistask == roottask)
      {
        fprintf(stdout, "Ends incorrectly.\n");
      }
      return EXIT_FAILURE;
    }
  
    return EXIT_SUCCESS;
  }
  else 
  {
    return EXIT_FAILURE;
  }
  
}