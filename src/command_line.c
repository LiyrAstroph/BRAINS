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

#include "brains.h"
#include "mygetopt.h"

int command_line_options(int argc, char** argv)
{
  int opt, opt_idx, info=EXIT_SUCCESS;
  extern int my_optind, my_opterr, my_optopt;
  extern char *my_optarg;
  extern int my_getopt_long(int argc, char *const *argv, const char *shortopts, 
                            const struct my_option *longopts, int *indexptr);

  /* cope with command options. */
  if(thistask == roottask)
  {
    static struct my_option long_options[] = 
    {
      {"para_name", my_no_argument, 0, 'n'},
      {"version", my_no_argument, 0, 'v'},
      {"help", my_no_argument, 0, 'h'},
      {"restart", my_no_argument, 0, 'r'},
      {"post_proc", my_no_argument, 0, 'p'},
      {"exam_prior", my_no_argument, 0, 'e'},
      {"recalc_info", my_no_argument, 0, 'c'},
      {"temperature", my_required_argument, 0, 't'},
      {"temp", my_required_argument, 0, 't'},
      {"seed", my_required_argument, 0, 's'},
      {"gravity", my_no_argument, 0, 'a'},
      {"load_prior", my_required_argument, 0, 'l'},
      {"para_value", my_required_argument, 0, 'u'},
      {0, 0, 0, 0}
    };

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
    parset.flag_force_run = 1;
    parset.flag_gravity = 0;
    parset.flag_load_prior = 0;
    parset.flag_para_value = 0;
    
    my_opterr = 0; /* surpress getopt error message. */
    my_optind = 0; /* in GNU, optind=0 reset getopt */

    while( (opt = my_getopt_long(argc, argv, "pt:rcs:ehvnfal:u:", long_options, &opt_idx)) != -1)
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
          parset.temperature = atof(my_optarg);
          printf("# Set a temperature %f.\n", parset.temperature);
          if(parset.temperature == 0.0)
          {
            printf("# Incorrect option -t %s.\n", my_optarg);
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
          parset.rng_seed = atoi(my_optarg);
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
        
        case 'l': /* Load parameter prior */
          parset.flag_load_prior = 1;
          strcpy(parset.prior_file, my_optarg);
          printf("# Load parameter prior from %s.\n", parset.prior_file);
          break;

        case 'a': /* use gravity's baseline */
          printf("# Use GRAVITY's 3C 273 baselines.\n");
          parset.flag_gravity = 1;
          break;
        
        case 'u': /* input parameter values */
          parset.flag_para_value = 1;
          strcpy(parset.para_value_file, my_optarg);
          printf("# Load parameter values from %s.\n", parset.para_value_file);
          break;
          
        case '?':
          printf("# Incorrect option -%c %s.\n", my_optopt, my_optarg);
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
        if(argv[my_optind] != NULL) /* parameter file is specified */
          strcpy(parset.param_file, argv[my_optind]); /* copy input parameter file */
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