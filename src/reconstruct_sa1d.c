/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_sa1d.c
 *  \brief reconstruct sa and 1d RM data.
 */

#ifdef SA

#include "brains.h"

void *best_model_sa1d;      /*!< best model */
void *best_model_std_sa1d;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess_sa1d()
{
}

void reconstruct_sa1d()
{
  int i, argc=0;
  char **argv;

  //configure restart of dnest
  argv = malloc(9*sizeof(char *));
  for(i=0; i<9; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restartsa_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartsa_dnest.txt");
  }
  if(parset.flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(parset.flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", parset.temperature);
  }
  if(parset.flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }
  
  //level-dependent sampling
  {
    strcpy(argv[argc++], "-l");
  }
  
  reconstruct_sa1d_init();
  
  dnest_sa1d(argc, argv);

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_sa1d();
  }

  reconstruct_sa1d_end();

  //clear up argv
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sa1d_init()
{
  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSSA1D");

  /* cloud sample related */
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));
}

void reconstruct_sa1d_end()
{
  free(clouds_weight);
  free(clouds_alpha);
  free(clouds_beta);
  free(clouds_tau);
}
#endif