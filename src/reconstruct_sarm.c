/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */


#ifdef SpecAstro

#include "brains.h"

void reconstruct_sarm()
{
  int i, argc=0;
  char **argv;

  //configure restart of dnest
  argv = malloc(11*sizeof(char *));
  for(i=0; i<11; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restartsarm_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartsarm_dnest.txt");
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
  
  // sample tag 
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "_sarm");

  reconstruct_sarm_init();
  
  sarm_smooth_init(n_vel_sarm_data, vel_sa_data, parset.sa_InstRes);
  
  dnest_sarm(argc, argv);

  sarm_smooth_end();
  reconstruct_sarm_end();

  //clear up argv
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sarm_init()
{

}

void reconstruct_sarm_end()
{

}
#endif