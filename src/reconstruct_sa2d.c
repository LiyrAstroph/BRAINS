/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_sa2d.c
 *  \brief reconstruct sa and 2d RM data.
 */

#ifdef SA

#include "brains.h"

void *best_model_sa2d;      /*!< best model */
void *best_model_std_sa2d;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess_sa2d()
{
}

void reconstruct_sa2d()
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
  
  reconstruct_sa2d_init();
  
  dnest_sa1d(argc, argv);

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_sa2d();
  }

  reconstruct_sa2d_end();

  //clear up argv
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sa2d_init()
{
  int i, j;
  double dT, Tspan;

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSSA2D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  
  /* RM setup */
  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  Tcon_min = Tcon_data[0] - time_back_set;
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 10.0);
  Tcon_max = fmax(Tcon_max, Tline_data[n_line_data -1]);  /* The time span should cover that of the emission line data */
  
  if(thistask == roottask)
    printf("Tcon_min_max: %f %f\n", Tcon_min - Tcon_data[0], Tcon_max - Tcon_data[n_con_data-1]);

  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  /* set Larr_rec */
  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }

  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  
  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  if(parset.time_back <= 0.0)
    Tline_min = fmax(Tline_min, Tcon_min + Tspan/2.0 - 10.0);

  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  Tline_max = fmin(Tline_max, Tcon_max);  /* The time span should be smaller than that of the continuum */

  if(thistask == roottask)
    printf("Tline_min_max: %f %f\n", Tline_min - Tline_data[0], Tline_max - Tline_data[n_line_data -1]);

  dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }

  double vel_max_set = Vline_data_ext[n_vel_data_ext -1], vel_min_set = Vline_data_ext[0];
  double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

  for(i=0; i<parset.n_vel_recon; i++)
  {
  	TransV[i] = vel_min_set + dVel*i;
  }

  Fcon_particles = malloc(parset.num_particles * sizeof(double *));
  Fcon_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fcon_particles[i] = malloc(parset.n_con_recon * sizeof(double));
    Fcon_particles_perturb[i] = malloc(parset.n_con_recon * sizeof(double));
  }

  TransTau_particles = malloc(parset.num_particles * sizeof(double *));
  TransTau_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    TransTau_particles[i] = malloc(parset.n_tau * sizeof(double));
    TransTau_particles_perturb[i] = malloc(parset.n_tau * sizeof(double));
  }

  Trans2D_at_veldata_particles = malloc(parset.num_particles * sizeof(double *));
  Trans2D_at_veldata_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans2D_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_data_ext * sizeof(double));
    Trans2D_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_data_ext * sizeof(double));
  }

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_line_data * n_vel_data_ext * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_line_data * n_vel_data_ext * sizeof(double));
  }

  con_q_particles = malloc(parset.num_particles * sizeof(double *));
  con_q_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    con_q_particles[i] = malloc(nq * sizeof(double));
    con_q_particles_perturb[i] = malloc(nq* sizeof(double));
  }
  
  /* SA setup */
  phase_sa_particles = malloc(parset.num_particles * sizeof(double *));
  phase_sa_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    phase_sa_particles[i] = malloc(n_base_sa_data * n_vel_sa_data * sizeof(double));
    phase_sa_particles_perturb[i] = malloc(n_base_sa_data * n_vel_sa_data * sizeof(double));
  }

  Fline_sa_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_sa_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_sa_particles[i] = malloc(n_vel_sa_data * sizeof(double));
    Fline_sa_particles_perturb[i] = malloc(n_vel_sa_data * sizeof(double));
  }

  prob_sa_particles = malloc(parset.num_particles * sizeof(double));
  prob_sa_particles_perturb = malloc(parset.num_particles * sizeof(double));

  /* cloud sample related */
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));

  if(parset.flag_save_clouds && thistask == roottask)
  {
    if(parset.n_cloud_per_task <= 10000)
      icr_cloud_save = 1;
    else
      icr_cloud_save = parset.n_cloud_per_task/10000;

    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }
}

void reconstruct_sa2d_end()
{
  int i;

  /* RM setup */
  free(Tline);
  free(Fline2d);
  free(TransV);
  free(Trans2D);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
  }
  free(Fcon_particles);
  free(Fcon_particles_perturb);
  
  free(Fline_at_data);
  
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans2D_at_veldata_particles[i]);
    free(Trans2D_at_veldata_particles_perturb[i]);
    free(TransTau_particles[i]);
    free(TransTau_particles_perturb[i]);
    free(Fline_at_data_particles[i]);
    free(Fline_at_data_particles_perturb[i]);

    free(con_q_particles[i]);
    free(con_q_particles_perturb[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(Trans2D_at_veldata_particles_perturb);
  free(TransTau_particles);
  free(TransTau_particles_perturb);
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);

  free(con_q_particles);
  free(con_q_particles_perturb);

  /* SA setup */
  for(i=0; i<parset.num_particles; i++)
  {
    free(phase_sa_particles[i]);
    free(phase_sa_particles_perturb[i]);

    free(Fline_sa_particles[i]);
    free(Fline_sa_particles_perturb[i]);
  }
  free(phase_sa_particles);
  free(phase_sa_particles_perturb);

  free(Fline_sa_particles);
  free(Fline_sa_particles_perturb);

  free(prob_sa_particles);
  free(prob_sa_particles_perturb);

  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
    free(par_prior_gaussian[i]);
  }
  free(par_range_model);
  free(par_prior_gaussian);
  free(par_prior_model);

  free(par_fix);
  free(par_fix_val);
  free(best_model_sa2d);
  free(best_model_std_sa2d);

  free(clouds_weight);
  free(clouds_alpha);
  free(clouds_beta);
  free(clouds_tau);

  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  if(thistask == roottask)
  {
    printf("Ends reconstruct_line2d.\n");
  }
}
#endif