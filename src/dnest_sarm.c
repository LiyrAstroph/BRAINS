/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifdef SpecAstro

#include "brains.h"

DNestFptrSet *fptrset_sarm;

int dnest_sarm(int argc, char **argv)
{
  int i;
  double logz_sarm;
  char dnest_data_dir[BRAINS_MAX_STR_LENGTH];

  set_sa_blr_model();

  /* RM */
  num_params_blr = 0;
  num_params_rm = parset.n_con_recon + num_params_var;;

  /* SARM */
  num_params_sa_blr = num_params_sa_blr_model + num_params_sa_extpar;
  num_params_sa = num_params_sa_blr;
  
  /* total */
  num_params_blr_tot = num_params_blr + num_params_sa_blr;
  num_params = num_params_sa + num_params_rm;

  idx_resp = num_params_blr_tot + num_params_drw + num_params_trend;
  idx_difftrend = idx_resp + num_params_resp;

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  par_range_model = malloc( num_params * sizeof(double *));
  par_prior_gaussian = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
  {
    par_range_model[i] = malloc(2*sizeof(double));
    par_prior_gaussian[i] = malloc(2*sizeof(double));
  }
  par_prior_model = malloc( num_params * sizeof(int));

  fptrset_sarm = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_sarm->from_prior = from_prior_sarm;
  fptrset_sarm->print_particle = print_particle_sarm;
  fptrset_sarm->restart_action = restart_action_sa;
  fptrset_sarm->accept_action = accept_action_sarm;
  fptrset_sarm->kill_action = kill_action_sarm;
  fptrset_sarm->perturb = perturb_sarm;
  fptrset_sarm->read_particle = read_particle_sarm;

  fptrset_sarm->log_likelihoods_cal_initial = log_likelihoods_cal_initial_sarm;
  fptrset_sarm->log_likelihoods_cal_restart = log_likelihoods_cal_restart_sarm;
  fptrset_sarm->log_likelihoods_cal = log_likelihoods_cal_sarm;

  set_par_range_sarm();
  set_par_fix_sa_blrmodel();

  for(i=num_params_blr + num_params_sa_blr_model; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* fix non-linear response */
  if(parset.flag_nonlinear !=1)
  {
    par_fix[idx_resp + 1] = 1;
    par_fix_val[idx_resp + 1] = 0.0;
  }

  /* fix FA, no need in SARM */
  par_fix[num_params_blr + num_params_sa_blr_model+2] = 1;
  par_fix_val[num_params_blr + num_params_sa_blr_model+2] = log(1.0);

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err != 1)
  {
    par_fix[num_params_blr_tot] = 1;
    par_fix_val[num_params_blr_tot] = log(1.0);
  }

  /* fix continuum variation parameter sigma and tau if flag_fixvar is true */
  if(parset.flag_fixvar == 1)
  {
    par_fix[num_params_blr_tot + 1] = 1;
    par_fix_val[num_params_blr_tot + 1] = var_param[1];
    par_fix[num_params_blr_tot + 2] = 1;
    par_fix_val[num_params_blr_tot + 2] = var_param[2];
  }

  print_par_names_sarm();

  force_update = parset.flag_force_update;
  if(parset.flag_para_name != 1)
  {
    strcpy(dnest_data_dir,parset.file_dir);
    strcat(dnest_data_dir, "/data/");
    logz_sarm = dnest(argc, argv, fptrset_sarm, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);
  }

  dnest_free_fptrset(fptrset_sarm);

  return 0;
}

void set_par_range_sarm()
{
  int i, j, i1, i2;

  /* SA BLR parameters */
  for(i=num_params_blr; i<num_params_blr + num_params_sa_blr_model; i++)
  {
    par_range_model[i][0] = sa_blr_range_model[i-num_params_blr][0];
    par_range_model[i][1] = sa_blr_range_model[i-num_params_blr][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  /* SA extra parameters */
  for(i=num_params_blr + num_params_sa_blr_model; i<num_params_blr_tot; i++)
  {
    par_range_model[i][0] = sa_extpar_range[i-(num_params_blr + num_params_sa_blr_model)][0];
    par_range_model[i][1] = sa_extpar_range[i-(num_params_blr + num_params_sa_blr_model)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* variability parameters */
   /* first systematic error */
  i = num_params_blr_tot;
  par_range_model[i][0] = var_range_model[0][0];
  par_range_model[i][1] = var_range_model[0][1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  for(i=num_params_blr_tot+1; i<num_params_drw + num_params_blr_tot; i++)
  {
    if(var_param_std[i-num_params_blr_tot] > 0.0)
    {
      par_range_model[i][0] = var_param[i-num_params_blr_tot] - 5.0 * var_param_std[i-num_params_blr_tot];
      par_range_model[i][1] = var_param[i-num_params_blr_tot] + 5.0 * var_param_std[i-num_params_blr_tot];

      /* make sure that the range lies within the initial range */
      par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[i-num_params_blr_tot][0]);
      par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[i-num_params_blr_tot][1]);

      par_prior_model[i] = GAUSSIAN;
      par_prior_gaussian[i][0] = var_param[i-num_params_blr_tot];
      par_prior_gaussian[i][1] = var_param_std[i-num_params_blr_tot];
    }
    else
    {
      par_range_model[i][0] = var_range_model[i-num_params_blr_tot][0];
      par_range_model[i][1] = var_range_model[i-num_params_blr_tot][1];

      par_prior_model[i] = UNIFORM;
      par_prior_gaussian[i][0] = 0.0;
      par_prior_gaussian[i][1] = 0.0;
    }
  }

  /* long-term trend of continuum */
  for(i=num_params_drw + num_params_blr_tot; i< num_params_drw + num_params_trend + num_params_blr_tot; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }

  /* response A and Ag */
  j = 0;
  i1 = idx_resp;
  i2 = idx_resp + num_params_resp;
  for(i=i1; i<i2; i++)
  {
    par_range_model[i][0] = resp_range[j][0];
    par_range_model[i][1] = resp_range[j][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;

    j++;
  }

  /* different trend in continuum and line */
  j = 0;
  i1 = idx_difftrend;
  i2 = idx_difftrend + num_params_difftrend;
  for(i= i1; i< i2; i++)
  {
    par_range_model[i][0] = var_range_model[4 + j][0];
    par_range_model[i][1] = var_range_model[4 + j][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;

    j++;
  }

  /* continuum ligth curve values */
  for(i=num_params_blr_tot+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[4+num_params_difftrend][0];
    par_range_model[i][1] = var_range_model[4+num_params_difftrend][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }

  return;
}


void print_par_names_sarm()
{
  if(thistask != roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH], str_name[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_sarm.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  strcpy(str_fmt, "%4d %-28s %10.6f %10.6f %4d %4d %15.6e\n");

  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-28s %10s %10s %4s %4s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val");

  i=-1;
  for(j=0; j<num_params_sa_blr_model; j++)
  {
    i++;
    strcpy(str_name, "\0");
    if(BLRmodel_sa_name!=NULL && BLRmodel_sa_name[i] != NULL)
    {
      strcpy(str_name, "SA BLR model ");
      strcat(str_name, BLRmodel_sa_name[i]);
    }
    else 
    {
      strcpy(str_name, "SA BLR model");
    }
    fprintf(fp, str_fmt, i, str_name, par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  for(j=0; j<num_params_sa_extpar; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "SA Extra Par", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  i++;
  fprintf(fp, str_fmt, i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  i++;
  fprintf(fp, str_fmt, i, "sigmad", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  i++;
  fprintf(fp, str_fmt, i, "taud", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  
  for(j=0; j<num_params_trend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  i++;
  fprintf(fp, str_fmt, i, "A", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);

  i++;
  fprintf(fp, str_fmt, i, "Ag", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
                            
  for(j=0; j<num_params_difftrend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "diff trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  for(j=0; j<parset.n_con_recon; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "time series", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }
  
  fclose(fp);
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_sarm(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    if(par_prior_model[i] == GAUSSIAN)
    {
      pm[i] = dnest_randn()*par_prior_gaussian[i][1] + par_prior_gaussian[i][0];
      dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
    }
    else
    {
      pm[i] = par_range_model[i][0] + dnest_rand() * (par_range_model[i][1] - par_range_model[i][0]);
    }
  }

  /* cope with fixed parameters */
  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }

  which_parameter_update = -1;

  return;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_initial_sarm(const void *model)
{
  double logL;
  logL = prob_initial_sarm(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_sarm(const void *model)
{
  double logL;
  logL = prob_initial_sarm(model);
  return logL;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_sarm(const void *model)
{
  double logL;
  logL = prob_sarm(model);
  return logL;
}

/*!
 * this function perturbs parameters.
 */
double perturb_sarm(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level, size_levels, count_saves; 

  /* 
   * fixed parameters need not to update 
   * perturb important parameters more frequently
   */
  do
  {
    rnd = dnest_rand();
    if(rnd < rnd_frac)
      which = dnest_rand_int(num_params_blr_tot + num_params_var);
    else
      which = dnest_rand_int(parset.n_con_recon) + num_params_blr_tot + num_params_var;

  }while(par_fix[which] == 1);
  
  which_parameter_update = which;
  
  /* level-dependent width */
  count_saves = dnest_get_count_saves();
  which_level_update = dnest_get_which_level_update();
  size_levels = dnest_get_size_levels();
  which_level = which_level_update > (size_levels-10)?(size_levels-10):which_level_update;

  if( which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = (limit2 - limit1);
  }
  else
  {
    limit1 = par_range_model[which][0];
    limit2 = par_range_model[which][1];
    width = (par_range_model[which][1] - par_range_model[which][0]);
  }
  width /= (2.35);

  if(par_prior_model[which] == GAUSSIAN)
  {
    logH -= (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    //dnest_wrap(&pm[which], limit1, limit2);
    logH += (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    //dnest_wrap(&pm[which], limit1, limit2);
  }

  return logH;
}

/*!
 * this function prints out parameters.
 */
void print_particle_sarm(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

/*!
 * This function read the particle from the file.
 */
void read_particle_sarm(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file sample_sa2d.txt.\n");
      exit(0);
    }
  }
  return;
}

void accept_action_sarm()
{
  int param;
  double *ptemp;

  // the parameter previously updated
  param = which_parameter_update;

  /* continuum parameter is updated */
  if( param >= num_params_blr_tot )
  {
    /* 
     *note that (response) Fline is also changed as long as Fcon is changed.
     *num_params_blr-the parameter is the systematic error of continuum.
     *the change of this parameter also changes continuum reconstruction.
     */

    ptemp = Fcon_rm_particles[which_particle_update];
    Fcon_rm_particles[which_particle_update] = Fcon_rm_particles_perturb[which_particle_update];
    Fcon_rm_particles_perturb[which_particle_update] = ptemp;
  }
  else if((which_parameter_update < num_params_blr_tot) && force_update != 1 )
  {
    /* RM and SA BLR parameter is updated 
     * Note a) that the (num_par_blr-1)-th parameter is systematic error of line.
     * when this parameter is updated, Trans2D, Fline, phase_sa, and Fline_sa are unchanged.
     *      b) Fline is always changed, except param = num_params_blr-1.
     */
    
    ptemp = TransTau_particles[which_particle_update];
    TransTau_particles[which_particle_update] = TransTau_particles_perturb[which_particle_update];
    TransTau_particles_perturb[which_particle_update] = ptemp;

    ptemp = Trans2D_at_veldata_particles[which_particle_update];
    Trans2D_at_veldata_particles[which_particle_update] = Trans2D_at_veldata_particles_perturb[which_particle_update];
    Trans2D_at_veldata_particles_perturb[which_particle_update] = ptemp;

    ptemp = Trans_alpha_at_veldata_particles[which_particle_update];
    Trans_alpha_at_veldata_particles[which_particle_update] = Trans_alpha_at_veldata_particles_perturb[which_particle_update];
    Trans_alpha_at_veldata_particles_perturb[which_particle_update] = ptemp;

    ptemp = Trans_beta_at_veldata_particles[which_particle_update];
    Trans_beta_at_veldata_particles[which_particle_update] = Trans_beta_at_veldata_particles_perturb[which_particle_update];
    Trans_beta_at_veldata_particles_perturb[which_particle_update] = ptemp;
  }
  
  return;
}

void kill_action_sarm(int i, int i_copy)
{
  memcpy(Fcon_rm_particles[i], Fcon_rm_particles[i_copy], parset.n_con_recon * sizeof(double));
  memcpy(TransTau_particles[i], TransTau_particles[i_copy], parset.n_tau*sizeof(double));
  memcpy(Trans2D_at_veldata_particles[i], Trans2D_at_veldata_particles[i_copy], parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
  memcpy(Trans_alpha_at_veldata_particles[i], Trans_alpha_at_veldata_particles[i_copy], parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
  memcpy(Trans_beta_at_veldata_particles[i], Trans_beta_at_veldata_particles[i_copy], parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
}
#endif