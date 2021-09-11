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

  /* fix FA */
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
    logz_sarm = dnest(argc, argv, fptrset_sarm, num_params, NULL, NULL, NULL, "data/", dnest_options_file, NULL, NULL);

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
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_sarm.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  strcpy(str_fmt, "%4d %-15s %10.6f %10.6f %4d %4d %15.6e\n");

  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-15s %10s %10s %4s %4s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val");

  i=-1;
  for(j=0; j<num_params_sa_blr_model; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "SA BLR model", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
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


#endif