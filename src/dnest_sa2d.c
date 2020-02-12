/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifdef SA

#include "brains.h"

void set_par_range_sa2d()
{
  int i;

  /* setup parameter range, BLR parameters first */
  for(i=0; i<num_params_blr_model; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* cope with narrow line */
  for(i=num_params_blr_model; i<num_params_blr_model + num_params_nlr; i++)
  {
    par_range_model[i][0] = nlr_range_model[i - num_params_blr_model][0];
    par_range_model[i][1] = nlr_range_model[i - num_params_blr_model][1];
    
    par_prior_model[i] = nlr_prior_model[i - num_params_blr_model];
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0; /* note that for logarithm prior of flux, this value is not used, so does not matter */
  }
  /* cope with spectral broadening */
  for(i=num_params_blr_model + num_params_nlr; i<num_params_blr_model + num_params_nlr + num_params_res; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  /* cope with line center */
  for(i=num_params_blr-num_params_linecenter-3; i< num_params_blr-3; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  
  for(i=num_params_blr-3; i<num_params_blr-1; i++)
  {
    par_range_model[i][0] = resp_range[i - (num_params_blr-3)][0];
    par_range_model[i][1] = resp_range[i - (num_params_blr-3)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* the last is systematic error */
  i = num_params_blr-1;
  par_range_model[i][0] = sys_err_line_range[0];
  par_range_model[i][1] = sys_err_line_range[1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  /* variability parameters */
   /* first systematic error */
  i = num_params_blr;
  par_range_model[i][0] = var_range_model[0][0];
  par_range_model[i][1] = var_range_model[0][1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  for(i=num_params_blr+1; i<num_params_drw + num_params_blr; i++)
  {
    if(var_param_std[i-num_params_blr] > 0.0)
    {
      par_range_model[i][0] = var_param[i-num_params_blr] - 5.0 * var_param_std[i-num_params_blr];
      par_range_model[i][1] = var_param[i-num_params_blr] + 5.0 * var_param_std[i-num_params_blr];

      /* make sure that the range lies within the initial range */
      par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[i-num_params_blr][0]);
      par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[i-num_params_blr][1]);

      par_prior_model[i] = GAUSSIAN;
      par_prior_gaussian[i][0] = var_param[i-num_params_blr];
      par_prior_gaussian[i][1] = var_param_std[i-num_params_blr];
    }
    else
    {
      par_range_model[i][0] = var_range_model[i-num_params_blr][0];
      par_range_model[i][1] = var_range_model[i-num_params_blr][1];

      par_prior_model[i] = UNIFORM;
      par_prior_gaussian[i][0] = 0.0;
      par_prior_gaussian[i][1] = 0.0;
    }
  }
  /* long-term trend of continuum */
  for(i=num_params_drw + num_params_blr; i< num_params_drw + num_params_trend + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  /* different trend in continuum and line */
  for(i=num_params_drw + num_params_trend + num_params_blr; i< num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr)][0];
    par_range_model[i][1] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* continuum ligth curve values */
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[4+num_params_difftrend][0];
    par_range_model[i][1] = var_range_model[4+num_params_difftrend][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }

  return;
}

/*!
 *  print names and prior ranges for parameters 
 *
 */
void print_par_names_sa2d()
{
  if(thistask != roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_model2d.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  printf("# Print parameter name in %s\n", fname);


  i=-1;
  for(j=0; j<num_params_blr_model; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "BLR model", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<num_params_nlr; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "narrow line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }
  for(j=0; j<num_params_res; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "line broaden", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }
  for(j=0; j<num_params_linecenter; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "line center", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }
  
  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "A", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "Ag", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sys_err_line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sigmad", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "taud", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  
  for(j=0; j<num_params_trend; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<num_params_difftrend; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "diff trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<parset.n_con_recon; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "time series", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }
  
  fclose(fp);
}
#endif