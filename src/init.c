/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file init.c
 *  \brief initialize the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <ctype.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_filter.h>

#include "brains.h"

/*! 
 * This function initialize the program.
 */
void init()
{
  int i, j;

  nq = 1 + parset.flag_trend;
  num_params_trend = nq;
  num_params_drw = 3; /* include systematic error */
  if(parset.flag_dim == 0)
  {
    num_params_resp = 0; /* no A and Ag for 0D */
  }
  else
  {
    num_params_resp = 2; /* response A and Ag */
  }

  if(parset.flag_trend_diff > 0)
  {
    num_params_difftrend = parset.flag_trend_diff; /* differences of the trends between continuum and line */
  }
  else
  {
    num_params_difftrend = 0;
  }
  
  num_params_var = num_params_drw + num_params_trend + num_params_resp + num_params_difftrend;

  /* number of parameters for narrow line, only valid for 2d RM. */
  num_params_nlr = 0;
  if(parset.flag_narrowline >= 2)
    num_params_nlr = 3;

  /* number of parameters for spectral broadening, only valid for 2d RM */
  num_params_res = 1;
  if(parset.flag_InstRes > 1)
    num_params_res = n_line_data;

  /* number of parameters for line center, only valid for 2d RM */
  num_params_linecenter = 0;
  if(parset.flag_linecenter >  0)
    num_params_linecenter = 1;
  else if(parset.flag_linecenter <  0)
    num_params_linecenter = n_line_data;
  
  switch(parset.flag_blrmodel)
  {
    case -1:  /* user defined analytical transfer function */
      BLRmodel_size = num_params_MyTransfun2d * sizeof(double);
      set_blr_range_model = set_par_range_mytransfun;
      break;
    case 0:
      BLRmodel_size = num_params_MyBLRmodel2d * sizeof(double);
      set_blr_range_model = set_blr_range_mymodel;
      break;
    case 1:
      BLRmodel_size = sizeof(BLRmodel1);
      set_blr_range_model = set_blr_range_model1;
      break;
    case 2:
      BLRmodel_size = sizeof(BLRmodel2);
      set_blr_range_model = set_blr_range_model2;
      break;
    case 3:
      BLRmodel_size = sizeof(BLRmodel3);
      set_blr_range_model = set_blr_range_model3;
      break;
    case 4:
      BLRmodel_size = sizeof(BLRmodel4);
      set_blr_range_model = set_blr_range_model4;
      break;
    case 5:
      BLRmodel_size = sizeof(BLRmodel5);
      set_blr_range_model = set_blr_range_model5;
      break;
    case 6:
      BLRmodel_size = sizeof(BLRmodel6);
      set_blr_range_model = set_blr_range_model6;
      break;
    case 7:
      BLRmodel_size = sizeof(BLRmodel7);
      set_blr_range_model = set_blr_range_model7;
      break;
    case 8:
      BLRmodel_size = sizeof(BLRmodel8);
      set_blr_range_model = set_blr_range_model8;
      break;
    case 9:
      BLRmodel_size = sizeof(BLRmodel9);
      set_blr_range_model = set_blr_range_model9;
      break;
    default:
      BLRmodel_size = sizeof(BLRmodel1);
      set_blr_range_model = set_blr_range_model1;
      break;
  }

  /* set bin type, 0: center, 1: left edge */
  if(parset.flag_bintype == 0)
  {
    bin_offset = 0.5;
  }
  else 
  {
    bin_offset = 0.0;
  }

  /* BLR sample along theta direction */
  if(parset.flag_theta_sample == 0)
  {
    theta_sample = theta_sample_outer;
  }
  else 
  {
    theta_sample = theta_sample_inner;
  }

#ifdef SpecAstro
  if(parset.flag_dim > 3 || parset.flag_dim < 0)
  {
    num_params_sa_extpar = sizeof(SAExtPar)/sizeof(double);
    switch(parset.flag_sa_blrmodel)
    {
      case 0:
        SABLRmodel_size = num_params_MyBLRmodel2d * sizeof(double);
        set_sa_blr_range_model = set_blr_range_mymodel;
        break;
      case 1:
        SABLRmodel_size = sizeof(SABLRmodel1);
        set_sa_blr_range_model = set_blr_range_model1;
        break;
      case 2:
        SABLRmodel_size = sizeof(SABLRmodel2);
        set_sa_blr_range_model = set_blr_range_model2;
        break;
      case 3:
        SABLRmodel_size = sizeof(SABLRmodel3);
        set_sa_blr_range_model = set_blr_range_model3;
        break;
      case 4:
        SABLRmodel_size = sizeof(SABLRmodel4);
        set_sa_blr_range_model = set_blr_range_model4;
        break;
      case 5:
        SABLRmodel_size = sizeof(SABLRmodel5);
        set_sa_blr_range_model = set_blr_range_model5;
        break;
      case 6:
        SABLRmodel_size = sizeof(SABLRmodel6);
        set_sa_blr_range_model = set_blr_range_model6;
        break;
      case 7:
        SABLRmodel_size = sizeof(SABLRmodel7);
        set_sa_blr_range_model = set_blr_range_model7;
        break;
      case 8:
        SABLRmodel_size = sizeof(SABLRmodel8);
        set_sa_blr_range_model = set_blr_range_model8;
        break;
      case 9:
        SABLRmodel_size = sizeof(SABLRmodel9);
        set_sa_blr_range_model = set_blr_range_model9;
        break;
      default:
        SABLRmodel_size = sizeof(SABLRmodel1);
        set_sa_blr_range_model = set_blr_range_model1;
        break;
    }
  }
#endif 
  
  if(parset.flag_dim != 3 && parset.flag_dim != 4)
  {
    /* set maximum continuum point */
    n_con_max = parset.n_con_recon;
    if(parset.flag_dim >=-1)
    {
      if(n_con_data > n_con_max)
        n_con_max = n_con_data;
    }
  }

  allocate_memory();

  /* initialize GSL */
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  gsl_blr = gsl_rng_alloc(gsl_T);

#ifndef Debug 
  if(parset.flag_rng_seed != 1)
  {
    gsl_rng_set(gsl_r, time(NULL)+thistask+1350); 
    gsl_rng_set(gsl_blr, time(NULL)+thistask+1350+50); 
  }
  else
  {
    gsl_rng_set(gsl_r, parset.rng_seed+thistask+1350); 
    gsl_rng_set(gsl_blr, parset.rng_seed+thistask+1350+50); 
  }
#else
  if(parset.flag_rng_seed != 1)
  {
    gsl_rng_set(gsl_r, 6666+thistask+1350); 
    gsl_rng_set(gsl_blr, 6666+thistask+1350+50); 
    printf("# debugging, task %d brains random seed %d.\n", thistask, 6666+thistask+1350);
  }
  else
  {
    gsl_rng_set(gsl_r, parset.rng_seed+thistask+1350); 
    gsl_rng_set(gsl_blr, parset.rng_seed+thistask+1350+50); 
  }
#endif

  /* default BH mass range */
  mass_range[0] = 0.1;
  mass_range[1] = 1.0e3;

  /* default rcloud_max_set */
  rcloud_max_set = 1.0e6;
  
  /* line profile  */
  if(parset.flag_dim == 3)
  {
    rcloud_min_set = 0.0;
    rcloud_max_set = 1.0e4;
    if(thistask == roottask)
    {
      printf("Rcloud_min_max_set: %.2f %.2f\n", rcloud_min_set, rcloud_max_set);
    }

    set_nlr_range_model();
    set_blr_range_model(blr_range_model);
  }
  
  /* not line profile or SA */
  if(parset.flag_dim != 3 && parset.flag_dim != 4)
  {
    gsl_acc = gsl_interp_accel_alloc();
    gsl_linear = gsl_interp_alloc(gsl_interp_linear, parset.n_con_recon);

    if(parset.flag_dim >=-1)
    {
      /* set Larr_data */
      for(i=0;i<n_con_data;i++)
      {
        Larr_data[i*nq + 0]=1.0;
        for(j=1; j<nq; j++)
          Larr_data[i*nq + j] = pow(Tcon_data[i], j);
      }
    }
  
    Tcad_data = 1.0;
    Tspan_data = 1.0e4;
    if(parset.flag_dim == 0)
    {
      /* set cadence and time span of data */
      Tspan_data = (Tcon_data[n_con_data -1] - Tcon_data[0]);
      Tspan_data_con = Tspan_data;
      if(Tspan_data < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect epochs in continuum, please check the input data.\n");
          exit(0);
        }
      }
      Tcad_data = Tspan_data;
      for(i=1; i< n_con_data; i++)
      {
        if(Tcad_data > Tcon_data[i] - Tcon_data[i-1])
          Tcad_data = Tcon_data[i] - Tcon_data[i-1];
      }
    }
  
    if( (parset.flag_dim > 0 && parset.flag_dim < 7) || parset.flag_dim == -1)
    {   
      /* set cadence and time span of data */
      Tspan_data_con = (Tcon_data[n_con_data -1] - Tcon_data[0]);
      Tspan_data = (Tline_data[n_line_data -1] - Tcon_data[0]);
      if(Tspan_data < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect epochs in continuum and line, please check the input data.\n");
          exit(0);
        }
      }
      Tcad_data = Tspan_data;
      for(i=1; i< n_con_data; i++)
      {
        if(Tcad_data > Tcon_data[i] - Tcon_data[i-1])
          Tcad_data = Tcon_data[i] - Tcon_data[i-1];
      }
      for(i=1; i< n_line_data; i++)
      {
        if(Tcad_data > Tline_data[i] - Tline_data[i-1])
          Tcad_data = Tline_data[i] - Tline_data[i-1];
      }
  
      /* set mediate time of continuum data */
      Tmed_data = 0.5*(Tcon_data[0] + Tcon_data[n_con_data-1]);
  
      for(i=0; i<num_params_difftrend; i++)
      {
        pow_Tcon_data[i] = (pow(Tcon_data[n_con_data-1]-Tmed_data, i+2) 
                          - pow(Tcon_data[0]-Tmed_data, i+2)) / (i+2) / Tspan_data_con;
      }
      
      /* set time back for continuum reconstruction */
      time_back_set = Tspan_data_con + (Tcon_data[0] - Tline_data[0]);
      time_back_set = fmax(2.0*Tcad_data, time_back_set);
  
      /* make rcloud_max and time_back consistent with each other, rcloud_max has a higher priority */
      double DT=Tcon_data[0] - time_back_set;
      if(parset.rcloud_max > 0.0)
      {
        DT = fmax(DT, Tline_data[0] - parset.rcloud_max*2.0);
      }
      else if(parset.time_back > 0.0) /* neglect when parset.rcloud_max is set */
      { 
        DT = fmax(DT, Tcon_data[0] - parset.time_back);
      }
      time_back_set = Tcon_data[0] - DT;
  
      /* set the range of cloud radial distribution */
      rcloud_min_set = 0.0;
      rcloud_max_set = Tspan_data/2.0;
      
      /* rcloud_max should smaller than  (Tl0 - Tc0)/2 */
      rcloud_max_set = fmin( rcloud_max_set,  (Tline_data[0] - Tcon_data[0] + time_back_set)/2.0 );
  
      if(parset.rcloud_max > 0.0)
        rcloud_max_set = fmin(rcloud_max_set, parset.rcloud_max);

      /* check whether n_con_recon is appropriate */
      double med_cad, med_cad_recon;
      //med_cad = get_mediate_cad(Tcon_data, n_con_data);
      med_cad = (Tcon_data[n_con_data-1] - Tcon_data[0])/(n_con_data -1);
      med_cad_recon = (Tcon_data[n_con_data-1] - Tcon_data[0] + time_back_set)/(parset.n_con_recon-1);
      if(med_cad_recon > 1.05*med_cad)
      {
        if(thistask == roottask)
        {
          // printf("# Too small NConRecon.\n" 
          //      "# Better to change it to %d or set RCloudMax/TimeBack to smaller than %f/%f.\n"
          //      "# Note that 1) too small RCloudMax/TimeBack might result in biases.\n"
          //      "            2) but too large NConRecon will seriously slow down the running.\n", 
          //      (int)(parset.n_con_recon * med_cad_recon/med_cad), 
          //      ((parset.n_con_recon-1) * med_cad - (Tcon_data[n_con_data-1] - Tcon_data[0]))/2.0,
          //      ((parset.n_con_recon-1) * med_cad - (Tcon_data[n_con_data-1] - Tcon_data[0])));
          // printf("\e[1;35m" "# Use '-f' option in command line if want to ignore this check.\n" "\e[0m");
           printf(
               "# Continuum sampling: %f; Reconstruction Sampling: %f.\n"
               "# Note that 1) too small RCloudMax/TimeBack might result in biases.\n"
               "            2) but too large NConRecon will seriously slow down the running.\n", 
               med_cad,  med_cad_recon);
        }
        if(parset.flag_force_run != 1)
        {
          exit(0);
        }
      }
    }

    if(thistask == roottask && parset.flag_dim >=-1)
    {
      printf("RM rcloud_min_max_set: %f %f\n", rcloud_min_set, rcloud_max_set);
    }
    
    // response range A and Ag
    resp_range[0][0] = log(0.1);
    resp_range[0][1] = log(10.0);
  
    resp_range[1][0] = -1.0;
    resp_range[1][1] =  3.0;

    /* set the range of continuum variation  */
    var_range_model[0][0] = log(1.0); /* systematic error in continuum */
    var_range_model[0][1] = log(1.0+10.0);
  
    var_range_model[1][0] = -15.0; /* log(sigma) */
    var_range_model[1][1] = -1.0; 
  
    var_range_model[2][0] = log(fmax(Tcad_data, 0.1)); /* log(tau) */
    var_range_model[2][1] = log(Tspan_data); 
  
    var_range_model[3][0] = -10.0; /* mean value or trend parameter values */
    var_range_model[3][1] =  10.0; 
  
    for(i=0; i<num_params_difftrend; i++)
    {
      var_range_model[4+i][0] = -1.0/pow(Tspan_data, i+1); /* slope of the trend in the differences between contiuum and line */
      var_range_model[4+i][1] =  1.0/pow(Tspan_data, i+1); 
    }
  
    var_range_model[4+num_params_difftrend][0] = -10.0; /* light curve values */
    var_range_model[4+num_params_difftrend][1] =  10.0; 
    
    // range for systematic error
    sys_err_line_range[0] = log(1.0);
    sys_err_line_range[1] = log(1.0+10.0);
    
    set_nlr_range_model();
    set_blr_range_model(blr_range_model);

  }
  
#ifdef SpecAstro
  /* SA */
  if(parset.flag_dim > 3)
  {
    if(parset.flag_dim >= 7) /* SARM, take into account SA line time */
    {
      double time_back_sarm_set;
      /* set cadence and time span of data */
      Tspan_data_con = (Tcon_data[n_con_data -1] - Tcon_data[0]);
      Tspan_data = (Tline_sarm_data[n_epoch_sarm_data -1] - Tcon_data[0]);
      if(Tspan_data < 0.0)
      {
        if(thistask == roottask)
        {
          printf("# Incorrect epochs in continuum and line, please check the input data.\n");
          exit(0);
        }
      }
      Tcad_data = Tspan_data;
      for(i=1; i< n_con_data; i++)
      {
        if(Tcad_data > Tcon_data[i] - Tcon_data[i-1])
          Tcad_data = Tcon_data[i] - Tcon_data[i-1];
      }
      for(i=1; i<n_epoch_sarm_data; i++)
      {
        if(Tcad_data > Tline_sarm_data[i] - Tline_sarm_data[i-1])
          Tcad_data = Tline_sarm_data[i] - Tline_sarm_data[i-1];
      }

      /* set time back for continuum reconstruction */
      time_back_sarm_set = Tspan_data_con + (Tcon_data[0] - Tline_sarm_data[0]);
      time_back_sarm_set = fmax(2.0*Tcad_data, time_back_sarm_set);

      /* make rcloud_max and time_back consistent with each other, rcloud_max has a higher priority */
      double DT=Tcon_data[0] - time_back_sarm_set;
      if(parset.rcloud_max > 0.0)
      {
        DT = fmax(DT, Tline_sarm_data[0] - parset.rcloud_max*2.0);
      }
      else if(parset.time_back > 0.0) /* neglect when parset.rcloud_max is set */
      { 
        DT = fmax(DT, Tcon_data[0] - parset.time_back);
      }
      time_back_sarm_set = Tcon_data[0] - DT;

      time_back_set = fmax(time_back_set, time_back_sarm_set);
  
      /* set the range of cloud radial distribution */
      rcloud_min_set = 0.0;
      rcloud_max_set = Tspan_data/2.0;
      
      /* rcloud_max should smaller than  (Tl0 - Tc0)/2 */
      rcloud_max_set = fmin( rcloud_max_set,  (Tline_sarm_data[0] - Tcon_data[0] + time_back_set)/2.0 );
    }
    
    if(parset.rcloud_max > 0.0)
    {
      rcloud_max_set = fmin(rcloud_max_set, parset.rcloud_max);
    }
    if(thistask == roottask)
      printf("SA rcloud_min_max set: %f %f\n", rcloud_min_set, rcloud_max_set);

    set_sa_blr_range_model(sa_blr_range_model);

    sa_extpar_range[0][0] = log(1.0);  /* DA */
    sa_extpar_range[0][1] = log(1000.0); 

    sa_extpar_range[1][0] = 0.0;         /* PA */
    sa_extpar_range[1][1] = 360.0;

    sa_extpar_range[2][0] = log(0.5);    /* FA  */
    sa_extpar_range[2][1] = log(2.0);

    sa_extpar_range[3][0] = -5.0*(wave_sa_data[1]-wave_sa_data[0]);    /* CO  */
    sa_extpar_range[3][1] =  5.0*(wave_sa_data[1]-wave_sa_data[0]);

    sa_extpar_range[4][0] = -200.0;    /* offset wrt to continuum  */
    sa_extpar_range[4][1] =  200.0;

    sa_extpar_range[5][0] = -200.0;    /* offset wrt to continuum  */
    sa_extpar_range[5][1] =  200.0;
  }
#endif 
}

/*!
 * This function allocates memory for variables used throughout the code. 
 */
void allocate_memory()
{
  int i;

  if(parset.flag_dim != 4)
  {
    if(parset.flag_dim != 3)
    {
      Tcon = malloc(parset.n_con_recon * sizeof(double));
      Fcerrs = malloc(parset.n_con_recon * sizeof(double));
    
      PSmat = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));
      //PNmat = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));
      USmat = malloc(parset.n_con_recon * n_con_data * sizeof(double));
      PSmat_data = malloc(n_con_data * n_con_data * sizeof(double));
      //PNmat_data = malloc(n_con_data * n_con_data * sizeof(double));
      PCmat_data = malloc(n_con_data * n_con_data * sizeof(double));
      IPCmat_data = malloc(n_con_data * n_con_data * sizeof(double));
      PQmat = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));
      PEmat1 = malloc(parset.n_con_recon * n_con_data * sizeof(double));
      PEmat2 = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));
          
      workspace = malloc((n_con_data*nq + 7*n_con_max + nq*nq + nq)*sizeof(double));
      workspace_uv = malloc(2*parset.n_con_recon*sizeof(double));
      Larr_data = malloc(n_con_data*nq*sizeof(double));
      Larr_rec = malloc(parset.n_con_recon*nq*sizeof(double));
      pow_Tcon_data = malloc(num_params_difftrend*sizeof(double));
      
      var_param = malloc(num_params_var * sizeof(double));
      var_param_std = malloc(num_params_var * sizeof(double));
    
      for(i=0; i<num_params_var; i++)
      {
        var_param[i] = var_param_std[i] = 0.0;
      }
    }

    blr_range_model = malloc(BLRmodel_size/sizeof(double) * sizeof(double *));
    for(i=0; i<BLRmodel_size/sizeof(double); i++)
    {
      blr_range_model[i] = malloc(2*sizeof(double));
    }
  }

  if(parset.flag_dim != 0 && parset.flag_dim != 3 && parset.flag_dim != 4)
  {
    /* used to smooth transfer function along delay axis */
    gauss_p = gsl_filter_gaussian_alloc(ngauss);
    hist_in = gsl_vector_alloc(parset.n_tau);
    hist_out = gsl_vector_alloc(parset.n_tau);
  }
  
#ifdef SpecAstro
  if(parset.flag_dim > 3 || parset.flag_dim < 0)
  {
    sa_extpar_range = malloc(num_params_sa_extpar * sizeof(double *));
    for(i=0; i<num_params_sa_extpar; i++)
    {
      sa_extpar_range[i]=malloc(2*sizeof(double));
    }

    sa_blr_range_model = malloc(SABLRmodel_size/sizeof(double) * sizeof(double *));
    for(i=0; i<SABLRmodel_size/sizeof(double); i++)
    {
      sa_blr_range_model[i] = malloc(2*sizeof(double));
    }
    
    if(parset.flag_sa_par_mutual != 0)
    {
      idx_sa_par_mutual = malloc(3*sizeof(int)); 
      idx_rm_par_mutual = malloc(3*sizeof(int));
    }
  }
#endif  

  return;
}

/*! 
 * This function free memory.
 */
void free_memory()
{
  int i;
  
  if(parset.flag_dim != 4)
  {
    if(parset.flag_dim != 3)
    {
      free(Tcon);
      free(Fcerrs);
    
      free(PSmat);
      //(PNmat);
      free(USmat);
      free(PSmat_data);
      //free(PNmat_data);
      free(PCmat_data);
      free(IPCmat_data);
      free(PQmat);
      free(PEmat1);
      free(PEmat2);
    
      free(workspace);
      free(workspace_uv);
      free(Larr_data);
      free(Larr_rec);
      free(pow_Tcon_data);
    
      free(var_param);
      free(var_param_std);

      gsl_interp_accel_free(gsl_acc);
      gsl_interp_free(gsl_linear);
    }

    for(i=0; i<BLRmodel_size/sizeof(double); i++)
    {
      free(blr_range_model[i]);
    }
    free(blr_range_model);
  }

  if(parset.flag_dim != 0 && parset.flag_dim != 3 && parset.flag_dim != 4)
  {
    gsl_filter_gaussian_free(gauss_p);
    gsl_vector_free(hist_in);
    gsl_vector_free(hist_out);
  }
  
#ifdef SpecAstro
  if(parset.flag_dim > 3 || parset.flag_dim < 0)
  {
    for(i=0; i<num_params_sa_extpar; i++)
    {
      free(sa_extpar_range[i]);
    }
    free(sa_extpar_range);

    for(i=0; i<SABLRmodel_size/sizeof(double); i++)
    {
      free(sa_blr_range_model[i]);
    }
    free(sa_blr_range_model);

    if(parset.flag_sa_par_mutual != 0)
    {
      free(idx_sa_par_mutual);
      free(idx_rm_par_mutual);
    }
  }
#endif  

  gsl_rng_free(gsl_r);
  gsl_rng_free(gsl_blr);
  return;
}

/*!
 * This function normalise the light curves to a scale of unity.
 */
void scale_con_line()
{
  int i, j;
  double ave_con, ave_line;
  
  con_scale = 1.0;
  line_scale = 1.0;

  ave_con = 0.0;
  ave_line = 0.0;
  for(i=0; i<n_con_data; i++)
  {
    ave_con += Fcon_data[i];
  }

  ave_con /= n_con_data;
  if(ave_con <= 0.0)
  {
    printf("Error in continuum data, a negative mean %f!\n", ave_con);
    exit(0);
  }
  con_scale = 1.0/ave_con;

  for(i=0; i<n_con_data; i++)
  {
    Fcon_data[i] *=con_scale;
    Fcerrs_data[i] *=con_scale;
  }

  if(thistask == roottask)
    printf("con scale: %e\t%e\n", con_scale, ave_con);
  
  con_error_mean *= con_scale;

  if(parset.flag_dim == 0)
  {
    return;
  }

  for(i=0; i<n_line_data; i++)
  {
    ave_line += Fline_data[i];
  }
  ave_line /=n_line_data;
  if(ave_line <= 0.0)
  {
    printf("Error in line data, a negative mean %f!\n", ave_line);
    exit(0);
  }
  line_scale = 1.0/ave_line;
  
  if(thistask == roottask)
    printf("line scale: %e\t%e\n", line_scale, ave_line);

  for(i=0; i<n_line_data; i++)
  {
    // note mask with error < 0.0
    if(Flerrs_data[i] > 0.0)
    {
      Fline_data[i] *= line_scale;
      Flerrs_data[i] *= line_scale;
    }

    if(parset.flag_dim==2 || parset.flag_dim == -1 || parset.flag_dim == 6)
      for(j=0; j<n_vel_data; j++)
      {
        // note mask with error < 0.0
        if(Flerrs2d_data[i*n_vel_data + j] > 0.0)
        {
          Fline2d_data[i*n_vel_data + j] *= line_scale;
          Flerrs2d_data[i*n_vel_data + j] *= line_scale;
        } 
      } 
  }
  line_error_mean *= line_scale;
  line_error_mean_sq *= (line_scale*line_scale);

  if( (parset.flag_dim==2 || parset.flag_dim == -1 || parset.flag_dim == 6) && parset.flag_narrowline!=0)
  {
    if(parset.flag_narrowline < 3)  /* fixed or Gaussian prior */
    {
      parset.flux_narrowline *= line_scale;
      parset.flux_narrowline_err *= line_scale;
    }
    else /* logarithm prior, scale the limit */
    {
      parset.flux_narrowline_low *= line_scale;
      parset.flux_narrowline_upp *= line_scale;
    }
  }
  return;
}

#ifdef SpecAstro
/*!
 * This function normalise the light curves to a scale of unity.
 */
void scale_con_line_sarm()
{
  int i, j;
  double ave_con, ave_line;
  
  con_scale = 1.0;
  line_sarm_scale = 1.0;

  ave_con = 0.0;
  ave_line = 0.0;
  for(i=0; i<n_con_data; i++)
  {
    ave_con += Fcon_data[i];
  }

  ave_con /= n_con_data;
  if(ave_con <= 0.0)
  {
    printf("Error in continuum data, a negative mean %f!\n", ave_con);
    exit(0);
  }
  con_scale = 1.0/ave_con;

  for(i=0; i<n_con_data; i++)
  {
    Fcon_data[i] *=con_scale;
    Fcerrs_data[i] *=con_scale;
  }

  if(thistask == roottask)
    printf("con scale: %e\t%e\n", con_scale, ave_con);
  
  con_error_mean *= con_scale;

  for(i=0; i<n_epoch_sarm_data; i++)
  {
    ave_line += Fline_sarm_data[i];
  }
  ave_line /=n_epoch_sarm_data;
  if(ave_line <= 0.0)
  {
    printf("Error in line data, a negative mean %f!\n", ave_line);
    exit(0);
  }
  
  line_sarm_scale = 1.0/ave_line;
  sarm_scale_ratio = con_scale/line_sarm_scale;  /* this is needed to cal flux ratio */
  
  if(thistask == roottask)
    printf("sarm line scale: %e\t%e\n", line_sarm_scale, ave_line);

  for(i=0; i<n_epoch_sarm_data; i++)
  {
    // note mask with error < 0.0
    if(Flerrs_sarm_data[i] > 0.0)
    {
      Fcon_sarm_data[i] *= line_sarm_scale;  /* note that sarm continuum also needs to scale */
      Fline_sarm_data[i] *= line_sarm_scale;
      Flerrs_sarm_data[i] *= line_sarm_scale;
    }

    for(j=0; j<n_vel_sarm_data; j++)
    {
      // note mask with error < 0.0
      if(Flerrs2d_sarm_data[i*n_vel_sarm_data + j] > 0.0)
      {
        Fline2d_sarm_data[i*n_vel_sarm_data + j] *= line_sarm_scale;
        Flerrs2d_sarm_data[i*n_vel_sarm_data + j] *= line_sarm_scale;
      } 
    } 
  }
  sarm_line_error_mean *= line_sarm_scale;
  return;
}
#endif

/*!
 * This function copes with parameter fixing.\n
 * Only fix BLR model parameters. 
 */
void set_par_fix_blrmodel()
{
  int i;
  char *pstr;
  
  npar_fix = 0;

  if(thistask == roottask)
  {
    pstr = parset.str_par_fix_val;
    // set the default value if not provided.
    for(i=strlen(parset.str_par_fix); i<num_params_blr_model; i++)
      parset.str_par_fix[i] = '0';

    for(i=0; i<num_params_blr_model; i++)
    {
      if(parset.str_par_fix[i] == '0')
      {
        par_fix[i] = 0;
        par_fix_val[i] = -DBL_MAX;  /* set to be the smallest value */
      }
      else if(parset.str_par_fix[i] == '1')
      {
        if(pstr == NULL)
        {
          printf("# %d-th BLR parameter value is not provided (counting from 0).\n", i);
          exit(0);
        }
        par_fix[i] = 1;
        sscanf(pstr, "%lf", &par_fix_val[i]);
        npar_fix++;
        printf("# %d-th BLR parameter fixed, value= %f.\n", i, par_fix_val[i]);
        pstr = strchr(pstr, ':'); /* values are separated by ":" */
        if(pstr!=NULL)
        {
          pstr++;
        }
      }
      else   // default value
      {
        par_fix[i] = 0;
        par_fix_val[i] = -DBL_MAX;
      }
    }
  }

  MPI_Bcast(par_fix, num_params_blr_model, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(par_fix_val, num_params_blr_model, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  return;
}

/*
 * set 1D BLR model and functions.
 */
void set_blr_model1d()
{
  switch(parset.flag_blrmodel)
  {
    case -1:
      num_params_blr_model = num_params_MyTransfun1d;
      transfun_1d_cal = transfun_1d_cal_mytransfun;
      BLRmodel_name = MyTFmodel_name;
      break;

    case 0: 
      num_params_blr_model = num_params_MyBLRmodel1d;
      gen_cloud_sample = gen_cloud_sample_mymodel;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = MyBLRmodel_name;
      break;

    case 1:
      num_params_blr_model = 6;
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;

    case 2:
      num_params_blr_model = 6;
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel2_name;
      break;

    case 3:
      num_params_blr_model = 6;
      gen_cloud_sample = gen_cloud_sample_model3;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel3_name;
      break;

    case 4:
      num_params_blr_model = 6;
      gen_cloud_sample = gen_cloud_sample_model3;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel4_name;
      break;

    case 5:
      num_params_blr_model = 9;
      gen_cloud_sample = gen_cloud_sample_model5;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel5_name;
      break;

    case 6:
      num_params_blr_model = 8;
      gen_cloud_sample = gen_cloud_sample_model6;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel6_name;
      break;

    case 7:
      num_params_blr_model = 13;
      gen_cloud_sample = gen_cloud_sample_model7;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel7_name;
      break;

    case 8:
      num_params_blr_model = 13;
      gen_cloud_sample = gen_cloud_sample_model8;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel8_name;
      break;
    
    case 9:
      num_params_blr_model = 5;
      gen_cloud_sample = gen_cloud_sample_model9;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel9_name;
      break;

    default:
      num_params_blr_model = 6;
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_1d_cal = transfun_1d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;
  }
  return;
}

/*
 * set 2D BLR model and functions.
 */
void set_blr_model2d()
{
  switch(parset.flag_blrmodel)
  {
    case -1:
      num_params_blr_model = num_params_MyTransfun2d;
      transfun_2d_cal = transfun_2d_cal_mytransfun;
      BLRmodel_name = MyTFmodel_name;
      break;

    case 0:
      num_params_blr_model = num_params_MyBLRmodel2d;
      gen_cloud_sample = gen_cloud_sample_mymodel;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = MyBLRmodel_name;
      break;

    case 1:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;

    case 2:
      num_params_blr_model = sizeof(BLRmodel2)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model2;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel2_name;
      break;

    case 3:
      num_params_blr_model = sizeof(BLRmodel3)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model3;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel3_name;
      break;

    case 4:
      num_params_blr_model = sizeof(BLRmodel4)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model4;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel4_name;
      break;

    case 5:
      num_params_blr_model = sizeof(BLRmodel5)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model5;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel5_name;
      break;
    
    case 6:
      num_params_blr_model = sizeof(BLRmodel6)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model6;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel6_name;
      break;
    
    case 7:
      num_params_blr_model = sizeof(BLRmodel7)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model7;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel7_name;
      break;

    case 8:
      num_params_blr_model = sizeof(BLRmodel8)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model8;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel8_name;
      break;
    
    case 9:
      num_params_blr_model = sizeof(BLRmodel9)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model9;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel9_name;
      break;

    default:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;
  }
  return;
}

/*!
 * set range of narrow line model parameters
 */
void set_nlr_range_model()
{
  if(num_params_nlr > 1)
  {
    if(parset.flag_narrowline == 2) /* Gaussian prior */
    {
      nlr_range_model[0][0] = -10.0;
      nlr_range_model[0][1] =  10.0;

      nlr_prior_model[0] = GAUSSIAN;
    }
    else
    {
      nlr_range_model[0][0] = log(parset.flux_narrowline_low); /* logrithmic prior */
      nlr_range_model[0][1] = log(parset.flux_narrowline_upp);

      nlr_prior_model[0] = UNIFORM;
    }

    nlr_range_model[1][0] = -10.0;
    nlr_range_model[1][1] =  10.0;

    nlr_prior_model[1] = GAUSSIAN;

    nlr_range_model[2][0] = -10.0;
    nlr_range_model[2][1] =  10.0;

    nlr_prior_model[2] = GAUSSIAN;
  }
  return;
}

/*!
 *  load prior ranges for parameters from fname 
 *
 */
void load_par_names(char *fname)
{
  if(thistask!= roottask)
    return;
  
  int i, begin, end;
  FILE *fp;
  char str[BRAINS_MAX_STR_LENGTH], buf[BRAINS_MAX_STR_LENGTH];
  int type, fix;
  double val_min, val_max, val, mean, std;

  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  printf("# Loading parameter prior from %s\n", fname);

  while(!feof(fp))
  {
    fgets(str, BRAINS_MAX_STR_LENGTH, fp);
    sscanf(str, "%s", buf);
    
    if(buf[0] == '#')
      continue;
    
    begin = 0;
    while (isspace((unsigned char) str[begin]))
      begin++;
    end = strlen(str) - 1;
    while ((end >= begin) && isspace((unsigned char) str[end]))
      end--;
    for (i = begin; i <= end; i++)
        str[i - begin] = str[i];
    str[i - begin] = '\0';
    if(strlen(str)==0)
      continue;
    
    /* format: %4d %-28s %10.6f %10.6f %4d %4d %15.6e %15.6e %15.6e*/
    if(sscanf(str,"%d %s %lf %lf %d %d %lf %lf %lf", &i, buf, &val_min, &val_max, &type, &fix, &val, &mean, &std)<9)
    {
      printf("Error in reading %s.\nThe line %s is problematic.\n", fname, str);
      exit(0);
    }
    
    if(i >= num_params)
    {
      printf("Error: the number %d exceed the total number %d of parameters.\n", i, num_params);
      exit(0);
    }

    par_range_model[i][0] = val_min;
    par_range_model[i][1] = val_max;
    par_prior_model[i] = type;
    par_fix[i] = fix;
    par_fix_val[i] = val;
    par_prior_gaussian[i][0] = mean;
    par_prior_gaussian[i][1] = std;
    //printf("%d %-28s %f %f %f %f\n", i, buf, val_min, val_max, mean, std);
  }
  
  fclose(fp);
}

/*!
 * set drw parameter ranges using continuum reconstruction results
 *
 */
void set_drw_par_range()
{
  int i;
  /* note omit of the continuum systematic error parameter */
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

  return;
}

/*!
 * set drw parameter ranges using continuum reconstruction results
 * and compare with the loaded priors
 *
 */
void set_drw_par_range_load()
{
  if(thistask != roottask)
    return;

  int i;
  double upp, low;
  /* note omit of the continuum systematic error parameter */
  for(i=num_params_blr_tot+1; i<num_params_drw + num_params_blr_tot; i++)
  {
    /* nozero std means the parameter is not fixed */
    if(var_param_std[i-num_params_blr_tot] > 0.0)
    {
      /* range from continuum reconstruction */
      low = var_param[i-num_params_blr_tot] - 5.0 * var_param_std[i-num_params_blr_tot];
      upp = var_param[i-num_params_blr_tot] + 5.0 * var_param_std[i-num_params_blr_tot];

      if(par_range_model[i][0] > upp || par_range_model[i][1] < low)
      {
        printf("# Error: the input drw parameter priors are inconsistent with those from continuum reconstructions.\n");
        printf("Input: [%f %f]; ContRecon: [%f %f]\n", par_range_model[i][0], par_range_model[i][1],
                                                       low, upp);
        exit(0);
      }

      /*compare with loaded priors */
      par_range_model[i][0] = fmax(par_range_model[i][0], low);
      par_range_model[i][1] = fmin(par_range_model[i][1], upp);

      par_prior_model[i] = GAUSSIAN;
      par_prior_gaussian[i][0] = var_param[i-num_params_blr_tot];
      par_prior_gaussian[i][1] = var_param_std[i-num_params_blr_tot];
    }
  }
  return;
}

/* 
 * get index of mbh from a BLR model.
 * 
 */
int get_idx_mbh_from_blrmodel(int blrmodel)
{
  int idx = -1;
  switch(blrmodel)
  {
    case 0:
      idx = offsetof(MyBLRmodel, mbh);
      break;

    case 1:
      idx = offsetof(BLRmodel1, mbh);
      break;

    case 2:
      idx = offsetof(BLRmodel2, mbh);
      break;
    
    case 3:
      idx = offsetof(BLRmodel3, mbh);
      break;

    case 4:
      idx = offsetof(BLRmodel4, mbh);
      break;
    
    case 5:
      idx = offsetof(BLRmodel5, mbh);
      break;

    case 6:
      idx = offsetof(BLRmodel6, mbh);
      break;
    
    case 7:
      idx = offsetof(BLRmodel7, mbh);
      break;
    
    case 8:
      idx = offsetof(BLRmodel8, mbh);
      break;
    
    case 9:
      idx = offsetof(BLRmodel9, mbh);
      break;
  }
  return idx / sizeof(double);
}

int get_idx_blrsize_from_blrmodel(int blrmodel)
{
  int idx = -1;
  switch(blrmodel)
  {
    case 0:
      idx = 0;   /* unknown, set to 0 by default*/
      break;

    case 1:
      idx = offsetof(BLRmodel1, mu);
      break;

    case 2:
      idx = offsetof(BLRmodel2, mu);
      break;
    
    case 3:
      idx = offsetof(BLRmodel3, Rin);
      break;

    case 4:
      idx = offsetof(BLRmodel4, Rin);
      break;
    
    case 5:
      idx = offsetof(BLRmodel5, mu);
      break;

    case 6:
      idx = offsetof(BLRmodel6, mu);
      break;
    
    case 7:
      idx = offsetof(BLRmodel7, mu);
      break;
    
    case 8:
      idx = offsetof(BLRmodel8, Rblr);
      break;
    
    case 9:
      idx = offsetof(BLRmodel9, mu);
      break;
  }
  return idx / sizeof(double);
}