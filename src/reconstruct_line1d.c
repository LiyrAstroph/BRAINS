/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_interp.h>

#include "dnest_line1d.h"
#include "allvars.h"
#include "proto.h"

void reconstruct_line1d()
{
  char *argv[]={""};

  reconstruct_line1d_init();

  dnest_line1d(0, argv);

  if(thistask == roottask)
  {
    calculate_con_from_model(best_model_line1d + 8 *sizeof(double));
    FILE *fp;
    char fname[200];
    int i;

    sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tcon[i], Fcon[i] / con_scale);
    }
    fclose(fp);

    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;
    beta_old_particles[which_particle_update] = -1.0;
    Fcon = Fcon_particles[which_particle_update];
    Trans1D = Trans1D_particles[which_particle_update];
    transfun_1d_cloud_direct(best_model_line1d);
    calculate_line_from_blrmodel(best_model_line1d, Tline, Fline, parset.n_line_recon);

    // output reconstructed line light curve
    sprintf(fname, "%s/%s", parset.file_dir, parset.pline_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_line_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tline[i], Fline[i] / line_scale);
    }
    fclose(fp);

    // output transfer function.
    sprintf(fname, "%s/%s", parset.file_dir, parset.tran_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_tau; i++)
    {
      fprintf(fp, "%f %f\n", TransTau[i], Trans1D[i]);
    }
    fclose(fp);
  }

  reconstruct_line1d_end();

}

void reconstruct_line1d_init()
{
  TransTau = malloc(parset.n_tau * sizeof(double));
  //Trans1D = malloc(parset.n_tau * sizeof(double));

  Fline_at_data = malloc(n_line_data * sizeof(double));

  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Flerrs = malloc(parset.n_line_recon * sizeof(double));

  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 20);
  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 20);

  int i;
  double dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }

  dTransTau = (parset.tau_max_set - parset.tau_min_set)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = parset.tau_min_set + dTransTau * i;
  }

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS1D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  Trans1D_particles = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans1D_particles[i] = malloc(parset.n_tau * sizeof(double));
  }

  // only record gamma-distribution random number of clouds
  clouds_particles = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    clouds_particles[i] = malloc(parset.n_cloud_per_task * sizeof(double));
  }

  beta_old_particles = malloc(parset.num_particles * sizeof(double));
  for(i=0; i<parset.num_particles; i++)
    beta_old_particles[i] = -1.0;

  Fcon_particles = malloc(parset.num_particles * sizeof(double *));
  Fcon_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fcon_particles[i] = malloc(parset.n_con_recon * sizeof(double));
    Fcon_particles_perturb[i] = malloc(parset.n_con_recon * sizeof(double));
  }

  model_old_particles = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    model_old_particles[i] = malloc(size_of_modeltype);
  }

  perturb_accept = malloc(parset.num_particles * sizeof(int));
}

void reconstruct_line1d_end()
{
  free(TransTau);
  free(Fline_at_data);

  free(Tline);
  free(Fline);
  free(Flerrs);

  int i;
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans1D_particles[i]);
    free(clouds_particles[i]);
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
    free(model_old_particles[i]);
  }
  free(Trans1D_particles);
  free(clouds_particles);
  free(Fcon_particles);
  free(Fcon_particles_perturb);
  free(model_old_particles);

  free(perturb_accept);

  free(beta_old_particles);

  if(thistask == roottask)
  {
    printf("Ends reconstruct_line1d.\n");
  }
}

double prob_line1d(const void *model)
{
  double prob = 0.0, fcon, var2, dy;
  int i;
  
  // if the previous perturb is accepted, store the previous Fcon at perturb stage;
  // otherwise, Fcon has no changes;
  if(perturb_accept[which_particle_update] == 1)
  {
    memcpy(Fcon_particles[which_particle_update], Fcon_particles_perturb[which_particle_update], 
      parset.n_con_recon*sizeof(double));
  }
  
  // only udate continuum reconstruction when the corresponding parameters are updated
  // or force to update (which_parameter_update = -1)
  if((which_parameter_update >= 8 ) || which_parameter_update == -1)
  {
    Fcon = Fcon_particles_perturb[which_particle_update];
    calculate_con_from_model(model + 8*sizeof(double));
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  }
  else
  {
    Fcon = Fcon_particles[which_particle_update];
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  }

  for(i=0; i<n_con_data; i++)
  {
    fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
    prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - ( 0.5*log(2.0*PI) + log(Fcerrs_data[i]) );
  }

  // only update transfer function when BLR model is changed
  // or force to update (which_parameter_update = -1)
  // Trans1D is a pointer to the transfer function
  if(which_parameter_update < 8 || which_parameter_update == -1)
  {
    Trans1D = Trans1D_particles[which_particle_update]; 
    transfun_1d_cloud_direct(model);
    //memcpy(Trans1D_particles[which_particle_update], Trans1D, parset.n_tau*sizeof(double));
  }
  else
  {
    Trans1D = Trans1D_particles[which_particle_update];
    //memcpy(Trans1D, Trans1D_particles[which_particle_update], parset.n_tau*sizeof(double));
  }

  calculate_line_from_blrmodel(model, Tline_data, Fline_at_data, n_line_data);
  for(i=0; i<n_line_data; i++)
  {
    dy = Fline_data[i] - Fline_at_data[i] ;
    var2 = Flerrs_data[i]*Flerrs_data[i];
    prob += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
  }

  return prob;
}