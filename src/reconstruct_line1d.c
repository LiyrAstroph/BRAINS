/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_interp.h>

#include "dnest_line1d.h"
#include "allvars.h"
#include "proto.h"

void reconstruct_line1d()
{
  reconstruct_line1d_init();

  dnest_line1d(0, " ");

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
  Trans1D = malloc(parset.n_tau * sizeof(double));

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
}

double prob_line1d(void *model)
{
  double prob = 0.0, fcon, var2, dy;
  int i;
  
  calculate_con_from_model(model + 8*sizeof(double));

  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

  for(i=0; i<n_con_data; i++)
  {
    fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
    prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - ( 0.5*log(2.0*PI) + log(Fcerrs_data[i]) );
  }

  // only update transfer function when BLR model is changed.
  transfun_1d_cloud_direct(model);
  calculate_line_from_blrmodel(model, Tline_data, Fline_at_data, n_line_data);
  for(i=0; i<n_line_data; i++)
  {
    dy = Fline_data[i] - Fline_at_data[i] ;
    var2 = Flerrs_data[i]*Flerrs_data[i];
    var2 = 1.0*var2;
    prob += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
  }

  return prob;
}

void reconstruct_line1d_end()
{
  free(TransTau);
  free(Trans1D);
  free(Fline_at_data);

  free(Tline);
  free(Fline);
  free(Flerrs);
}