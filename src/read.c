/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file read.c
 *  \brief read configuration file and data files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*!
 * read parameter set from parameter file.
 */
void read_parset()
{
  if(thistask == roottask)
  {
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    nt = 0;
    strcpy(tag[nt], "FileDir");
    addr[nt] = &parset.file_dir;
    id[nt++] = STRING;

    strcpy(tag[nt], "ContinnumFile");
    addr[nt] = &parset.continuum_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "LineFile");
    addr[nt] = &parset.line_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "Line2DFile");
    addr[nt] = &parset.line2d_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "FlagDim");
    addr[nt] = &parset.flag_dim;
    id[nt++] = INT;

    strcpy(tag[nt], "NConRecon");
    addr[nt] = &parset.n_con_recon;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagTrend");
    addr[nt] = &parset.flag_trend;
    id[nt++] = INT;

    strcpy(tag[nt], "ConConstructFileOut");
    addr[nt] = &parset.pcon_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "NLineRecon");
    addr[nt] = &parset.n_line_recon;
    id[nt++] = INT;

    strcpy(tag[nt], "NVelRecon");
    addr[nt] = &parset.n_vel_recon;
    id[nt++] = INT;

    strcpy(tag[nt], "LineConstructFileOut");
    addr[nt] = &parset.pline_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "TranFileOut");
    addr[nt] = &parset.tran_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "Line2DConstructFileOut");
    addr[nt] = &parset.pline2d_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "Line2DDataConstructFileOut");
    addr[nt] = &parset.pline2d_data_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "Tran2DFileOut");
    addr[nt] = &parset.tran2d_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "Tran2DDataFileOut");
    addr[nt] = &parset.tran2d_data_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "NCloudPerCore");
    addr[nt] = &parset.n_cloud_per_task;
    id[nt++] = INT;

    strcpy(tag[nt], "NVPerCloud");
    addr[nt] = &parset.n_vel_per_cloud;
    id[nt++] = INT;

    strcpy(tag[nt], "NTau");
    addr[nt] = &parset.n_tau;
    id[nt++] = INT;

    strcpy(tag[nt], "TauMinSet");
    addr[nt] = &parset.tau_min_set;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "TauMaxSet");
    addr[nt] = &parset.tau_max_set;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "FlagCloudsOut");
    addr[nt] = &parset.flag_save_clouds;
    id[nt++] = INT;

    strcpy(tag[nt], "CloudsFileOut");
    addr[nt] = &parset.cloud_out_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "InstRes");
    addr[nt] = &parset.InstRes;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "FlagNarrowLine");
    addr[nt] = &parset.flag_narrowline;
    id[nt++] = INT;

    strcpy(tag[nt], "FluxNarrowLine");
    addr[nt] = &parset.flux_narrowline;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "WidthNarrowLine");
    addr[nt] = &parset.width_narrowline;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ShiftNarrowLine");
    addr[nt] = &parset.shift_narrowline;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "BLRParFix");
    addr[nt] = &parset.str_par_fix;
    id[nt++] = STRING;

    strcpy(tag[nt], "BLRParFixVal");
    addr[nt] = &parset.str_par_fix_val;
    id[nt++] = STRING;

    strcpy(tag[nt], "FlagFixVar");
    addr[nt] = &parset.flag_fixvar;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagBLRModel");
    addr[nt] = &parset.flag_blrmodel;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagTrendDiff");
    addr[nt] = &parset.flag_trend_diff;
    id[nt++] = INT;


    /* default values */
    parset.flag_dim = 0;
    parset.flag_trend = 0;
    parset.flag_narrowline = 0;
    parset.flag_fixvar = 0;
    parset.flag_blrmodel = 1;
    parset.flag_trend_diff = 0;
    
    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='%')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
        {
          j = i;
          tag[i][0] = 0;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(id[j])
        {
          case DOUBLE:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy(addr[j], buf2);
            break;
          case INT:
            *((int *)addr[j]) = (int) atof(buf2);
            break;
        }
      }
      else
      {
        fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                      parset.param_file, buf1);
        exit(0);
      }
    }
    fclose(fparam);

    if(parset.flag_dim > 2 || parset.flag_dim < -1)
    {
      fprintf(stderr, "# Error in flag_dim: value %d is not allowed.\nPlease specify a value in [-1-2].\n", parset.flag_dim);
      exit(0);
    }

    if(parset.flag_trend > 2 || parset.flag_trend < 0)
    {
      fprintf(stderr, "# Error in flag_trend: value %d is not allowed.\nPlease specify a value in [0-1].\n", parset.flag_trend);
      exit(0);
    }

    if(parset.flag_narrowline == 1 && parset.flag_dim == 2)
    {
      printf("# add narrow-line: flux=%e, width=%fkm/s, shift=%fkm/s.\n", parset.flux_narrowline, 
        parset.width_narrowline, parset.shift_narrowline);
      parset.width_narrowline /= VelUnit;
      parset.shift_narrowline /= VelUnit;
    }

    if(parset.flag_blrmodel == 3 || parset.flag_blrmodel == 4)
    {
      parset.n_vel_per_cloud = 1;
    }

    if(parset.flag_dim < 1)
    {
      parset.flag_trend_diff = 0;
    }

    if(parset.flag_dim < 0) // create mock data
    {
      // set large values
      parset.n_cloud_per_task = fmax(2.0e5, parset.n_cloud_per_task);
      parset.n_vel_per_cloud = fmax(10.0, parset.n_vel_per_cloud);
      printf("# set NCloudPerCore and NVPerCloud: %d %d\n", parset.n_cloud_per_task, parset.n_vel_per_cloud);
    }
  }
  
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}

/*!
 * read dataset.
 */
void read_data()
{
  FILE *fp;
  int i;
  char buf[200], fname[200];

  // first need to determine the number of data points 
  if(thistask == roottask)
  {
    
    int count;
    
    // continuum file
    sprintf(fname, "%s/%s", parset.file_dir, parset.continuum_file);
    fp = fopen(fname, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    // count the number of lines
    count = 0;
    while(1)
    {
      fgets(buf, 200, fp);
      if(feof(fp)!=0)
        break;
      count++;
    }
    fclose(fp);
    n_con_data = count;
    
    printf("continuum data points: %d\n", n_con_data);

    if(parset.flag_dim == 1)
    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.line_file);
    // emission flux line
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      // count the number of lines
      count = 0;
      while(1)
      {
        fgets(buf, 200, fp);
        if(feof(fp)!=0)
          break;
        count++;
      }
      fclose(fp);
      n_line_data = count;
      printf("line data points: %d\n", n_line_data);
    }

    if(parset.flag_dim == 2)
    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.line2d_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      fscanf(fp, "%s %d %d\n", buf, &n_line_data, &n_vel_data);
      fclose(fp);
      printf("line2d data points: %d %d\n", n_line_data, n_vel_data);
    }
  }

  MPI_Bcast(&n_con_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  if(parset.flag_dim == 1)
    MPI_Bcast(&n_line_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  if(parset.flag_dim == 2)
  {
    MPI_Bcast(&n_line_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&n_vel_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  }

  // now allocate memory for data
  allocate_memory_data();

  // now read data
  if(thistask == roottask)
  {
    // continuum data
    sprintf(fname, "%s/%s", parset.file_dir, parset.continuum_file);
    fp = fopen(fname, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    for(i=0; i<n_con_data; i++)
    {
      fscanf(fp, "%lf %lf %lf \n", &Tcon_data[i], &Fcon_data[i], &Fcerrs_data[i]);
    }
    fclose(fp);

  }

  MPI_Bcast(Tcon_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(Fcon_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(Fcerrs_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  // read line
  if(parset.flag_dim == 1)
  {
    if(thistask == roottask)
    {
    // line flux data
      sprintf(fname, "%s/%s", parset.file_dir, parset.line_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<n_line_data; i++)
      {
        fscanf(fp, "%lf %lf %lf \n", &Tline_data[i], &Fline_data[i], &Flerrs_data[i]);
      }
      fclose(fp);

      if(Tline_data[0] - Tcon_data[n_con_data-1] > 0.0)
      {
        fprintf(stderr, "# Error: No time overlap between ontinuum and line time series.\n");
        fprintf(stderr, "# Error: continuum, %f-%f; line, %f-%f.\n", Tcon_data[0], Tcon_data[n_con_data-1], 
          Tline_data[0], Tline_data[n_line_data-1]);
        exit(-1);
      }

    // cal mean line error
      line_error_mean = 0.0;
      for(i=0; i<n_line_data;i++)
      {
        line_error_mean += Flerrs_data[i];
      }
      line_error_mean /= n_line_data;
    }

    MPI_Bcast(Tline_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fline_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Flerrs_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

    MPI_Bcast(&line_error_mean, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }

  // read 2d line data
  if(parset.flag_dim == 2)
  {
    if(thistask == roottask)
    {
      int j;

      sprintf(fname, "%s/%s", parset.file_dir, parset.line2d_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      fgets(buf, 200, fp);

      for(i=0; i<n_line_data; i++)
      {
        for(j=0; j<n_vel_data; j++)
        {
          fscanf(fp, "%lf%lf%lf%lf\n", &Vline_data[j], &Tline_data[i], 
                 &Fline2d_data[i*n_vel_data + j], &Flerrs2d_data[i*n_vel_data + j]);

          Vline_data[j] /= VelUnit;
        }

        fscanf(fp, "\n");
      }
      fclose(fp);

      if(Tline_data[0] - Tcon_data[n_con_data-1] > 0.0)
      {
        fprintf(stderr, "# Error: No time overlap between ontinuum and line time series.\n");
        fprintf(stderr, "# Error: continuum, %f-%f; line, %f-%f.\n", Tcon_data[0], Tcon_data[n_con_data-1], 
          Tline_data[0], Tline_data[n_line_data-1]);
        exit(-1);
      }
      
      // cal mean line error
      line_error_mean = 0.0;
      for(i=0; i<n_line_data;i++)
        for(j=0; j<n_vel_data; j++)
        {
          line_error_mean += Flerrs2d_data[i];
        }

      line_error_mean /= (n_line_data*n_vel_data);
    }

    // broadcast 2d data
    MPI_Bcast(Vline_data, n_vel_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Tline_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fline2d_data, n_line_data*n_vel_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Flerrs2d_data, n_line_data*n_vel_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

    MPI_Bcast(&line_error_mean, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    // each task calculates line fluxes
    cal_emission_flux();
  }

  return;
}

/*!
 * allocate data memory. 
 */
void allocate_memory_data()
{
  Tcon_data = malloc(n_con_data * sizeof(double));
  Fcon_data = malloc(n_con_data * sizeof(double));
  Fcerrs_data = malloc(n_con_data * sizeof(double));

  if(parset.flag_dim == 1)
  {
    Tline_data = malloc(n_line_data * sizeof(double));
    Fline_data = malloc(n_line_data * sizeof(double));
    Flerrs_data = malloc(n_line_data * sizeof(double));
  }

  if(parset.flag_dim == 2)
  {
    Vline_data = malloc(n_vel_data * sizeof(double));
    Tline_data = malloc(n_line_data * sizeof(double));
    Fline_data = malloc(n_line_data * sizeof(double));
    Flerrs_data = malloc(n_line_data * sizeof(double));

    Fline2d_data = malloc(n_line_data * n_vel_data * sizeof(double *));
    Flerrs2d_data = malloc(n_line_data * n_vel_data * sizeof(double *));
  }
}

/*!
 * fee data memory.
 */
void free_memory_data()
{
  free(Tcon_data);
  free(Fcon_data);
  free(Fcerrs_data);

  if(parset.flag_dim == 1)
  {
    free(Tline_data);
    free(Fline_data);
    free(Flerrs_data);
  }

  if(parset.flag_dim == 2)
  {
    free(Vline_data);
    free(Tline_data);
    free(Fline_data);
    free(Flerrs_data);

    free(Fline2d_data);
    free(Flerrs2d_data);
  }
}

/*! 
 * calculate the integrated emission line flux.
 */
void cal_emission_flux()
{
  int i, j;
  double dV;
  
  // assume that velocity grid is equally spaced 
  dV = (Vline_data[n_vel_data-1]-Vline_data[0])/(n_vel_data-1);

// using trapezoid formula.
  for(j=0; j<n_line_data; j++)
  { 
    Fline_data[j] = Fline2d_data[j*n_vel_data + 0]/2.0;
    Flerrs_data[j] = (Flerrs2d_data[j*n_vel_data + 0] * Flerrs2d_data[j*n_vel_data + 0])/2.0;
    for(i=1; i<n_vel_data-1; i++)
    {
      Fline_data[j] += Fline2d_data[j*n_vel_data + i];
      Flerrs_data[j] += Flerrs2d_data[j*n_vel_data + i]*Flerrs2d_data[j*n_vel_data + i];
    }
    Fline_data[j] += Fline2d_data[j*n_vel_data + n_vel_data-1]/2.0;
    Flerrs_data[j] += (Flerrs2d_data[j*n_vel_data + n_vel_data-1]*Flerrs2d_data[j*n_vel_data + n_vel_data-1])/2.0;

    Fline_data[j] *= dV;
    Flerrs_data[j] *= dV*dV;
    Flerrs_data[j] = sqrt(Flerrs_data[j]);
  }
}

/*!
 * get number of particles from the option file.
 */
void get_num_particles(char *fname)
{
  FILE *fp;
  char buf[BRAINS_MAX_STR_LENGTH], buf1[BRAINS_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  sscanf(buf, "%d", &parset.num_particles);
  fclose(fp);
}

/*!
 * get file name of posterior sample. 
 */
void get_posterior_sample_file(char *fname, char *samplefile)
{
  FILE *fp;
  char buf[BRAINS_MAX_STR_LENGTH], buf1[BRAINS_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.new_level_interval);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.save_interval);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.thread_steps);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_levels);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.lambda);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.beta);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_saves);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_info_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.levels_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sampler_state_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
  sscanf(buf, "%s", samplefile);
  fclose(fp);
}
