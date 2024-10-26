#include <stdio.h>
#include <stdlib.h>

#include "brains.h"

/*!
 *  \file blr_range.c
 *  \brief setting parameter ranges of blr models.
 */

/*
 * setup BLR model parameter range. 
 */

// model 1
void set_blr_range_model1(double **blr_range)
{
  int i;

  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 0.999;
  //inc
  blr_range[i][0] = 0.0;  // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];

  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //lambda
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.5;
  //q
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  //blr_range_model[2][1] = fmin(blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 2
void set_blr_range_model2(double **blr_range)
{
  int i;

  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 0.999;
  //inc
  blr_range[i][0] = 0.0;  // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //sigr
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //sigtheta
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  return;
}

// model 3
void set_blr_range_model3(double **blr_range)
{
  int i;
  
  i = 0;
  //Rin
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //F
  blr_range[i][0] = log(1.0);
  blr_range[i++][1] = log(1.0e2);
  //alpha
  blr_range[i][0] = -3.0;
  blr_range[i++][1] = 3.0;
  //inc
  blr_range[i][0] = 0.0;  // in cosine
  blr_range[i++][1] = 1.0;
  //opn 
  blr_range[i][0] = 0.0;   // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //q
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(blr_range_model[3][1] + blr_range_model[4][1]));
  return;
}

// model 4
void set_blr_range_model4(double **blr_range)
{
  int i;
  
  i = 0;
  //Rin
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //F
  blr_range[i][0] = log(1.0);
  blr_range[i++][1] = log(1.0e2);
  //alpha
  blr_range[i][0] = -3.0;
  blr_range[i++][1] = 3.0;
  //inc
  blr_range[i][0] = 0.0; // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = sqrt(2.0)/2.0;
  //q
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(blr_range_model[3][1] + blr_range_model[4][1]));
  return;
}

// model 5
void set_blr_range_model5(double **blr_range)
{
  int i;
  
  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //Fin
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //Fout
  blr_range[i][0] = log(1.0);
  blr_range[i++][1] = log(10.0);
  //alpha
  blr_range[i][0] = 1.0;
  blr_range[i++][1] = 3.0;
  //inc
  blr_range[i][0] = 0.0; // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in degree
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //beta
  blr_range[i][0] = 1.0;
  blr_range[i++][1] = 5.0;
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //fellip
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //fflow
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //sigr_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //sigr_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //theta_rot
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 90.0;
  //sig_turb
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);

  //rcloud_max_set = fmax(rcloud_max_set, exp(blr_range_model[2][1] + blr_range_model[4][1]));

  return;
}


// model 6
void set_blr_range_model6(double **blr_range)
{
  int i;
  
  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //inc
  blr_range[i][0] = 0.0; //in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //gamma
  blr_range[i][0] = 1.0;
  blr_range[i++][1] = 5.0;
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //fellip
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //fflow
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //sigr_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //sigr_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //theta_rot
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 90.0;
  //sig_turb
  blr_range[i][0] = log(0.0001);
  blr_range[i++][1] = log(0.1);

  //rcloud_max_set = 44.85;
  //blr_range_model[2][1] = fmin(blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 7
void set_blr_range_model7(double **blr_range)
{
  int i;
  

  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //inc
  blr_range[i][0] = 0.0;  // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //gamma
  blr_range[i][0] = 1.0;
  blr_range[i++][1] = 5.0;
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  //fsh
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //mu_un
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta_un
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F_un
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //opn_un
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 90.0;

  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];

  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);
  //fellip
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //fflow
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //sigr_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_circ
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //sigr_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);
  //sigthe_rad
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(1.0);
  //theta_rot
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 90.0;

  //fellip_un
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //fflow_un
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;

  //sig_turb
  blr_range[i][0] = log(0.001);
  blr_range[i++][1] = log(0.1);

  return;
}

// model 8
void set_blr_range_model8(double **blr_range)
{
  int i;
  
  i = 0;
  //theta_min
  blr_range[i][0] = 20.0;     // in degree
  blr_range[i++][1] = 90.0;
  //dtheta_max
  blr_range[i][0] = 0.0;      // in degree
  blr_range[i++][1] = 90.0;
  //r_min
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.25);
  //fr_max
  blr_range[i][0] = 1.0;  //
  blr_range[i++][1] = 10.0;
  //gamma
  blr_range[i][0] = 0.0;  //
  blr_range[i++][1] = 3.0;
  //alpha
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 3.0;
  //lambda
  blr_range[i][0] = -3.0;
  blr_range[i++][1] = 0.0;
  //k
  blr_range[i][0] = -0.5;
  blr_range[i++][1] = 0.5;
  //xi
  blr_range[i][0] = 0.0;
  blr_range[i++][1] = 1.0;
  //Rv
  blr_range[i][0] = log(10.0);
  blr_range[i++][1] = log(50.0);
  //Rblr
  blr_range[i][0] = log(10.0);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //inc
  blr_range[i][0] = 0.0;   // cos(inc)
  blr_range[i++][1] = 1.0;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);

  return;
}

// model 9
void set_blr_range_model9(double **blr_range)
{
  int i;

  i = 0;
  //mu
  blr_range[i][0] = log(0.1);
  blr_range[i++][1] = log(rcloud_max_set*0.5);
  //beta
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 2.0;
  //F
  blr_range[i][0] = 0.001;
  blr_range[i++][1] = 0.999;
  //inc
  blr_range[i][0] = 0.0;  // in cosine
  blr_range[i++][1] = 1.0;
  //opn
  blr_range[i][0] = 0.0;  // in rad
  blr_range[i++][1] = 90.0;
  //eta0
  blr_range[i][0]   = eta_range_model[0][0];  // 
  blr_range[i++][1] = eta_range_model[0][1];
  //eta1
  blr_range[i][0]   = eta_range_model[1][0];  // 
  blr_range[i++][1] = eta_range_model[1][1];
  //eta_alpha
  blr_range[i][0]   = eta_range_model[2][0];  // 
  blr_range[i++][1] = eta_range_model[2][1];
  //mbh
  blr_range[i][0] = log(mass_range[0]);
  blr_range[i++][1] = log(mass_range[1]);

  return;
}


/*==============================================================================*/
/* the following functions are deprecated */

/*
 * setup  SA BLR model parameter range. 
 */
// model 1
void set_sa_blr_range_model1()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //lambda
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.5;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //sa_blr_range_model[2][1] = fmin(sa_blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 2
void set_sa_blr_range_model2()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //sigr
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigtheta
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  return;
}

// model 3
void set_sa_blr_range_model3()
{
  int i;
  
  i = 0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn 
  sa_blr_range_model[i][0] = 0.0;   // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[3][1] + sa_blr_range_model[4][1]));
  return;
}

// model 4
void set_sa_blr_range_model4()
{
  int i;
  
  i = 0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = sqrt(2.0)/2.0;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[3][1] + sa_blr_range_model[4][1]));
  return;
}

// model 5
void set_sa_blr_range_model5()
{
  int i;
  
  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //Fin
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //Fout
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(10.0);
  //alpha
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //beta
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;
  //sig_turb
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[2][1] + sa_blr_range_model[4][1]));

  return;
}


// model 6
void set_sa_blr_range_model6()
{
  int i;
  
  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; //in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //gamma
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;
  //sig_turb
  sa_blr_range_model[i][0] = log(0.0001);
  sa_blr_range_model[i++][1] = log(0.1);

  //rcloud_max_set = 44.85;
  //sa_blr_range_model[2][1] = fmin(sa_blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 7
void set_sa_blr_range_model7()
{
  int i;
  

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //gamma
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //fsh
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //mu_un
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta_un
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //opn_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;

  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];

  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;

  //fellip_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //sig_turb
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);

  return;
}

// model 8
void set_sa_blr_range_model8()
{
  int i;
  
  i = 0;
  //theta_min
  sa_blr_range_model[i][0] = 20.0;     // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //dtheta_max
  sa_blr_range_model[i][0] = 0.0;      // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //r_min
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.25);
  //fr_max
  sa_blr_range_model[i][0] = 1.0;  //
  sa_blr_range_model[i++][1] = 10.0;
  //gamma
  sa_blr_range_model[i][0] = 0.0;  //
  sa_blr_range_model[i++][1] = 3.0;
  //alpha
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 3.0;
  //lambda
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 0.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //Rv
  sa_blr_range_model[i][0] = log(10.0);
  sa_blr_range_model[i++][1] = log(50.0);
  //Rblr
  sa_blr_range_model[i][0] = log(10.0);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //inc
  sa_blr_range_model[i][0] = 0.0;   // cos(inc)
  sa_blr_range_model[i++][1] = 1.0;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);

  return;
}

// model 9
void set_sa_blr_range_model9()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //eta0
  sa_blr_range_model[i][0]   = eta_range_model[0][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[0][1];
  //eta1
  sa_blr_range_model[i][0]   = eta_range_model[1][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[1][1];
  //eta_alpha
  sa_blr_range_model[i][0]   = eta_range_model[2][0];  // 
  sa_blr_range_model[i++][1] = eta_range_model[2][1];
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);

  return;
}