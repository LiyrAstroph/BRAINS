/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file proto.h
 *  \brief declare proto types.
 */

void begin_run();
void end_run();
void read_parset();

void read_data();
void allocate_memory_data();
void free_memory_data();
void allocate_memory();
void free_memory();


void init();
void scale_con_line();
void cal_emission_flux();
void get_num_particles(char *fname);
void get_posterior_sample_file(char *fname, char *samplefile);
void set_par_fix(int num_params_blr);

void print_help();

void (*set_blr_range_model)();
void set_blr_range_model1();
void set_blr_range_model2();
void set_blr_range_model3();
void set_blr_range_model4();
void set_blr_range_model5();

/* continuum reconstruction */
int dnest_con(int argc, char **argv);
void reconstruct_con_init();
void reconstruct_con();
void reconstruct_con_end();
double prob_con_variability(const void *model);
double prob_con_variability_initial(const void *model);
void set_covar_Pmat(double sigma, double tau, double alpha);
void set_covar_Pmat_data(double sigma, double tau, double alpha,double syserr);
void set_covar_Umat(double sigma, double tau, double alpha);
void calculate_con_from_model(const void *model);
void reconstruct_con_from_varmodel(double sigma_hat, double tau, double alpha, double syserr);
void create_con_from_random(double sigma_hat, double tau, double alpha, double syserr);
void postprocess_con();
void set_par_range_con();

/* 1d line reconstruction */
int dnest_line1d(int argc, char **argv);
void reconstruct_line1d();
void reconstruct_line1d_init();
void reconstruct_line1d_end();
double prob_line1d(const void *model);
double prob_initial_line1d(const void *model);
double prob_restart_line1d(const void *model);

void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl);
void (*transfun_1d_cloud_direct)(const void *pm, int flag_save);

void transfun_1d_cloud_direct_model1(const void *pm, int flag_save);
void transfun_1d_cloud_direct_model3(const void *pm, int flag_save);
void transfun_1d_cloud_direct_model5(const void *pm, int flag_save);

void postprocess1d();
void set_par_range_model1d();

void restart_clouds_1d(int iflag);

/* 2d line reconstruction */
int dnest_line2d(int argc, char **argv);
double prob_line2d(const void *model);
double prob_initial_line2d(const void *model);
double prob_restart_line2d(const void *model);
void reconstruct_line2d();
void reconstruct_line2d_init();
void reconstruct_line2d_end();
void line_gaussian_smooth_2D_FFT(const double *transv, double *fl2d, int nl, int nv);
void (*transfun_2d_cloud_direct)(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv);
void postprocess2d();
void set_par_range_model2d();

void transfun_2d_cloud_direct_model1(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cloud_direct_model2(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cloud_direct_model3(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cloud_direct_model4(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cloud_direct_model5(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);

void smooth_init(int nv, const double *transv);
void smooth_end();

void restart_clouds_2d(int iflag);

/* sim */
void sim();
void sim_init();
void sim_end();

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);

/* matrix operations */
void inverse_mat(double *a, int n, int *info);
double det_mat(double *a, int n, int *info);
double lndet_mat(double *a, int n, int *info);
void display_mat(double *a, int m, int n);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k);
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
void Chol_decomp_L(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);
void test_mathfun();