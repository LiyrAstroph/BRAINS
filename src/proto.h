/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
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

/* continuum reconstruction */
int dnest_con(int argc, char **argv);
void reconstruct_con_init();
void reconstruct_con();
void reconstruct_con_end();
double prob_con_variability(const void *model);
void set_covar_Pmat(double sigma, double tau, double alpha);
void set_covar_Pmat_data(double sigma, double tau, double alpha);
void set_covar_Umat(double sigma, double tau, double alpha);
void calculate_con_from_model(const void *model);
void reconstruct_con_from_varmodel(double sigma, double tau, double alpha);

/* 1d line reconstruction */
int dnest_line1d(int argc, char **argv);
void reconstruct_line1d();
void reconstruct_line1d_init();
void reconstruct_line1d_end();
double prob_line1d(const void *model);
void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl);
void transfun_1d_cloud_direct(const void *pm);


/* 2d line reconstruction */
int dnest_line2d(int argc, char **argv);
double prob_line2d(const void *model);
void reconstruct_line2d();
void reconstruct_line2d_init();
void reconstruct_line2d_end();
void line_gaussian_smooth_2D_FFT(const double *transv, double *fl2d, int nl, int nv);
void transfun_2d_cloud_direct(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv);

void smooth_init(int nv);
void smooth_end();

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
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
void Chol_decomp_L(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);