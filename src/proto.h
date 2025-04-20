/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
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

double get_mediate_cad(double *tcon, int ncon);

void init();
void scale_con_line();
void cal_emission_flux();
void get_num_particles(char *fname);
void set_par_fix_blrmodel();
void load_par_names(char *fname);
void set_drw_par_range();
void set_drw_par_range_load();

void print_help();
void fprint_param();
int search_pardict(char *tag);
int check_parset_isset();
int check_time_sorted(double *time_series, int n);
int check_equal_bin(double *x, int n);

void (*set_blr_range_model)(double **blr_range);
void set_blr_range_model1(double **blr_range);
void set_blr_range_model2(double **blr_range);
void set_blr_range_model3(double **blr_range);
void set_blr_range_model4(double **blr_range);
void set_blr_range_model5(double **blr_range);
void set_blr_range_model6(double **blr_range);
void set_blr_range_model7(double **blr_range);
void set_blr_range_model8(double **blr_range);
void set_blr_range_model9(double **blr_range);
void set_nlr_range_model();

/* continuum reconstruction */
int dnest_con(int argc, char **argv);
void reconstruct_con_init();
void reconstruct_con();
void reconstruct_con_end();
double prob_con_variability(const void *model);
double prob_con_variability_initial(const void *model);
double prob_con_variability_semiseparable(const void *model);
double prob_con_variability_initial_semiseparable(const void *model);
void set_covar_Pmat(double sigma, double tau, double alpha);
void set_covar_Pmat_data(double sigma, double tau, double alpha,double syserr);
void set_covar_Umat(double sigma, double tau, double alpha);
void calculate_con_from_model(const void *model);
void calculate_con_from_model_semiseparable(const void *model);
void reconstruct_con_from_varmodel(double sigma_hat, double tau, double alpha, double syserr);
void create_con_from_random(double sigma_hat, double tau, double alpha, double syserr);
void postprocess_con();
void set_par_range_con();
void print_par_names_con();
void calculate_con_rm(const void *pm);
double interp_con_rm(double tp);

void set_blr_model1d();
void set_blr_model2d();

/* 1d line reconstruction */
int dnest_line1d(int argc, char **argv);
void reconstruct_line1d();
void reconstruct_line1d_init();
void reconstruct_line1d_end();
double prob_line1d(const void *model);
double prob_initial_line1d(const void *model);
double prob_restart_line1d(const void *model);

void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl);
void (*transfun_1d_cal)(const void *pm, int flag_save);
void transfun_1d_cal_cloud(const void *pm, int flag_save);
void transfun_1d_cal_with_sample();

void postprocess1d();
void set_par_range_model1d();
void print_par_names_model1d();

/* 2d line reconstruction */
int dnest_line2d(int argc, char **argv);
double prob_line2d(const void *model);
double prob_initial_line2d(const void *model);
double prob_restart_line2d(const void *model);
void reconstruct_line2d();
void reconstruct_line2d_init();
void reconstruct_line2d_end();
void line_gaussian_smooth_2D_FFT(const double *transv, double *fl2d, int nl, int nv, const void *pm);
void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv);
void postprocess2d();
void set_par_range_model2d();
void print_par_names_model2d();

void (*gen_cloud_sample)(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model1(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model2(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model3(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model4(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model5(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model6(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model7(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model8(const void *pm, int flag_type, int flag_save);
void gen_cloud_sample_model9(const void *pm, int flag_type, int flag_save);

void (*transfun_2d_cal)(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cal_cloud(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void transfun_2d_cal_with_sample(double *transv, double *trans2d, int n_vel);

void smooth_init(int nv, const double *transv);
void smooth_end();
void smooth_test();

int get_idx_mbh_from_blrmodel(int blrmodel);
int get_idx_blrsize_from_blrmodel(int blrmodel);

/* line profile fitting */
void postprocesslp();
void reconstruct_lp();
void reconstruct_lp_init();
void reconstruct_lp_end();
void postprocesslp();

int dnest_lp(int argc, char **argv);
double prob_lp(const void *model);
void set_par_range_lp();
void print_par_names_lp();
void line_gaussian_smooth_FFT(const double *transv, double *fl, int nv, const void *pm);
void cal_line_profile_with_sample(const void *pm, double *vel, double *fv, int nvel);


/* sim */
void sim();
void sim_init();
void sim_end();
void set_par_value_sim(double *pm, int flag_model);
void print_par_value_sim(double *pm, int np);
void read_param_value(double *pm, char *fname);

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);

void smooth_transfer_function_tau(double *trans1d);
void smooth_transfer_function2d_tau(double *trans2d, int n_vel);

#ifdef SpecAstro

void set_sa_blr_model();
void set_idx_par_mutual();
void set_par_fix_sa_blrmodel();

/* SA */
void reconstruct_sa();
void reconstruct_sa_init();
void reconstruct_sa_end();

int dnest_sa(int argc, char **argv);
double prob_sa(const void *model);
void set_par_range_sa();
void print_par_names_sa();
void set_par_range_sa1d();
void print_par_names_sa1d();
void set_par_range_sa2d();
void print_par_names_sa2d();
void calculate_sa_from_blrmodel(const void *pm, int flag_save);
void postprocess_sa();

/* SA + 1D RM */
void reconstruct_sa1d();
void reconstruct_sa1d_init();
void reconstruct_sa1d_end();

int dnest_sa1d(int argc, char **argv);
double prob_sa1d(const void *model);
double prob_initial_sa1d(const void *model);
void set_par_range_sa1d();
void print_par_names_sa1d();
void set_par_range_sa1d();
void print_par_names_sa1d();
void calculate_sa1d_from_blrmodel(const void *pm);
void calculate_sa_with_sample(const void *pm);
void calculate_sa_sim_with_sample(const void *pm, double *vel_sa, int n_vel_sa, double *base_sa, int n_base_sa,
                                  double *p_sa, double *F_sa);
void postprocess_sa1d();

void calculate_sa_transfun_from_blrmodel(const void *pm, int flag_save);

/* SA + 2D RM */
void reconstruct_sa2d();
void reconstruct_sa2d_init();
void reconstruct_sa2d_end();

int dnest_sa2d(int argc, char **argv);
double prob_sa2d(const void *model);
double prob_initial_sa2d(const void *model);
void set_par_range_sa2d();
void print_par_names_sa2d();
void set_par_range_sa2d();
void print_par_names_sa2d();
void calculate_sa2d_from_blrmodel(const void *pm);
void postprocess_sa2d();

void calculate_sa_transfun2d_from_blrmodel(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);

/* SA general */
void (*gen_sa_cloud_sample)(const void *pm, int flag_type, int flag_save);

void (*set_sa_blr_range_model)(double **blr_range);
// void set_sa_blr_range_model1();
// void set_sa_blr_range_model2();
// void set_sa_blr_range_model3();
// void set_sa_blr_range_model4();
// void set_sa_blr_range_model5();
// void set_sa_blr_range_model6();
// void set_sa_blr_range_model7();
// void set_sa_blr_range_model8();
// void set_sa_blr_range_model9();

/* SA smooth */
void sa_smooth_init(int n_v_sa, const double *v_sa, double sigV);
void sa_smooth_end();
void sa_smooth_run(double *v_sa, double *F_sa, int n_v_sa, double *p_sa, int n_base_sa);

/* SARM */
int dnest_sarm(int argc, char **argv);
double prob_sarm(const void *model);
double prob_initial_sarm(const void *model);
void set_par_range_sarm();
void print_par_names_sarm();
void reconstruct_sarm();
void reconstruct_sarm_init();
void reconstruct_sarm_end();
void cal_emission_flux_sarm();
void scale_con_line_sarm();
void postprocess_sarm();
void transfun_sarm_cal_cloud(const void *pm, double *transv, double *trans2d, double *trans_sarm_alpha, double *trans_sarm_beta, 
                             int n_vel, int flag_save);
void transfun_sarm_cal_with_sample(double *transv, double *trans2d, double *trans_sarm_alpha, double *trans_sarm_beta, int n_vel);

void calculate_sarm_with_sample(const void *pm);
void calculate_sarm_sim_with_sample(const void *pm, double *tline_sarm, double *vel_sa, double *trans2d, 
                                    double *trans_alpha, double *trans_beta,
                                    int n_sa_vel, int n_sarm_line, double *base_sarm, 
                                    int n_sarm_base, double *phase_sarm, double *Fline_sarm, 
                                    double *momentum_alpha, double *momentum_beta,
                                    double *photocenter_alpha, double *photocenter_beta);

void sarm_smooth_init(int n_v_sarm, const double *v_sarm, double sigV);
void sarm_smooth_end();
void sarm_smooth_run(const void *pm, double *v_sarm, double *F_sarm, int n_v_sarm, int n_line_sarm, double *p_sarm, int n_base_sarm);
#endif