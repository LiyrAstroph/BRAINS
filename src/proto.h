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

void print_help();
void fprint_param();
int search_pardict(char *tag);
int check_parset_isset();
int check_time_sorted(double *time_series, int n);

void (*set_blr_range_model)();
void set_blr_range_model1();
void set_blr_range_model2();
void set_blr_range_model3();
void set_blr_range_model4();
void set_blr_range_model5();
void set_blr_range_model6();
void set_blr_range_model7();
void set_blr_range_model8();
void set_blr_range_model9();

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

int get_idx_mbh_from_blrmodel();

/* sim */
void sim();
void sim_init();
void sim_end();
void set_par_value_sim(double *pm, int flag_model);

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);

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

void (*set_sa_blr_range_model)();
void set_sa_blr_range_model1();
void set_sa_blr_range_model2();
void set_sa_blr_range_model3();
void set_sa_blr_range_model4();
void set_sa_blr_range_model5();
void set_sa_blr_range_model6();
void set_sa_blr_range_model7();
void set_sa_blr_range_model8();
void set_sa_blr_range_model9();

/* SA smooth */
void sa_smooth_init(int n_v_sa, const double *v_sa, double sigV);
void sa_smooth_end();
void sa_smooth_run(double *v_sa, double *F_sa, int n_v_sa, double *p_sa, int n_base_sa);
#endif