//#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <string>
#include <mpi.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>

using std::cout;
using std::cin;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;


typedef struct
{
	double x, y, z;
	int e[8];
	int ID;
	int n_n[7];
	int corner_ID;
	int val;
	int loc;
	int proc;
	int local;
	int global;
	int req;
} MNODE;

typedef struct
{
	int e[10];
} TMP_NODE;

typedef struct
{
	int n_n[7];
} singul;

typedef struct
{
	double u, v, e;
} RESIDUAL;

typedef struct
{
	int ID;
} COR_NODE;


typedef struct
{
	double zeta_x, zeta_y, zeta_z, eta_x, eta_y, eta_z, xi_x, xi_y, xi_z;
	double zeta_xip, zeta_yip, zeta_zip, eta_xip, eta_yip, eta_zip, xi_xip, xi_yip, xi_zip, zeta_xim, zeta_yim, zeta_zim, eta_xim, eta_yim, eta_zim, xi_xim, xi_yim, xi_zim;
	double zeta_xjp, zeta_yjp, zeta_zjp, eta_xjp, eta_yjp, eta_zjp, xi_xjp, xi_yjp, xi_zjp, zeta_xjm, zeta_yjm, zeta_zjm, eta_xjm, eta_yjm, eta_zjm, xi_xjm, xi_yjm, xi_zjm;
	double zeta_xkp, zeta_ykp, zeta_zkp, eta_xkp, eta_ykp, eta_zkp, xi_xkp, xi_ykp, xi_zkp, zeta_xkm, zeta_ykm, zeta_zkm, eta_xkm, eta_ykm, eta_zkm, xi_xkm, xi_ykm, xi_zkm;
} TRANSFORMED;

typedef struct
{
	double x_zeta, x_eta, x_xi, y_zeta, y_eta, y_xi, z_zeta, z_eta, z_xi;
} JACOB;

typedef struct
{
	double ip;
} Qip;

typedef struct
{
	double ip, im, jp, jm, kp, km;
} DETERM;

typedef struct
{
	double zeta_1, zeta_2, zeta_3;
	double eta_1, eta_2, eta_3;
	double xi_1, xi_2, xi_3;
} JAC;

typedef struct
{
	double zeta_1, zeta_2, zeta_3, zeta_4, zeta_5, zeta_6;
	double eta_1, eta_2, eta_3, eta_4, eta_5, eta_6;
	double xi_1, xi_2, xi_3, xi_4, xi_5, xi_6;
} TEMP_JAC;

typedef struct
{
	int connect[9], element_neighbour[9];
	int type;
	int location;
	int corner;
} ELEM;

typedef struct
{
	int check;
} BOUND;

typedef struct
{
	int proc;
	int n_proc;
	int global;
} SUB_DOM;

typedef struct
{
	int proc;
	int node;
	int i;
	int j;
	int k;
} DOM_TIP;


extern char inputfilename[];
extern char filename[];
extern double del_zeta;
extern double del_eta;
extern double Mach;
extern double Free_t;
extern double Free_rho;
extern double Free_mu;
extern double Free_pres;
extern double Free_kine;
extern double Free_sound;
extern double univ_gas;
extern double char_length;
//extern float gamma;
extern double Pr;
extern double del_t;
extern double Epsilon;
extern int restart_num;
extern double Reyl;
extern int iterations;
extern char file_ext[];
extern double back_pressure;
extern int AIRFOIL;
extern int OTHER_GEOM;

extern int itera1, itera2, itera;
extern int slaves[];
extern int myrank, size;
extern int i, j, NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL, inl, out, wal, bou, NE, IBCODE, all_bou, bound, memory, memory_elem, nodes, no_of_nodes, vn, n[3];
extern int out_old, inl_old, wal_old, bou_old, out_mem, out_nl, inl_mem, inl_nl, wal_mem, wal_nl, bou_mem, bou_nl;
extern int inl_elem, out_elem, wal_elem, bou_elem;
extern int temp2, temp3, temp4, temp5, temp6, temp7;
extern int gar1, gar2, element, elem1, n1, count, orient[2], line1, shift, len, len2, shift2, gar3, gar4, gar5, value;
extern int P0, P1, P2, P3, P4, P5, P6, P7;
extern int g_node, g_elem, all_bou_node;
extern int gambit, icemcfd;
extern int corner_element[20], corner_node[20], cor, wal_node, initial, inl_node, out_node, corn_node, bou_node, temp, restart, position;
extern char garbage1[], filename2[];
extern double angle[], max_angle, denominator, determinant, temp_d;
extern double temp_u, temp_v, temp_w, temp_p, temp_t, temp_e, temp_mu, temp_rho;
extern int sd_gh, sd_node, h, gh, temp1, no_of_tip_send, no_of_tip_recv;
extern double Qi_iplus_half_pos[6], Qi_iminus_half_pos[6], Qi_iplus_half_neg[6], Qi_iminus_half_neg[6], Qi_iplus_half_f[6], Qi_iminus_half_f[6];
extern double Qj_iplus_half_pos[6], Qj_iminus_half_pos[6], Qj_iplus_half_neg[6], Qj_iminus_half_neg[6], Qj_iplus_half_E[6], Qj_iminus_half_E[6];
extern double Qk_iplus_half_pos[6], Qk_iminus_half_pos[6], Qk_iplus_half_neg[6], Qk_iminus_half_neg[6], Qk_iplus_half_G[6], Qk_iminus_half_G[6];
extern double Qi_iplus_half_pos_char[6], Qi_iminus_half_pos_char[6], Qi_iplus_half_neg_char[6], Qi_iminus_half_neg_char[6];
extern double Qj_iplus_half_pos_char[6], Qj_iminus_half_pos_char[6], Qj_iplus_half_neg_char[6], Qj_iminus_half_neg_char[6];
extern double Qk_iplus_half_pos_char[6], Qk_iminus_half_pos_char[6], Qk_iplus_half_neg_char[6], Qk_iminus_half_neg_char[6];
extern double h_F_ip[6], h_F_im[6], h_E_ip[6], h_E_im[6], h_G_ip[6], h_G_im[6], d2F_d2z_ip[6], d2F_d2z_im[6], d2E_d2e_ip[6], d2E_d2e_im[6], d2G_d2x_ip[6], d2G_d2x_im[6];
extern double F_ip_h[6], F_im_h[6], E_ip_h[6], E_im_h[6], G_ip_h[6], G_im_h[6], d4F_d4z_ip[6], d4F_d4z_im[6], d4E_d4e_ip[6], d4E_d4e_im[6], d4G_d4x_ip[6], d4G_d4x_im[6];
extern double F_ip_pos[6], F_ip_neg[6], F_im_pos[6], F_im_neg[6], E_ip_pos[6], E_ip_neg[6], E_im_pos[6], E_im_neg[6], G_ip_pos[6], G_ip_neg[6], G_im_pos[6], G_im_neg[6];
extern double F_jp_pos[6], F_jp_neg[6], F_jm_pos[6], F_jm_neg[6], E_jp_pos[6], E_jp_neg[6], E_jm_pos[6], E_jm_neg[6], G_jp_pos[6], G_jp_neg[6], G_jm_pos[6], G_jm_neg[6];
extern double F_kp_pos[6], F_kp_neg[6], F_km_pos[6], F_km_neg[6], E_kp_pos[6], E_kp_neg[6], E_km_pos[6], E_km_neg[6], G_kp_pos[6], G_kp_neg[6], G_km_pos[6], G_km_neg[6];
extern double F_ip_pos_comp[6], F_ip_neg_comp[6], F_im_pos_comp[6], F_im_neg_comp[6], E_ip_pos_comp[6], E_ip_neg_comp[6], E_im_pos_comp[6], E_im_neg_comp[6], G_ip_pos_comp[6], G_ip_neg_comp[6], G_im_pos_comp[6], G_im_neg_comp[6];
extern double F_jp_pos_comp[6], F_jp_neg_comp[6], F_jm_pos_comp[6], F_jm_neg_comp[6], E_jp_pos_comp[6], E_jp_neg_comp[6], E_jm_pos_comp[6], E_jm_neg_comp[6], G_jp_pos_comp[6], G_jp_neg_comp[6], G_jm_pos_comp[6], G_jm_neg_comp[6];
extern double F_kp_pos_comp[6], F_kp_neg_comp[6], F_km_pos_comp[6], F_km_neg_comp[6], E_kp_pos_comp[6], E_kp_neg_comp[6], E_km_pos_comp[6], E_km_neg_comp[6], G_kp_pos_comp[6], G_kp_neg_comp[6], G_km_pos_comp[6], G_km_neg_comp[6];
extern double F_ip[6], F_im[6], E_ip[6], E_im[6], G_ip[6], G_im[6], F_jp[6], F_jm[6], E_jp[6], E_jm[6], G_jp[6], G_jm[6], F_kp[6], F_km[6], E_kp[6], E_km[6], G_kp[6], G_km[6];
extern vector<double> alpha_u_ip, alpha_v_ip, alpha_w_ip, alpha_u_im, alpha_v_im, alpha_w_im;
extern double dF_W[6], dF_C[6], dE_W[6], dE_C[6], dG_W[6], dG_C[6];
extern vector<double> alpha_u_jp, alpha_v_jp, alpha_w_jp, alpha_u_jm, alpha_v_jm, alpha_w_jm, alpha_u_kp, alpha_v_kp, alpha_w_kp, alpha_u_km, alpha_v_km, alpha_w_km, alpha_u, alpha_v, alpha_w;
extern double eigen_Qip[8], eigen_Qim[8], eigen_Qjp[8], eigen_Qjm[8], eigen_Qkp[8], eigen_Qkm[8], eigen_Qinp[8], eigen_Qinm[8], eigen_Qjnp[8], eigen_Qjnm[8], eigen_Qknp[8], eigen_Qknm[8];
extern double zeta_xip, zeta_xim, zeta_yip, zeta_yim, zeta_zip, zeta_zim, eta_xjp, eta_xjm, eta_yjp, eta_yjm, eta_zjp, eta_zjm, xi_xkp, xi_xkm, xi_ykp, xi_ykm, xi_zkp, xi_zkm;
extern double e_inf, max_div_um, max_div_vm, max_div_wm, max_div_em;
extern vector<int> d;
extern vector<MNODE> d_node, tmp_node;
extern vector<SUB_DOM> sd_gh_node;
extern vector<int> c, tmp_c, recv_c, sd_inlet_node, sd_outlet_node, sd_wall_node, sd_boundary_node, proc_node;
extern vector<vector<int>> b, gb, loc_dat, recv_b;
extern int neigh_pro, glob, loca, sd_inl_node, sd_wal_node, sd_out_node, sd_bou_node;
extern string line, line2;
extern string line_n;
extern vector<JACOB> t_jacobian;
extern vector<TRANSFORMED> t_metric;
extern vector<DOM_TIP> tip;
extern vector<DOM_TIP> tip_recv;
extern vector<double> t_det;
extern ifstream inFile;
extern ofstream outFile;
extern vector<double> det, DUCROS;
extern vector<int> boundary_elements, all_boundary_nodes, temp_boundary_nodes;
extern vector<int> inlet_node, outlet_node, wall_node, boundary_node;
extern vector<double> tauzz, tauxx, tauzx, tauex, tauze, tauee, qz, qe, qx;
extern vector<double> roe_rho_ip, roe_u_ip, roe_v_ip, roe_w_ip, roe_h_ip, roe_rho_im, roe_u_im, roe_v_im, roe_w_im, roe_h_im, roe_rho_jp, roe_u_jp, roe_v_jp, roe_w_jp, roe_h_jp, roe_rho_jm, roe_u_jm, roe_v_jm, roe_w_jm, roe_h_jm, roe_rho_kp, roe_u_kp, roe_v_kp, roe_w_kp, roe_h_kp, roe_rho_km, roe_u_km, roe_v_km, roe_w_km, roe_h_km;
extern vector<double> Qi_iplus_half, Qi_iminus_half, Qj_iplus_half, Qj_iminus_half, Qk_iplus_half, Qk_iminus_half;
extern vector<double> dF, dE, dG, dFv, dEv, dGv, final_U, roe_R, roe_a, del_cfl, v_dash1;
extern vector<double> TAU_SGS_XY, TAU_SGS_YZ, TAU_SGS_XZ, TAU_SGS_XX, TAU_SGS_YY, TAU_SGS_ZZ, H_SGS_X, H_SGS_Y, H_SGS_Z, D_SGS_X, D_SGS_Y, D_SGS_Z;
extern vector<vector<double>> u, v, w, rho, p, t, mu, a, e, ki;
extern vector<vector<double>> Qi_iminus, Qi_iplus, Qi_iminus_n, Qi_iplus_n, Qi_half_p, Qi_half_np, Qi_half_m, Qi_halfn_m;
extern vector<vector<double>> Qj_iminus, Qj_iplus, Qj_iminus_n, Qj_iplus_n, Qj_half_p, Qj_half_np, Qj_half_m, Qj_halfn_m;
extern vector<vector<double>> Qk_iminus, Qk_iplus, Qk_iminus_n, Qk_iplus_n, Qk_half_p, Qk_half_np, Qk_half_m, Qk_halfn_m;
extern vector<vector<double>> IS_Qim, IS_Qinm, IS_Qip, IS_Qinp, IS_Qjm, IS_Qjnm, IS_Qjp, IS_Qjnp, IS_Qkm, IS_Qknm, IS_Qkp, IS_Qknp, w_Qip, w_Qinp, w_Qim, w_Qinm, w_Qjp, w_Qjnp, w_Qjm, w_Qjnm, w_Qkp, w_Qknp, w_Qkm, w_Qknm, W_Qip, W_Qinp, W_Qim, W_Qinm, W_Qjp, W_Qjnp, W_Qjm, W_Qjnm, W_Qkp, W_Qknp, W_Qkm, W_Qknm;
extern vector<vector<double>> E, F, G, Ev, Fv, Gv, F1, E1, G1, Fv1, Ev1, Gv1;
extern vector<vector<double>> r_eigen_Qip, r_eigen_Qjp, r_eigen_Qkp, r_eigen_Qim, r_eigen_Qjm, r_eigen_Qkm, l_eigen_Qip, l_eigen_Qjp, l_eigen_Qkp, l_eigen_Qim, l_eigen_Qjm, l_eigen_Qkm;
extern vector<vector<vector<double>>> U, L;
extern vector<MNODE> node;
extern vector<singul> singular;
extern vector<vector<RESIDUAL>> diver;
extern vector<JACOB> jacobian;
extern vector<ELEM> CD;
extern vector<ELEM> wq;
extern vector<TRANSFORMED> metric;
extern vector<DETERM> deter;
extern vector<vector<Qip>> Q;
extern vector<JAC> jac;
extern vector<TEMP_JAC> jaco;
extern vector<TMP_NODE> temp_node;
extern vector<BOUND> temp_check;
extern vector<int> all_neigh_pro;
extern vector<int> temp_a;
extern double temp_u, temp_v, temp_w, temp_p, temp_t, temp_e, temp_mu, temp_rho;
extern vector<double> max_div_u, max_div_v, max_div_e;
extern vector<int> sd_count, sd_gh_list;
extern int iter, temp_nodes, glob;
extern vector<int> neigh_pro_list;
extern vector<vector<int>> all_c_list;
extern int iteration;

int inputfile_preprocessor();

void calculate_mesh();

void restart_file();

void boundary_nodes_count();

void read_mesh_file();

void node_neighbour();

void node_location();

void nodefile();

void metric_file();

void subdomain_neighbour();

void proc_request();

void metric_term();

void writepltfile();

void write_restart_file();

void writepltfile_jacob();

void viscousflux_variables();

void roe_average();

void intialise(int RK);

void solver();

void weno_solver_0();

void weno_solver_1();

void weno_solver_2();

void weno_solver_3();

void weno_solver_4();

void weno_solver_5();


