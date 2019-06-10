/**********************************************SRI ANJINEYA***********************************************************/
/****************************************OM GAM GAM GANAPATHAYE NAMAHA************************************************/
/****************************************************OM***************************************************************/
/**********************************************************************************************************************/
/**********************************************Developed By************************************************************/
/*********************************G SRINIVASAN AND PROF K.P.SINHAMAHAPATRA*********************************************/
/**********************************INDIAN INSTITUTE OF TECHNOLOGY KHARAGPUR********************************************/

//#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <string>
#include "functions.h"
#include <mpi.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using std::cout;
using std::cin;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;

char inputfilename[] = "wedge";
char filename[] = "wedge.neu";
double Epsilon = 10e-8;
double del_zeta = 1.0;
double del_eta = 1.0;
double Mach = 1.5;
double Free_t = 273.0;
double Free_rho = 1.225;
double Free_mu = 1.7894e-5;
double Free_pres = 101325.0;
double Free_kine = 1.4;
double Free_sound = 340.28;
double univ_gas = 287.0;
double char_length = 1.0;
double Pr = 0.72;
double del_t = 0.005;
int restart_num = 1;
double Reyl = 500.0;
int iterations = 20000000;
char file_ext[] = "cfx5";             /**********neu if file made in gambit******cfx5 if file made in icemcfd***********/
double back_pressure = 0.0;
int AIRFOIL_GEOM = 0;
int OTHER_GEOM = 0;

/*********************************Global variables***********************************************************************/
int IMAX, JMAX, KMAX;
vector<vector<vector<int>>> ordered_node;
int itera1, itera2, itera;
int slaves[1];
int myrank, size;
int i, j, NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL, inl, out, wal, bou, NE, IBCODE, all_bou, bound, memory, memory_elem, nodes, no_of_nodes, vn, n[3];
int out_old, inl_old, wal_old, bou_old, out_mem, out_nl, inl_mem, inl_nl, wal_mem, wal_nl, bou_mem, bou_nl;
int inl_elem, out_elem, wal_elem, bou_elem;
int temp1, temp2, temp3, temp4, temp5, temp6, temp7;
int gar1, gar2, element, elem1, n1, count, orient[2], line1, shift, len, len2, shift2, gar3, gar4, gar5, value;
int P0, P1, P2, P3, P4, P5, P6, P7;
int g_node, g_elem, all_bou_node;
int gambit, icemcfd;
int corner_element[20], corner_node[20], cor, wal_node, initial, inl_node, out_node, corn_node, bou_node, temp, restart, position;
int sd_gh, sd_node, h, gh, no_of_tip_send, no_of_tip_recv;

char garbage1[50], filename2[50];

double angle[4], max_angle, denominator, determinant, temp_d;
double temp_u, temp_v, temp_w, temp_p, temp_t, temp_e, temp_mu, temp_rho;
double Qi_iplus_half_pos[6], Qi_iminus_half_pos[6], Qi_iplus_half_neg[6], Qi_iminus_half_neg[6], Qi_iplus_half_f[6], Qi_iminus_half_f[6];
double Qj_iplus_half_pos[6], Qj_iminus_half_pos[6], Qj_iplus_half_neg[6], Qj_iminus_half_neg[6], Qj_iplus_half_E[6], Qj_iminus_half_E[6];
double Qk_iplus_half_pos[6], Qk_iminus_half_pos[6], Qk_iplus_half_neg[6], Qk_iminus_half_neg[6], Qk_iplus_half_G[6], Qk_iminus_half_G[6];
double Qi_iplus_half_pos_char[6], Qi_iminus_half_pos_char[6], Qi_iplus_half_neg_char[6], Qi_iminus_half_neg_char[6];
double Qj_iplus_half_pos_char[6], Qj_iminus_half_pos_char[6], Qj_iplus_half_neg_char[6], Qj_iminus_half_neg_char[6];
double Qk_iplus_half_pos_char[6], Qk_iminus_half_pos_char[6], Qk_iplus_half_neg_char[6], Qk_iminus_half_neg_char[6];
double h_F_ip[6], h_F_im[6], h_E_ip[6], h_E_im[6], h_G_ip[6], h_G_im[6], d2F_d2z_ip[6], d2F_d2z_im[6], d2E_d2e_ip[6], d2E_d2e_im[6], d2G_d2x_ip[6], d2G_d2x_im[6];
double F_ip_h[6], F_im_h[6], E_ip_h[6], E_im_h[6], G_ip_h[6], G_im_h[6], d4F_d4z_ip[6], d4F_d4z_im[6], d4E_d4e_ip[6], d4E_d4e_im[6], d4G_d4x_ip[6], d4G_d4x_im[6];
double F_ip_pos[6], F_ip_neg[6], F_im_pos[6], F_im_neg[6], E_ip_pos[6], E_ip_neg[6], E_im_pos[6], E_im_neg[6], G_ip_pos[6], G_ip_neg[6], G_im_pos[6], G_im_neg[6];
double F_jp_pos[6], F_jp_neg[6], F_jm_pos[6], F_jm_neg[6], E_jp_pos[6], E_jp_neg[6], E_jm_pos[6], E_jm_neg[6], G_jp_pos[6], G_jp_neg[6], G_jm_pos[6], G_jm_neg[6];
double F_kp_pos[6], F_kp_neg[6], F_km_pos[6], F_km_neg[6], E_kp_pos[6], E_kp_neg[6], E_km_pos[6], E_km_neg[6], G_kp_pos[6], G_kp_neg[6], G_km_pos[6], G_km_neg[6];
double F_ip_pos_comp[6], F_ip_neg_comp[6], F_im_pos_comp[6], F_im_neg_comp[6], E_ip_pos_comp[6], E_ip_neg_comp[6], E_im_pos_comp[6], E_im_neg_comp[6], G_ip_pos_comp[6], G_ip_neg_comp[6], G_im_pos_comp[6], G_im_neg_comp[6];
double F_jp_pos_comp[6], F_jp_neg_comp[6], F_jm_pos_comp[6], F_jm_neg_comp[6], E_jp_pos_comp[6], E_jp_neg_comp[6], E_jm_pos_comp[6], E_jm_neg_comp[6], G_jp_pos_comp[6], G_jp_neg_comp[6], G_jm_pos_comp[6], G_jm_neg_comp[6];
double F_kp_pos_comp[6], F_kp_neg_comp[6], F_km_pos_comp[6], F_km_neg_comp[6], E_kp_pos_comp[6], E_kp_neg_comp[6], E_km_pos_comp[6], E_km_neg_comp[6], G_kp_pos_comp[6], G_kp_neg_comp[6], G_km_pos_comp[6], G_km_neg_comp[6];
double F_ip[6], F_im[6], E_ip[6], E_im[6], G_ip[6], G_im[6], F_jp[6], F_jm[6], E_jp[6], E_jm[6], G_jp[6], G_jm[6], F_kp[6], F_km[6], E_kp[6], E_km[6], G_kp[6], G_km[6];
double dF_W[6], dF_C[6], dE_W[6], dE_C[6], dG_W[6], dG_C[6];
vector<double> alpha_u_jp, alpha_v_jp, alpha_w_jp, alpha_u_jm, alpha_v_jm, alpha_w_jm, alpha_u_kp, alpha_v_kp, alpha_w_kp, alpha_u_km, alpha_v_km, alpha_w_km, alpha_u, alpha_v, alpha_w;
double eigen_Qip[8], eigen_Qim[8], eigen_Qjp[8], eigen_Qjm[8], eigen_Qkp[8], eigen_Qkm[8], eigen_Qinp[8], eigen_Qinm[8], eigen_Qjnp[8], eigen_Qjnm[8], eigen_Qknp[8], eigen_Qknm[8];
double zeta_xip, zeta_xim, zeta_yip, zeta_yim, zeta_zip, zeta_zim, eta_xjp, eta_xjm, eta_yjp, eta_yjm, eta_zjp, eta_zjm, xi_xkp, xi_xkm, xi_ykp, xi_ykm, xi_zkp, xi_zkm;
double e_inf, max_div_um, max_div_vm, max_div_wm, max_div_em;

vector<double> alpha_u_ip, alpha_v_ip, alpha_w_ip, alpha_u_im, alpha_v_im, alpha_w_im;
vector<int> d;
vector<MNODE> d_node, tmp_node;
vector<SUB_DOM> sd_gh_node;
vector<int> c, tmp_c, recv_c, sd_inlet_node, sd_outlet_node, sd_wall_node, sd_boundary_node, proc_node;
vector<vector<int>> b, gb, loc_dat, recv_b;
int neigh_pro, glob, loca, sd_inl_node, sd_wal_node, sd_out_node, sd_bou_node;
string line, line2;
string line_n;
vector<JACOB> t_jacobian;
vector<TRANSFORMED> t_metric;
vector<DOM_TIP> tip;
vector<DOM_TIP> tip_recv;
vector<double> t_det;
ifstream inFile;
ofstream outFile;
vector<double> det, DUCROS;
vector<int> boundary_elements, all_boundary_nodes, temp_boundary_nodes;
vector<int> inlet_node, outlet_node, wall_node, boundary_node;
vector<double> tauzz, tauxx, tauzx, tauex, tauze, tauee, qz, qe, qx;
vector<double> roe_rho_ip, roe_u_ip, roe_v_ip, roe_w_ip, roe_h_ip, roe_rho_im, roe_u_im, roe_v_im, roe_w_im, roe_h_im, roe_rho_jp, roe_u_jp, roe_v_jp, roe_w_jp, roe_h_jp, roe_rho_jm, roe_u_jm, roe_v_jm, roe_w_jm, roe_h_jm, roe_rho_kp, roe_u_kp, roe_v_kp, roe_w_kp, roe_h_kp, roe_rho_km, roe_u_km, roe_v_km, roe_w_km, roe_h_km;
vector<double> Qi_iplus_half, Qi_iminus_half, Qj_iplus_half, Qj_iminus_half, Qk_iplus_half, Qk_iminus_half;
vector<double> dF, dE, dG, dFv, dEv, dGv, final_U, roe_R, roe_a, del_cfl, v_dash1;
vector<double> TAU_SGS_XY, TAU_SGS_YZ, TAU_SGS_XZ, TAU_SGS_XX, TAU_SGS_YY, TAU_SGS_ZZ, H_SGS_X, H_SGS_Y, H_SGS_Z, D_SGS_X, D_SGS_Y, D_SGS_Z;
vector<vector<double>> u, v, w, rho, p, t, mu, a, e, ki;
vector<vector<double>> Qi_iminus, Qi_iplus, Qi_iminus_n, Qi_iplus_n, Qi_half_p, Qi_half_np, Qi_half_m, Qi_halfn_m;
vector<vector<double>> Qj_iminus, Qj_iplus, Qj_iminus_n, Qj_iplus_n, Qj_half_p, Qj_half_np, Qj_half_m, Qj_halfn_m;
vector<vector<double>> Qk_iminus, Qk_iplus, Qk_iminus_n, Qk_iplus_n, Qk_half_p, Qk_half_np, Qk_half_m, Qk_halfn_m;
vector<vector<double>> IS_Qim, IS_Qinm, IS_Qip, IS_Qinp, IS_Qjm, IS_Qjnm, IS_Qjp, IS_Qjnp, IS_Qkm, IS_Qknm, IS_Qkp, IS_Qknp, w_Qip, w_Qinp, w_Qim, w_Qinm, w_Qjp, w_Qjnp, w_Qjm, w_Qjnm, w_Qkp, w_Qknp, w_Qkm, w_Qknm, W_Qip, W_Qinp, W_Qim, W_Qinm, W_Qjp, W_Qjnp, W_Qjm, W_Qjnm, W_Qkp, W_Qknp, W_Qkm, W_Qknm;
vector<vector<double>> E, F, G, Ev, Fv, Gv, F1, E1, G1, Fv1, Ev1, Gv1;
vector<vector<double>> r_eigen_Qip, r_eigen_Qjp, r_eigen_Qkp, r_eigen_Qim, r_eigen_Qjm, r_eigen_Qkm, l_eigen_Qip, l_eigen_Qjp, l_eigen_Qkp, l_eigen_Qim, l_eigen_Qjm, l_eigen_Qkm;
vector<vector<vector<double>>> U, L;
vector<MNODE> node;
vector<singul> singular;
vector<vector<RESIDUAL>> diver;
vector<JACOB> jacobian;
vector<ELEM> CD;
vector<ELEM> wq;
vector<TRANSFORMED> metric;
vector<DETERM> deter;
vector<vector<Qip>> Q;
vector<JAC> jac;
vector<TEMP_JAC> jaco;
vector<TMP_NODE> temp_node;
vector<BOUND> temp_check;
vector<int> all_neigh_pro;
vector<int> temp_a;
vector<double> max_div_u, max_div_v, max_div_e;
vector<int> sd_count, sd_gh_list;
int iter, temp_nodes;
vector<int> neigh_pro_list;
vector<vector<int>> all_c_list;
int iteration;
/*********************************************************************************************************************/


/*********************************************************************************************************************/

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int k, memsize, memsize1;
	char *buffer, *buffer2;
	char *Bcast_buffer1 = NULL;
	char *Bcast_buffer2 = NULL;
	char *Bcast_buffer3 = NULL;
	char *Bcast_buffer4 = NULL;
	slaves[0] = 0;
	MPI_Group group_a, group_b;
	MPI_Comm comm_slaves;

	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	MPI_Comm_group(MPI_COMM_WORLD, &group_a);
	MPI_Group_excl(group_a, 1, slaves, &group_b);
	MPI_Comm_create(MPI_COMM_WORLD, group_b, &comm_slaves);

	if (myrank == 0)
	{
	//	inputfile_preprocessor();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int m, gar;

	d.resize(size + 1);
	initial = 0;

	/************************************************Restart File Num****************************************************/
	
	calculate_mesh();
	
	//AIRFOIL_calculate_mesh();

	restart_file();

	if (myrank == 0)
	{
		cout << "restart_num = " << restart_num << "\n";
		int restart;
		/****************************************element file************************************************************/
		ifstream fp1, fp2;
		i = 0;
		system("rm -rf node_neighbour.neu elem_neighbour.neu");
		/****************************************************************************************************************/
		
		if (strcmp(file_ext, "neu") == 0)
		{
			gambit = 1;
			icemcfd = 0;
		}

		if (strcmp(file_ext, "cfx5") == 0)
		{
			gambit = 0;
			icemcfd = 1;
		}

		cout << "Number of nodal points =  " << NUMNP << "\nNumber of Elements = " << NELEM << "\n";

		/************Obtaining number of boundary nodes*******************************************/
		
		boundary_nodes_count();

		nodes = no_of_nodes;
		/*****************************************************************************************/

		node.resize((no_of_nodes * 6) + (NUMNP + 10));
		temp_node.resize((no_of_nodes * 6) + (NUMNP + 10));
		CD.resize((no_of_nodes * 6) + (NELEM + 10));
		nodes = ((no_of_nodes * 15) + (NUMNP + 1));
		
		if (AIRFOIL_GEOM == 1)
		{
			ordered_node.resize(IMAX+10);
			for (i =1; i<= IMAX+9; i++)
			{
				ordered_node[i].resize(JMAX+10);
			}
			for (i =1; i<= IMAX+9; i++)
			{
				for (j =1; j<= JMAX+9; j++)
				{
					ordered_node[i][j].resize(KMAX+10);
				}
			}
		}
		
		for (i = 0; i <= NUMNP; i++)
		{
			node[i].x = 0.0;
			node[i].y = 0.0;
			node[i].z = 0.0;

			node[i].e[0] = 0;
			node[i].e[1] = 0;
			node[i].e[2] = 0;
			node[i].e[3] = 0;
			node[i].e[4] = 0;
			node[i].e[5] = 0;
			node[i].e[6] = 0;
			node[i].e[7] = 0;

			temp_node[i].e[0] = 0;
			temp_node[i].e[1] = 0;
			temp_node[i].e[2] = 0;
			temp_node[i].e[3] = 0;
			temp_node[i].e[4] = 0;
			temp_node[i].e[5] = 0;
			temp_node[i].e[6] = 0;
			temp_node[i].e[7] = 0;
			temp_node[i].e[8] = 0;

			node[i].n_n[0] = 0;
			node[i].n_n[1] = 0;
			node[i].n_n[2] = 0;
			node[i].n_n[3] = 0;
			node[i].n_n[4] = 0;
			node[i].n_n[5] = 0;

			node[i].ID = 0;
			node[i].corner_ID = 0;
			node[i].val = 0;
			node[i].loc = 0;
			node[i].proc = 0;
			node[i].local = 0;
			node[i].global = 0;
			node[i].req = 0;
		}

		for (i = 0; i <= NELEM; i++)
		{
			CD[i].connect[0] = 0;
			CD[i].connect[1] = 0;
			CD[i].connect[2] = 0;
			CD[i].connect[3] = 0;
			CD[i].connect[4] = 0;
			CD[i].connect[5] = 0;
			CD[i].connect[6] = 0;
			CD[i].connect[7] = 0;
		}
		/***********************Reading x and y coordinates of each node from mesh file*****************************/

		read_mesh_file();
		
		corn_node = cor;
		wal_node = wal;
		inl_node = inl;
		out_node = out;
		bou_node = bou;
		/*************************************************************************************************************************/

		/*********************************************Looking for node_NEIGHBOUR.neu FILE******************************************/
		/**************Rearranging node numbers of each element & Finding Elements connected to each node**************************/
		cout << "****************************************************\n";
				
		node_neighbour();
		/**********************************************************************************************************/			
		temp_node.clear();
		vector<TMP_NODE>().swap(temp_node);
		
		node_location();
		
		nodefile();
		/**********************************************************************************************************/
		all_bou_node = inl + wal + bou + out + corn_node;
		/**********************************************Communication***********************************************/
		
		
		position = 0;
		MPI_Pack_size((NUMNP + 10) * 3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		Bcast_buffer1 = new char [memsize];
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Pack(&node[i].x, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].y, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].z, 1, MPI_DOUBLE, Bcast_buffer1, memsize, &position, MPI_COMM_WORLD);

		}

		position = 0;
		MPI_Pack_size((NUMNP + 10) * 10, MPI_INT, MPI_COMM_WORLD, &memsize);
		Bcast_buffer2 = new char[memsize];
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Pack(&node[i].n_n[0], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[1], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[2], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[3], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[4], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].n_n[5], 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].loc, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].ID, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].proc, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&node[i].corner_ID, 1, MPI_INT, Bcast_buffer2, memsize, &position, MPI_COMM_WORLD);
		}
		printf("done broadcasting\n");

		/*********************************************************************************************************/
		ofstream fnode;
		ifstream inFile;
		ofstream outFile;
		fnode.open("metis_input.txt");

		fnode << NELEM << " 3\n";
		for (i = 1; i <= NELEM; i++)
		{
			fnode << CD[i].connect[4] << "\t" << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[7] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\n";
		}
		fnode.close();

		cout << "started domain decomposition using metis\n";
		sprintf(filename2, "nohup ./partdmesh metis_input.txt %d > metis.out", size - 1);
		system(filename2);
		for (i=1; i<=250; i++)
		{
			if (i%50 != 0)
			{
				printf("*");
			}
			if (i%50 == 0)
			{
				printf("\n");
			}
		}
		sprintf(filename2, "mv metis_input.txt.npart.%d nodes_proc.txt", size - 1);
		system(filename2);
		printf("completed decomposition using metis\n");


		inFile.open("nodes_proc.txt");
		for (i = 1; i <= NUMNP; i++)
		{
			getline(inFile, line2 , '\n');
			istringstream streamA(line2);
			streamA >> node[i].proc;
			node[i].proc = node[i].proc + 1;
		}
		inFile.close();

		for (i = 1; i<size; i++)
		{
			sprintf(filename2, "nodes_proc_%d.txt", i);
			fnode.open(filename2);
			for (j = 1; j <= NUMNP; j++)
			{
				fnode << node[j].proc << "\n";
			}
			fnode.close();
		}

		/*********************************************Writing nodefile************************************************/
		nodefile();
		/*************************************************************************************************************/

		printf("\n\n\n");

		/**********************************************************************************************
		storing all boundary nodes in a single array....
		this array will be freed before solver is executed
		***********************************************************************************************/
		
		bound = 0;
		all_bou_node = inl + out + wal + bou;
		all_boundary_nodes.resize(all_bou_node + 10);
		temp_check.resize(NUMNP + 10);

		for (i = 0; i<inl; i++)
		{
			if (node[inlet_node[i]].corner_ID == 0 && temp_check[inlet_node[i]].check == 0)
			{
				all_boundary_nodes[bound] = inlet_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i = 0; i<out; i++)
		{
			if (node[outlet_node[i]].corner_ID == 0 && temp_check[outlet_node[i]].check == 0)
			{
				all_boundary_nodes[bound] = outlet_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i = 0; i<wal; i++)
		{
			if (node[wall_node[i]].corner_ID == 0 && temp_check[wall_node[i]].check == 0)
			{
				all_boundary_nodes[bound] = wall_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i = 0; i<bou; i++)
		{
			if (node[boundary_node[i]].corner_ID == 0 && temp_check[boundary_node[i]].check == 0)
			{
				all_boundary_nodes[bound] = boundary_node[i];
				temp_check[all_boundary_nodes[bound]].check = 1;
				bound++;
			}
		}
		for (i = 0; i < corn_node; i++)
		{
			if (temp_check[corner_node[i]].check == 0)
			{
				all_boundary_nodes[bound] = corner_node[i];
				bound++;
			}
		}
		all_bou_node = bound;
		
		/***********************jacobian and Metrics calculation**************************************/
		//all_bou=bound;
		memory = (NUMNP + 1) + (3 * all_bou_node);
		memory_elem = ((NELEM + 1) + (3 * all_bou_node));
		metric.resize(memory);
		jacobian.resize(memory);
		deter.resize(memory);
		jac.resize(memory);
		jaco.resize(memory);
		det.resize(memory);

		/***************************************************************************************************/
		metric_term();

		jac.clear();
		jaco.clear();
		temp_check.clear();

		metric_file();

		/***************SEGREGATING AND WRITING BOUNDARY NODES FOR CORRESPONDING PROCESSOR******************/

		for (i = 1; i<size; i++)
		{
			inl = 0;
			out = 0;
			bou = 0;
			wal = 0;
			for (j = 0; j<inl_node; j++)
			{
				if (node[inlet_node[j]].proc == i)
				{
					inl++;
				}
			}

			sprintf(filename2, "inlet_%d.txt", i);
			outFile.open(filename2);
			for (j = 0; j<inl_node; j++)
			{
				if (node[inlet_node[j]].proc == i)
				{
					outFile << inlet_node[j] << "\n";
				}
			}
			outFile.close();

			for (j = 0; j<out_node; j++)
			{
				if (node[outlet_node[j]].proc == i)
				{
					out++;
				}
			}
			sprintf(filename2, "outlet_%d.txt", i);
			outFile.open(filename2);
			for (j = 0; j<out_node; j++)
			{
				if (node[outlet_node[j]].proc == i)
				{
					outFile << outlet_node[j] << "\n";
				}
			}
			outFile.close();

			for (j = 0; j<wal_node; j++)
			{
				if (node[wall_node[j]].proc == i)
				{
					wal++;
				}
			}
			sprintf(filename2, "wall_%d.txt", i);
			outFile.open(filename2);
			for (j = 0; j<wal_node; j++)
			{
				if (node[wall_node[j]].proc == i)
				{
					outFile << wall_node[j] << "\n";
				}
			}
			outFile.close();

			for (j = 0; j<bou_node; j++)
			{
				if (node[boundary_node[j]].proc == i)
				{
					bou++;
				}
			}
			sprintf(filename2, "boundary_%d.txt", i);
			outFile.open(filename2);
			for (j = 0; j<bou_node; j++)
			{
				if (node[boundary_node[j]].proc == i)
				{
					outFile << boundary_node[j] << "\n";
				}
			}
			outFile.close();

			MPI_Send(&inl, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			MPI_Send(&out, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			MPI_Send(&wal, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			MPI_Send(&bou, 1, MPI_INT, i, i, MPI_COMM_WORLD);
		}


		/*******************************give boundary ID to nodes******************************************/
		for (i = 0; i<out_node; i++)
		{
			node[outlet_node[i]].ID = 2;
		}
		for (i = 0; i<inl_node; i++)
		{
			node[inlet_node[i]].ID = 1;
		}
		for (i = 0; i<wal_node; i++)
		{
			if (node[i].loc != 28)
			{
				node[wall_node[i]].ID = 3;
			}
		}
		for (i = 0; i<bou_node; i++)
		{
			node[boundary_node[i]].ID = 4;
		}
		

		MPI_Pack_size((NUMNP + 10) * 5, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		Bcast_buffer3 = new char[memsize];
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Pack(&jacobian[i].x_zeta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_eta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].x_xi, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_zeta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].y_eta, 1, MPI_DOUBLE, Bcast_buffer3, memsize, &position, MPI_COMM_WORLD);
		}

		MPI_Pack_size((NUMNP + 10) * 4, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		Bcast_buffer4 = new char[memsize];
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Pack(&jacobian[i].y_xi, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_zeta, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_eta, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&jacobian[i].z_xi, 1, MPI_DOUBLE, Bcast_buffer4, memsize, &position, MPI_COMM_WORLD);

		}
		/****************************************************************************************************/
		//	iterations = 100000;
		restart = 2;


	}
	MPI_Barrier(MPI_COMM_WORLD);
	/***********************************Broadcast NUMNP and NELEM**********************************/
	MPI_Bcast(&NUMNP, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&NELEM, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&restart, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Pack_size((NUMNP + 10) * 3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
	if (myrank >0)
	{
		Bcast_buffer1 = new char[memsize];
	}
	MPI_Bcast(Bcast_buffer1, (NUMNP + 10) * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		delete[] Bcast_buffer1;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Pack_size((NUMNP + 10) * 10, MPI_INT, MPI_COMM_WORLD, &memsize);
	if (myrank >0)
	{
		Bcast_buffer2 = new char[memsize];
	}
	MPI_Bcast(Bcast_buffer2, (NUMNP + 10) * 10, MPI_INT, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		delete[] Bcast_buffer2;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Pack_size((NUMNP + 10) * 5, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
	if (myrank >0)
	{
		Bcast_buffer3 = new char[memsize];
	}
	MPI_Bcast(Bcast_buffer3, (NUMNP + 10) * 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		delete[] Bcast_buffer3;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Pack_size((NUMNP + 10) * 4, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
	if (myrank >0)
	{
		Bcast_buffer4 = new char[memsize];
	}
	MPI_Bcast(Bcast_buffer4, (NUMNP + 10) * 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (myrank == 0)
	{
		delete[] Bcast_buffer4;
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0)
	{
		printf("done proper broadcasting\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	/********************************************************************************************/
	/********************************************************************************************/
	if (myrank>0)
	{
		tip.resize(5000);
		tip_recv.resize(5000);

		MPI_Recv(&inl_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&out_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&wal_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);
		MPI_Recv(&bou_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);

		d_node.resize((NUMNP + 10));
		for (i = 0; i <= NUMNP; i++)
		{
			d_node[i].x = 0;
			d_node[i].y = 0;
			d_node[i].z = 0;
			d_node[i].e[0] = 0;
			d_node[i].e[1] = 0;
			d_node[i].e[2] = 0;
			d_node[i].e[3] = 0;
			d_node[i].e[4] = 0;
			d_node[i].e[5] = 0;
			d_node[i].e[6] = 0;
			d_node[i].e[7] = 0;
			d_node[i].ID = 0;
			d_node[i].n_n[0] = 0;
			d_node[i].n_n[1] = 0;
			d_node[i].n_n[2] = 0;
			d_node[i].n_n[3] = 0;
			d_node[i].n_n[4] = 0;
			d_node[i].n_n[5] = 0;
			d_node[i].corner_ID = 0;
			d_node[i].val = 0;
			d_node[i].loc = 0;
			d_node[i].proc = 0;
			d_node[i].local = 0;
			d_node[i].global = 0;
			d_node[i].req = 0;
		}
		MPI_Pack_size((NUMNP + 10) * 3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].x, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].y, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer1, memsize, &position, &d_node[i].z, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		}
		delete[] Bcast_buffer1;

		MPI_Pack_size((NUMNP + 10) * 10, MPI_INT, MPI_COMM_WORLD, &memsize);
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[0], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[1], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[2], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[3], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[4], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].n_n[5], 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].loc, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].ID, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].proc, 1, MPI_INT, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer2, memsize, &position, &d_node[i].corner_ID, 1, MPI_INT, MPI_COMM_WORLD);
		}
		delete Bcast_buffer2;

		/******************************************************************************************************/
		/*******************READING PROC LIST AND COUNTING NODES IN SUBDOMAIN**********************************/
		sd_node = 0;
		sprintf(filename2, "nodes_proc_%d.txt", myrank);
		inFile.open(filename2);
		for (i = 1; i <= NUMNP; i++)
		{
			getline(inFile, line2, '\n');
			istringstream streamA(line2);
			streamA >> d_node[i].proc;
			if (d_node[i].proc == myrank)
			{
				sd_node++;
			}
		}
		inFile.close();
		sd_node = sd_node;
		tmp_node.resize(sd_node + 10);

		for (i = 0; i <= sd_node + 8; i++)
		{
			tmp_node[i].x = 0;
			tmp_node[i].y = 0;
			tmp_node[i].z = 0;
			tmp_node[i].e[0] = 0;
			tmp_node[i].e[1] = 0;
			tmp_node[i].e[2] = 0;
			tmp_node[i].e[3] = 0;
			tmp_node[i].e[4] = 0;
			tmp_node[i].e[5] = 0;
			tmp_node[i].e[6] = 0;
			tmp_node[i].e[7] = 0;
			tmp_node[i].ID = 0;
			tmp_node[i].n_n[0] = 0;
			tmp_node[i].n_n[1] = 0;
			tmp_node[i].n_n[2] = 0;
			tmp_node[i].n_n[3] = 0;
			tmp_node[i].n_n[4] = 0;
			tmp_node[i].n_n[5] = 0;
			//	tmp_node[i].location=0;
			//	tmp_node[i].corner=0;
			tmp_node[i].corner_ID = 0;
			tmp_node[i].val = 0;
			tmp_node[i].loc = 0;
			//	tmp_node[i].weno=0;
			tmp_node[i].proc = 0;
			tmp_node[i].local = 0;
			tmp_node[i].global = 0;
			tmp_node[i].req = 0;
		}

		j = 1;
		for (i = 1; i <= NUMNP; i++)
		{
			if (d_node[i].proc == myrank)
			{
				tmp_node[j].global = i;
				d_node[i].local = j;
				j++;
			}
		}
		/***************************************************************************************************************/
		h = sd_node;
		k = 0;
		/**********************************counting subdomain ghost nodes***********************************************/
		h = sd_node;
		k = 0;
		for (i = 1; i <= sd_node; i++)
		{
			if (d_node[tmp_node[i].global].corner_ID != 0 || d_node[tmp_node[i].global].loc != 0)
			{
				k++;
			}
			for (j = 0; j <= 5; j++)
			{
				if (d_node[tmp_node[i].global].n_n[j] != 0 && d_node[d_node[tmp_node[i].global].n_n[j]].proc != myrank && d_node[d_node[tmp_node[i].global].n_n[j]].proc != 0)
				{
					h = h + 6;
				}
			}
		}
		nodes = h;
		sd_gh = h - sd_node;

		/***********************************************************************************************/
		tmp_c.resize(size + 1);
		c.resize(size + 1);
		b.resize(size + 1);
		proc_node.resize(size + 1);
		gb.resize(size + 1);
		loc_dat.resize(size + 1);
		recv_b.resize(size + 1);
		recv_c.resize(size + 1);


		for (i = 0; i <= size; i++)
		{
			b[i].resize(sd_gh + 10);
			gb[i].resize(sd_gh + 10);
			loc_dat[i].resize(sd_gh + 10);
			recv_b[i].resize(sd_gh + 10);
		}
		node.resize(nodes + 10 + (k * 4));
		singular.resize(nodes + 10 + (k * 4));
		sd_gh_node.resize(nodes + 10 + (k * 4));

		all_bou_node = k;
		temp_boundary_nodes.resize(k + 10);

		for (i = 0; i <= nodes + 7 + (k * 4); i++)
		{
			node[i].x = 0;
			node[i].y = 0;
			node[i].z = 0;
			node[i].e[0] = 0;
			node[i].e[1] = 0;
			node[i].e[2] = 0;
			node[i].e[3] = 0;
			node[i].e[4] = 0;
			node[i].e[5] = 0;
			node[i].e[6] = 0;
			node[i].e[7] = 0;
			node[i].ID = 0;
			node[i].n_n[0] = 0;
			node[i].n_n[1] = 0;
			node[i].n_n[2] = 0;
			node[i].n_n[3] = 0;
			node[i].n_n[4] = 0;
			node[i].n_n[5] = 0;
			node[i].corner_ID = 0;
			node[i].val = 0;
			node[i].loc = 0;
			node[i].proc = 0;
			node[i].local = 0;
			node[i].global = 0;
			node[i].req = 0;
			singular[i].n_n[0] = 0;
			singular[i].n_n[1] = 0;
			singular[i].n_n[2] = 0;
			singular[i].n_n[3] = 0;
			singular[i].n_n[4] = 0;
			singular[i].n_n[5] = 0;
		}
		/**********************************************************************************************/
		j = 1;
		for (i = 1; i <= NUMNP; i++)
		{
			if (d_node[i].proc == myrank)
			{
				node[j].global = i;
				node[j].loc = d_node[i].loc;
				node[j].val = d_node[i].val;
				node[j].corner_ID = d_node[i].corner_ID;
				d_node[i].local = j;
				j++;
			}
		}
		tmp_node.clear();
		

		/**************************Recieving nodal data of this processors domain****************************************/
		h = sd_node + 1;
		gh = sd_node;
		k = 0;
		temp = 0;
		for (i = 1; i <= sd_node; i++)
		{
			node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
			node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
			node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
			node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
			node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
			node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
			node[i].val = d_node[node[i].global].val;
			node[i].loc = d_node[node[i].global].loc;
			node[i].corner_ID = d_node[node[i].global].corner_ID;
			node[i].ID = d_node[node[i].global].ID;
			node[i].proc = d_node[node[i].global].proc;
			if (d_node[node[i].global].corner_ID != 0 || d_node[node[i].global].loc != 0)
			{
				temp_boundary_nodes[temp] = i;
				temp++;
			}
		}

		/***************GLOBAL AND LOCAL NUMBERING OF NODES AND ITS NEIGHBOURS*************************/
		/**********DOMAIN NEIGHBOUR AND SUBDOMAIN BOUNDARY POINTS PROCS IDENTIFICATION*****************/

		subdomain_neighbour();


		/***********************************************************************************************/

		for (i = 1; i <= nodes; i++)
		{
			node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
			node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
			node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
			node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
			node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
			node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
			node[i].loc = d_node[node[i].global].loc;
			node[i].x = d_node[node[i].global].x;
			node[i].y = d_node[node[i].global].y;
			node[i].z = d_node[node[i].global].z;
			node[i].val = d_node[node[i].global].val;
			node[i].corner_ID = d_node[node[i].global].corner_ID;
			node[i].ID = d_node[node[i].global].ID;
		}

		/*****************SEGREGATING NODES OF OTHER PROCESSORS*******************************************/
		for (i = 0; i <= size; i++)
		{
			tmp_c[i] = 0;
		}
		for (j = sd_node + 1; j< h; j++)
		{
			tmp_c[sd_gh_node[j].proc]++;
		}
		//proc_node = NULL;
		j = 0;
		for (i = 1; i<size; i++)
		{
			if (tmp_c[i] != 0 && i != myrank)
			{
				c[j] = i;
				proc_node[c[j]] = tmp_c[i];
				d[j] = j;
				j++;
			}
		}
		neigh_pro = j;

		/***********************Arranging list of neighbour procs in acending order********************/
		for (i = 0; i<neigh_pro; i++)
		{
			for (j = 0; j<neigh_pro; j++)
			{
				if (c[j] > c[j + 1] && neigh_pro >(j + 1))
				{
					temp = c[j + 1];
					c[j + 1] = c[j];
					c[j] = temp;
					temp = d[j + 1];
					d[j + 1] = d[j];
					d[j] = temp;
				}
			}
		}
		/***********************************************************************************************/

		MPI_Send(&neigh_pro, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);

		position = 0;
		MPI_Pack_size(neigh_pro, MPI_INT, MPI_COMM_WORLD, &memsize);
		buffer = new char[memsize];
		for (i = 0; i<neigh_pro; i++)
		{
			MPI_Pack(&c[i], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
		}
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		delete[] buffer;

		all_neigh_pro.resize(size + 1);
		all_c_list.resize(size + 1);

		for (i = 1; i <= size - 1; i++)
		{
			MPI_Recv(&all_neigh_pro[i], 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, &status);

			all_c_list[i].resize(all_neigh_pro[i] + 1);
			position = 0;
			MPI_Pack_size(all_neigh_pro[i], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Recv(buffer, memsize, MPI_PACKED, 0, myrank, MPI_COMM_WORLD, &status);
			for (j = 1; j <= all_neigh_pro[i]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &glob, 1, MPI_INT, MPI_COMM_WORLD);
				all_c_list[i][j] = glob;
			}
			delete[] buffer;
		}

		for (i = 1; i <= size - 1; i++)
		{
			for (j = 1; j <= all_neigh_pro[i]; j++)
			{
				if (all_c_list[i][j] == myrank && proc_node[i] == 0 && i != myrank)
				{
					neigh_pro++;
					c[neigh_pro - 1] = i;
					proc_node[i] = 0;
				}
			}
		}


		for (i = 0; i<neigh_pro; i++)
		{
			for (j = 0; j<neigh_pro; j++)
			{
				if (c[j] > c[j + 1] && neigh_pro >(j + 1))
				{
					temp = c[j + 1];
					c[j + 1] = c[j];
					c[j] = temp;
					temp = d[j + 1];
					d[j + 1] = d[j];
					d[j] = temp;
				}
			}
		}

		for (j = 0; j< neigh_pro; j++)
		{
			k = 0;
			for (i = sd_node + 1; i <= gh; i++)
			{
				if (sd_gh_node[i].proc == c[j])
				{
					b[j][k] = sd_gh_node[i].global;
					k++;
				}
			}
		}

		/****************************data send to master for proc request checking*********************/
		temp = 0;
		for (i = 0; i<neigh_pro; i++)
		{
			temp = temp + proc_node[c[i]];
		}
		MPI_Send(&sd_gh, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);
		position = 0;
		MPI_Pack_size(temp, MPI_INT, MPI_COMM_WORLD, &memsize);
		buffer = new char[memsize];
		for (i = 0; i<neigh_pro; i++)
		{
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Pack(&b[i][j], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
			}
		}
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		delete[] buffer;

		/**********************SEND RECIEVE LOCAL NODE NEIGHBOUR DATA**********************************/
		position = 0;
		no_of_tip_recv = 0;
		k = 0;
		for (i = 0; i<neigh_pro; i++)
		{
			position = 0;
			MPI_Sendrecv(&proc_node[c[i]], 1, MPI_INT, c[i], c[i], &recv_c[c[i]], 1, MPI_INT, c[i], myrank, MPI_COMM_WORLD, &status);
			MPI_Pack_size(proc_node[c[i]], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Pack_size(recv_c[c[i]], MPI_INT, MPI_COMM_WORLD, &memsize1);
			buffer2 = new char[memsize1];                              /***********carefull with buffer1 and buffer2******************/
			position = 0;
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Pack(&b[i][j], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
				recv_b[i][j] = d_node[b[i][j]].local;
				if (d_node[b[i][j]].loc >= 27 && d_node[b[i][j]].loc <= 32)
				{
					tip_recv[++k].node = d_node[b[i][j]].local;
					no_of_tip_recv++;
				}
			}
			MPI_Sendrecv(buffer, position, MPI_PACKED, c[i], c[i], buffer2, memsize1, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Unpack(buffer2, memsize1, &position, &glob, 1, MPI_INT, MPI_COMM_WORLD);
				gb[i][j] = glob;    /*************Neigh processor's needs local node number of nodes with these global numbers frm current procs ***********/
			}
			delete[] buffer;
			delete[] buffer2;
		}
		position = 0;
		no_of_tip_send = 0;
		k = 0;
		for (i = 0; i<neigh_pro; i++)
		{
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				loc_dat[i][j] = d_node[gb[i][j]].local;    /************List of node numbers to be sent to neighbour processor********/
				if (d_node[gb[i][j]].loc >= 27 && d_node[gb[i][j]].loc <= 32)
				{
					no_of_tip_send++;
					tip[++k].i = i;
					tip[k].j = j;
					tip[k].node = d_node[gb[i][j]].local;
					tip[k].proc = d_node[gb[i][j]].proc;
				}
			}
		}

		for (i = 0; i<neigh_pro; i++)
		{
			MPI_Pack_size(recv_c[c[i]], MPI_INT, MPI_COMM_WORLD, &memsize1);
			buffer2 = new char[memsize1];                                          /***********carefull with buffer1 and buffer2******************/
			MPI_Pack_size(proc_node[c[i]], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			position = 0;
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&loc_dat[i][j], 1, MPI_INT, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}

			MPI_Sendrecv(buffer2, memsize1, MPI_PACKED, c[i], c[i], buffer, memsize, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &loca, 1, MPI_INT, MPI_COMM_WORLD);
				b[i][j] = loca;
			}
			delete[] buffer;
			delete[] buffer2;
		}


		/**************************Obtaining boundary nodes data from master to slave*************************/
		sprintf(filename2, "inlet_%d.txt", myrank);
		inFile.open(filename2);
		inlet_node.resize(inl_node + 1);
		i = 0;
		while (inFile.good())
		{
			getline(inFile, line, '\n');
			istringstream streamA(line);
			streamA >> temp7;
			inlet_node[i] = d_node[temp7].local;
			i++;
		}
		inFile.close();

		sprintf(filename2, "outlet_%d.txt", myrank);
		inFile.open(filename2);
		outlet_node.resize(out_node + 1);
		i = 0;
		while (inFile.good())
		{
			getline(inFile, line, '\n');
			istringstream streamA(line);
			streamA >> temp7;
			outlet_node[i] = d_node[temp7].local;
			i++;
		}
		inFile.close();

		sprintf(filename2, "boundary_%d.txt", myrank);
		inFile.open(filename2);
		boundary_node.resize(bou_node + 1);
		i = 0;
		while (inFile.good())
		{
			getline(inFile, line, '\n');
			istringstream streamA(line);
			streamA >> temp7;
			boundary_node[i] = d_node[temp7].local;
			i++;
		}
		inFile.close();

		sprintf(filename2, "wall_%d.txt", myrank);
		inFile.open(filename2);
		wall_node.resize(wal_node + 1);
		i = 0;
		while (inFile.good())
		{
			getline(inFile, line, '\n');
			istringstream streamA(line);
			streamA >> temp7;
			wall_node[i] = d_node[temp7].local;
			i++;
		}
		inFile.close();

		/************************sub_domain boundary conditions*****************************/
		sd_inl_node = 0;
		sd_out_node = 0;
		sd_wal_node = 0;
		sd_bou_node = 0;
		for (i = sd_node + 1; i<nodes; i++)
		{
			if (node[i].ID == 1)
			{
				sd_inl_node++;
			}
			if (node[i].ID == 2)
			{
				sd_out_node++;
			}
			if (node[i].ID == 3 || node[i].ID == 10)
			{
				sd_wal_node++;
			}
			if (node[i].ID == 4)
			{
				sd_bou_node++;
			}
		}

		sd_inlet_node.resize(sd_inl_node + 1);
		sd_outlet_node.resize(sd_out_node + 1);
		sd_wall_node.resize(sd_wal_node + 1);
		sd_boundary_node.resize(sd_bou_node + 1);

		temp = 0;
		for (i = sd_node + 1; i <= nodes; i++)
		{
			if (node[i].ID == 1)
			{
				sd_inlet_node[temp++] = i;
			}
		}
		temp = 0;
		for (i = sd_node + 1; i <= nodes; i++)
		{
			if (node[i].ID == 2)
			{
				sd_outlet_node[temp++] = i;
			}
		}
		temp = 0;
		for (i = sd_node + 1; i <= nodes; i++)
		{
			if (node[i].ID == 3 || node[i].ID == 10)
			{
				sd_wall_node[temp++] = i;
			}
		}
		temp = 0;
		for (i = sd_node + 1; i <= nodes; i++)
		{
			if (node[i].ID == 4)
			{
				sd_boundary_node[temp++] = i;
			}
		}
		/*************************************************************************************/
		temp = all_bou_node;
		all_bou_node = all_bou_node + sd_inl_node + sd_out_node + sd_wal_node + sd_bou_node;
		all_boundary_nodes.resize(all_bou_node + k + 10);
		for (i = 0; i<all_bou_node; i++)
		{
			all_boundary_nodes[i] = temp_boundary_nodes[i];
		}
		for (i = 0; i<sd_inl_node; i++)
		{
			all_boundary_nodes[temp++] = sd_inlet_node[i];
		}
		for (i = 0; i<sd_out_node; i++)
		{
			all_boundary_nodes[temp++] = sd_outlet_node[i];
		}
		for (i = 0; i<sd_wal_node; i++)
		{
			all_boundary_nodes[temp++] = sd_wall_node[i];
		}
		for (i = 0; i<sd_bou_node; i++)
		{
			all_boundary_nodes[temp++] = sd_boundary_node[i];
		}

		/*************************************************************************************/
		metric.resize(nodes + 10 + (all_bou_node * 5));
		jacobian.resize(nodes + 10 + (all_bou_node * 5));
		deter.resize(nodes + 10 + (all_bou_node * 5));
		det.resize(nodes + 10 + (all_bou_node * 5));
		Q.resize(nodes + 10 + (all_bou_node * 5));

		for (i = 0; i<(nodes + 10 + (all_bou_node * 5)); i++)
		{
			Q[i].resize(10);
		}

		t_metric.resize(NUMNP + 10);
		t_jacobian.resize(NUMNP + 10);
		t_det.resize(NUMNP + 10);

		MPI_Pack_size((NUMNP + 10) * 5, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].x_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].y_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer3, memsize, &position, &t_jacobian[i].y_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		}
		delete[] Bcast_buffer3;

		MPI_Pack_size((NUMNP + 10) * 4, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		position = 0;
		for (i = 1; i <= NUMNP; i++)
		{
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].y_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_zeta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_eta, 1, MPI_DOUBLE, MPI_COMM_WORLD);
			MPI_Unpack(Bcast_buffer4, memsize, &position, &t_jacobian[i].z_xi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
		}
		delete[] Bcast_buffer4;


		for (i = 1; i <= NUMNP; i++)
		{
			t_det[i] = (((t_jacobian[i].x_zeta)*((t_jacobian[i].y_eta*t_jacobian[i].z_xi) - (t_jacobian[i].y_xi*t_jacobian[i].z_eta))) + ((t_jacobian[i].x_eta)*((t_jacobian[i].y_xi*t_jacobian[i].z_zeta) - (t_jacobian[i].y_zeta*t_jacobian[i].z_xi))) + ((t_jacobian[i].x_xi)*((t_jacobian[i].y_zeta*t_jacobian[i].z_eta) - (t_jacobian[i].y_eta*t_jacobian[i].z_zeta))));
			t_metric[i].zeta_x = (1.0 / t_det[i])*(t_jacobian[i].y_eta*t_jacobian[i].z_xi - t_jacobian[i].y_xi*t_jacobian[i].z_eta);
			t_metric[i].zeta_y = (1.0 / t_det[i])*(t_jacobian[i].z_eta*t_jacobian[i].x_xi - t_jacobian[i].z_xi*t_jacobian[i].x_eta);
			t_metric[i].zeta_z = (1.0 / t_det[i])*(t_jacobian[i].x_eta*t_jacobian[i].y_xi - t_jacobian[i].x_xi*t_jacobian[i].y_eta);
			t_metric[i].eta_x = (1.0 / t_det[i])*(t_jacobian[i].y_xi*t_jacobian[i].z_zeta - t_jacobian[i].y_zeta*t_jacobian[i].z_xi);
			t_metric[i].eta_y = (1.0 / t_det[i])*(t_jacobian[i].z_xi*t_jacobian[i].x_zeta - t_jacobian[i].z_zeta*t_jacobian[i].x_xi);
			t_metric[i].eta_z = (1.0 / t_det[i])*(t_jacobian[i].x_xi*t_jacobian[i].y_zeta - t_jacobian[i].x_zeta*t_jacobian[i].y_xi);
			t_metric[i].xi_x = (1.0 / t_det[i])*(t_jacobian[i].y_zeta*t_jacobian[i].z_eta - t_jacobian[i].y_eta*t_jacobian[i].z_zeta);
			t_metric[i].xi_y = (1.0 / t_det[i])*(t_jacobian[i].z_zeta*t_jacobian[i].x_eta - t_jacobian[i].z_eta*t_jacobian[i].x_zeta);
			t_metric[i].xi_z = (1.0 / t_det[i])*(t_jacobian[i].x_zeta*t_jacobian[i].y_eta - t_jacobian[i].x_eta*t_jacobian[i].y_zeta);

		}

		/****************************READING TRANSFORMATION DATA***********************************/
		i = 1;
		for (i = 1; i <= nodes; i++)
		{
			jacobian[i].x_zeta = t_jacobian[node[i].global].x_zeta;
			jacobian[i].x_eta = t_jacobian[node[i].global].x_eta;
			jacobian[i].x_xi = t_jacobian[node[i].global].x_xi;
			jacobian[i].y_zeta = t_jacobian[node[i].global].y_zeta;
			jacobian[i].y_eta = t_jacobian[node[i].global].y_eta;
			jacobian[i].y_xi = t_jacobian[node[i].global].y_xi;
			jacobian[i].z_zeta = t_jacobian[node[i].global].z_zeta;
			jacobian[i].z_eta = t_jacobian[node[i].global].z_eta;
			jacobian[i].z_xi = t_jacobian[node[i].global].z_xi;
			metric[i].zeta_x = t_metric[node[i].global].zeta_x;
			metric[i].eta_x = t_metric[node[i].global].eta_x;
			metric[i].xi_x = t_metric[node[i].global].xi_x;
			metric[i].zeta_y = t_metric[node[i].global].zeta_y;
			metric[i].eta_y = t_metric[node[i].global].eta_y;
			metric[i].xi_y = t_metric[node[i].global].xi_y;
			metric[i].zeta_z = t_metric[node[i].global].zeta_z;
			metric[i].eta_z = t_metric[node[i].global].eta_z;
			metric[i].xi_z = t_metric[node[i].global].xi_z;
			det[i] = t_det[node[i].global];
		}
		/****************************Free arrays******************************************************/
		t_jacobian.clear();
		t_metric.clear();
		t_det.clear();
		temp_boundary_nodes.clear();
		sd_gh_node.clear();
		
		/**************************Memory Reallocation for ghost points*******************************/
		g_node = sd_node + sd_gh + 1;
		/*********************ghost elements and nodes creation***************************************/

		for (i = 0; i<all_bou_node; i++)
		{
			if (node[all_boundary_nodes[i]].corner_ID == 0)
			{
				if (node[all_boundary_nodes[i]].loc == 4)
				{
					j = 3; /********no node on WEST**************/
					k = 1;
					m = 4;
					gar5 = 1;
				}

				if (node[all_boundary_nodes[i]].loc == 3)
				{
					j = 2; /********no element on SOUTH**************/
					k = 0;
					m = 3;
					gar5 = 2;
				}

				if (node[all_boundary_nodes[i]].loc == 2)
				{
					j = 1; /********no element on EAST**************/
					k = 3;
					m = 2;
					gar5 = 1;
				}

				if (node[all_boundary_nodes[i]].loc == 1)
				{
					j = 0; /********no element on NORTH**************/
					k = 2;
					m = 1;
					gar5 = 2;
				}

				if (node[all_boundary_nodes[i]].loc == 6)
				{
					j = 4; /********no element on NORTH**************/
					k = 5;
					m = 1;
					gar5 = 2;
				}

				if (node[all_boundary_nodes[i]].loc == 5)
				{
					j = 5; /********no element on NORTH**************/
					k = 4;
					m = 1;
					gar5 = 2;
				}

				if (node[all_boundary_nodes[i]].loc > 14 && node[all_boundary_nodes[i]].loc < 27)
				{
					gar5 = 4;
				}

				if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
				{
					gar5 = 3;
				}

				if ((node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 56))
				{
					gar5 = 5;
				}

				switch (gar5)
				{
				case 1:
					node[all_boundary_nodes[i]].n_n[j] = g_node;
					node[g_node].n_n[k] = all_boundary_nodes[i];
					metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
					det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;

					node[g_node - 1].n_n[j] = g_node;
					node[g_node].n_n[k] = g_node - 1;
					metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
					det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;

					node[g_node - 1].n_n[j] = g_node;
					node[g_node].n_n[k] = g_node - 1;
					metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
					det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;
					break;

				case 2:
					node[all_boundary_nodes[i]].n_n[j] = g_node;
					node[g_node].n_n[k] = all_boundary_nodes[i];
					metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
					det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;

					node[g_node - 1].n_n[j] = g_node;
					node[g_node].n_n[k] = g_node - 1;
					metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
					det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;

					node[g_node - 1].n_n[j] = g_node;
					node[g_node].n_n[k] = g_node - 1;
					metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
					metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
					metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
					metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
					metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
					metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
					metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
					metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
					metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
					det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
					node[g_node].loc = 100;
					g_node++;
					g_elem++;
					break;

				case 3:
					if (node[all_boundary_nodes[i]].loc == 27 || node[all_boundary_nodes[i]].loc == 28)
					{
						temp6 = 3;
					}
					if (node[all_boundary_nodes[i]].loc == 29)
					{
						temp6 = 4;
					}
					if (node[all_boundary_nodes[i]].loc == 30)
					{
						temp6 = 4;
					}

					for (temp = 0; temp<temp6; temp++)
					{
						if (temp == 0 && (node[all_boundary_nodes[i]].loc == 28 || node[all_boundary_nodes[i]].loc == 29 || node[all_boundary_nodes[i]].loc == 30))
						{
							j = 0; /********no element on NORTH**************/
							k = 2;
							m = 1;
						}
						if (temp == 1 && (node[all_boundary_nodes[i]].loc == 28 || node[all_boundary_nodes[i]].loc == 29 || node[all_boundary_nodes[i]].loc == 30))
						{
							j = 2; /********no element on SOUTH**************/
							k = 0;
							m = 3;
						}
						if (temp == 0 && node[all_boundary_nodes[i]].loc == 27)
						{
							j = 4; /********no element on NORTH**************/
							k = 5;
							m = 1;
						}
						if (temp == 1 && node[all_boundary_nodes[i]].loc == 27)
						{
							j = 5; /********no element on SOUTH**************/
							k = 4;
							m = 3;
						}
						if (temp == 2)
						{
							j = 1; /********no element on EAST**************/
							k = 3;
							m = 2;
						}
						if (temp == 3 && node[all_boundary_nodes[i]].loc == 29)
						{
							j = 5;
							k = 4;
						}
						if (temp == 3 && node[all_boundary_nodes[i]].loc == 30)
						{
							j = 4;
							k = 5;
						}
						if (temp <= 1 || temp == 3)
						{
							if (temp == 0)
							{
								singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
								singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
							}
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						}
						if (temp == 2)
						{
							singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
							singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						}
					}
					break;

				case 4:
					temp1 = 0;
					temp2 = 0;
					temp3 = 0;
					temp4 = 0;
					temp5 = 0;
					temp6 = 0;
					for (temp = 0; temp<2; temp++)
					{
						temp7 = 0;
						if (node[all_boundary_nodes[i]].n_n[3] == 0 && temp1 == 0 && temp7 == 0)
						{
							j = 3; /********no node on WEST**************/
							k = 1;
							m = 4;
							temp1++;
							temp7++;
						}
						else if (node[all_boundary_nodes[i]].n_n[2] == 0 && temp2 == 0 && temp7 == 0)
						{
							j = 2; /********no element on SOUTH**************/
							k = 0;
							m = 3;
							temp2++;
							temp7++;
						}
						else if (node[all_boundary_nodes[i]].n_n[1] == 0 && temp3 == 0 && temp7 == 0)
						{
							j = 1; /********no element on EAST**************/
							k = 3;
							m = 2;
							temp3++;
							temp7++;
						}
						else if (node[all_boundary_nodes[i]].n_n[0] == 0 && temp4 == 0 && temp7 == 0)
						{
							j = 0; /********no element on NORTH**************/
							k = 2;
							m = 1;
							temp4++;
							temp7++;
						}
						else if (node[all_boundary_nodes[i]].n_n[4] == 0 && temp5 == 0 && temp7 == 0)
						{
							j = 4; /********no element on Z-plus**************/
							k = 5;
							m = 1;
							temp5++;
							temp7++;
						}
						else if (node[all_boundary_nodes[i]].n_n[5] == 0 && temp6 == 0 && temp7 == 0)
						{
							j = 5; /********no element on Z-minus**************/
							k = 4;
							m = 1;
							temp6++;
							temp7++;
						}

						node[all_boundary_nodes[i]].n_n[j] = g_node;
						node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
						det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
					}
					break;

				case 5:
					if ((node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 36) || (node[all_boundary_nodes[i]].loc >= 41 && node[all_boundary_nodes[i]].loc <= 44) || (node[all_boundary_nodes[i]].loc >= 49 && node[all_boundary_nodes[i]].loc <= 52))
					{
						temp6 = 2;
					}
					if ((node[all_boundary_nodes[i]].loc >= 37 && node[all_boundary_nodes[i]].loc <= 40) || (node[all_boundary_nodes[i]].loc >= 45 && node[all_boundary_nodes[i]].loc <= 48))
					{
						temp6 = 3;
					}
					if (node[all_boundary_nodes[i]].loc >= 53 && node[all_boundary_nodes[i]].loc <= 56)
					{
						temp6 = 1;
					}

					for (temp = 0; temp<temp6; temp++)
					{
						if (temp == 0 && ((node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 36) || (node[all_boundary_nodes[i]].loc >= 37 && node[all_boundary_nodes[i]].loc <= 40)))
						{
							j = 2;
							k = 0;
							m = 1;
						}
						if (temp == 0 && ((node[all_boundary_nodes[i]].loc >= 41 && node[all_boundary_nodes[i]].loc <= 44) || (node[all_boundary_nodes[i]].loc >= 45 && node[all_boundary_nodes[i]].loc <= 48)))
						{
							j = 0;
							k = 2;
							m = 1;
						}
						if (temp == 0 && (node[all_boundary_nodes[i]].loc == 49 || node[all_boundary_nodes[i]].loc == 52))
						{
							j = 5;
							k = 4;
							m = 2;
						}
						if (temp == 0 && (node[all_boundary_nodes[i]].loc == 50 || node[all_boundary_nodes[i]].loc == 51))
						{
							j = 4;
							k = 5;
							m = 2;
						}
						if (temp == 1 && (node[all_boundary_nodes[i]].loc == 33 || node[all_boundary_nodes[i]].loc == 43 || node[all_boundary_nodes[i]].loc == 51 || node[all_boundary_nodes[i]].loc == 52 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 47 || node[all_boundary_nodes[i]].loc == 48))
						{
							j = 3;
							k = 1;
							m = 2;
						}
						if (temp == 1 && (node[all_boundary_nodes[i]].loc == 35 || node[all_boundary_nodes[i]].loc == 41 || node[all_boundary_nodes[i]].loc == 49 || node[all_boundary_nodes[i]].loc == 50 || node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 45 || node[all_boundary_nodes[i]].loc == 46))
						{
							j = 1;
							k = 3;
							m = 2;
						}
						if (temp == 1 && (node[all_boundary_nodes[i]].loc == 34 || node[all_boundary_nodes[i]].loc == 44))
						{
							j = 5;
							k = 4;
							m = 2;
						}
						if (temp == 1 && (node[all_boundary_nodes[i]].loc == 36 || node[all_boundary_nodes[i]].loc == 42))
						{
							j = 4;
							k = 5;
							m = 2;
						}
						if (temp == 2 && (node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 48 || node[all_boundary_nodes[i]].loc == 45))
						{
							j = 5;
							k = 4;
							m = 2;
						}
						if (temp == 2 && (node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 46 || node[all_boundary_nodes[i]].loc == 47))
						{
							j = 4;
							k = 5;
							m = 2;
						}

						if (temp <= 1 || temp == 3)
						{
							//	if (temp == 0)
							{
								singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
								singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
							}
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						}
						if (temp == 2)
						{
							singular[all_boundary_nodes[i]].n_n[j] = node[all_boundary_nodes[i]].n_n[j];
							singular[all_boundary_nodes[i]].n_n[k] = node[all_boundary_nodes[i]].n_n[k];
							node[all_boundary_nodes[i]].n_n[j] = g_node;
							node[g_node].n_n[k] = all_boundary_nodes[i];
							metric[g_node].zeta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_z;
							det[g_node] = det[singular[all_boundary_nodes[i]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;

							node[g_node - 1].n_n[j] = g_node;
							node[g_node].n_n[k] = g_node - 1;
							metric[g_node].zeta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
							metric[g_node].eta_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
							metric[g_node].xi_x = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
							metric[g_node].zeta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
							metric[g_node].eta_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
							metric[g_node].xi_y = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
							metric[g_node].zeta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
							metric[g_node].eta_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
							metric[g_node].xi_z = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
							det[g_node] = det[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
							node[g_node].loc = 100;
							g_node++;
							g_elem++;
						}
					}
					break;
				}
			}

			else if (node[all_boundary_nodes[i]].corner_ID != 0)
			{
				for (gar = 0; gar<3; gar++)
				{
					if (node[all_boundary_nodes[i]].corner_ID == 1)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 2)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 3)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 4)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;
						}
					}

					if (node[all_boundary_nodes[i]].corner_ID == 5)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 6)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 7)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 2;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 8)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 2;
						}
					}

					switch (gar5)
					{
					case 1:
						node[all_boundary_nodes[i]].n_n[j] = g_node;
						node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
						det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
						break;

					case 2:
						node[all_boundary_nodes[i]].n_n[j] = g_node;
						node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[g_node].zeta_x = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[all_boundary_nodes[i]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[all_boundary_nodes[i]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[all_boundary_nodes[i]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[all_boundary_nodes[i]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[all_boundary_nodes[i]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[all_boundary_nodes[i]].n_n[k]].xi_z;
						det[g_node] = det[node[all_boundary_nodes[i]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;

						node[g_node - 1].n_n[j] = g_node;
						node[g_node].n_n[k] = g_node - 1;
						metric[g_node].zeta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_x;
						metric[g_node].eta_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_x;
						metric[g_node].xi_x = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_x;
						metric[g_node].zeta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_y;
						metric[g_node].eta_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_y;
						metric[g_node].xi_y = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_y;
						metric[g_node].zeta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_z;
						metric[g_node].eta_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_z;
						metric[g_node].xi_z = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_z;
						det[g_node] = det[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]];
						node[g_node].loc = 100;
						g_node++;
						g_elem++;
						break;
					}
				}
			}
		}

		nodes = g_node;
		/**************************************evaluating metrics at half points***********************************/
		for (i = 1; i <= sd_node; i++)
		{
			metric[i].zeta_xip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_x) - 25.0*(metric[node[i].n_n[3]].zeta_x) + 150.0*(metric[i].zeta_x) + 150.0*(metric[node[i].n_n[1]].zeta_x) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_x) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_x));
			metric[i].zeta_yip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_y) - 25.0*(metric[node[i].n_n[3]].zeta_y) + 150.0*(metric[i].zeta_y) + 150.0*(metric[node[i].n_n[1]].zeta_y) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_y) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_y));
			metric[i].zeta_zip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_z) - 25.0*(metric[node[i].n_n[3]].zeta_z) + 150.0*(metric[i].zeta_z) + 150.0*(metric[node[i].n_n[1]].zeta_z) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_z) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].zeta_z));
			metric[i].eta_xip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_x) - 25.0*(metric[node[i].n_n[3]].eta_x) + 150.0*(metric[i].eta_x) + 150.0*(metric[node[i].n_n[1]].eta_x) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_x) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_x));
			metric[i].eta_yip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_y) - 25.0*(metric[node[i].n_n[3]].eta_y) + 150.0*(metric[i].eta_y) + 150.0*(metric[node[i].n_n[1]].eta_y) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_y) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_y));
			metric[i].eta_zip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_z) - 25.0*(metric[node[i].n_n[3]].eta_z) + 150.0*(metric[i].eta_z) + 150.0*(metric[node[i].n_n[1]].eta_z) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_z) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_z));
			metric[i].xi_xip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_x) - 25.0*(metric[node[i].n_n[3]].xi_x) + 150.0*(metric[i].xi_x) + 150.0*(metric[node[i].n_n[1]].xi_x) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_x) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_x));
			metric[i].xi_yip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_y) - 25.0*(metric[node[i].n_n[3]].xi_y) + 150.0*(metric[i].xi_y) + 150.0*(metric[node[i].n_n[1]].xi_y) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_y) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_y));
			metric[i].xi_zip = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_z) - 25.0*(metric[node[i].n_n[3]].xi_z) + 150.0*(metric[i].xi_z) + 150.0*(metric[node[i].n_n[1]].xi_z) - 25.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_z) + 3.0*(metric[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_z));
			deter[i].ip = (1.0 / 256.0)*(3.0*(det[node[node[i].n_n[3]].n_n[3]]) - 25.0*(det[node[i].n_n[3]]) + 150.0*(det[i]) + 150.0*(det[node[i].n_n[1]]) - 25.0*(det[node[node[i].n_n[1]].n_n[1]]) + 3.0*(det[node[node[node[i].n_n[1]].n_n[1]].n_n[1]]));

			metric[i].zeta_xim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_x) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_x) + 150.0*(metric[node[i].n_n[3]].zeta_x) + 150.0*(metric[i].zeta_x) - 25.0*(metric[node[i].n_n[1]].zeta_x) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_x));
			metric[i].zeta_yim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_y) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_y) + 150.0*(metric[node[i].n_n[3]].zeta_y) + 150.0*(metric[i].zeta_y) - 25.0*(metric[node[i].n_n[1]].zeta_y) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_y));
			metric[i].zeta_zim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].zeta_z) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].zeta_z) + 150.0*(metric[node[i].n_n[3]].zeta_z) + 150.0*(metric[i].zeta_z) - 25.0*(metric[node[i].n_n[1]].zeta_z) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].zeta_z));
			metric[i].eta_xim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_x) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_x) + 150.0*(metric[node[i].n_n[3]].eta_x) + 150.0*(metric[i].eta_x) - 25.0*(metric[node[i].n_n[1]].eta_x) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_x));
			metric[i].eta_yim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_y) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_y) + 150.0*(metric[node[i].n_n[3]].eta_y) + 150.0*(metric[i].eta_y) - 25.0*(metric[node[i].n_n[1]].eta_y) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_y));
			metric[i].eta_zim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_z) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].eta_z) + 150.0*(metric[node[i].n_n[3]].eta_z) + 150.0*(metric[i].eta_z) - 25.0*(metric[node[i].n_n[1]].eta_z) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].eta_z));
			metric[i].xi_xim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_x) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_x) + 150.0*(metric[node[i].n_n[3]].xi_x) + 150.0*(metric[i].xi_x) - 25.0*(metric[node[i].n_n[1]].xi_x) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_x));
			metric[i].xi_yim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_y) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_y) + 150.0*(metric[node[i].n_n[3]].xi_y) + 150.0*(metric[i].xi_y) - 25.0*(metric[node[i].n_n[1]].xi_y) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_y));
			metric[i].xi_zim = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_z) - 25.0*(metric[node[node[i].n_n[3]].n_n[3]].xi_z) + 150.0*(metric[node[i].n_n[3]].xi_z) + 150.0*(metric[i].xi_z) - 25.0*(metric[node[i].n_n[1]].xi_z) + 3.0*(metric[node[node[i].n_n[1]].n_n[1]].xi_z));
			deter[i].im = (1.0 / 256.0)*(3.0*(det[node[node[node[i].n_n[3]].n_n[3]].n_n[3]]) - 25.0*(det[node[node[i].n_n[3]].n_n[3]]) + 150.0*(det[node[i].n_n[3]]) + 150.0*(det[i]) - 25.0*(det[node[i].n_n[1]]) + 3.0*(det[node[node[i].n_n[1]].n_n[1]]));

			metric[i].zeta_xjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_x) - 25.0*(metric[node[i].n_n[2]].zeta_x) + 150.0*(metric[i].zeta_x) + 150.0*(metric[node[i].n_n[0]].zeta_x) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_x) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_x));
			metric[i].zeta_yjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_y) - 25.0*(metric[node[i].n_n[2]].zeta_y) + 150.0*(metric[i].zeta_y) + 150.0*(metric[node[i].n_n[0]].zeta_y) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_y) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_y));
			metric[i].zeta_zjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_z) - 25.0*(metric[node[i].n_n[2]].zeta_z) + 150.0*(metric[i].zeta_z) + 150.0*(metric[node[i].n_n[0]].zeta_z) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_z) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_z));
			metric[i].eta_xjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_x) - 25.0*(metric[node[i].n_n[2]].eta_x) + 150.0*(metric[i].eta_x) + 150.0*(metric[node[i].n_n[0]].eta_x) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_x) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_x));
			metric[i].eta_yjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_y) - 25.0*(metric[node[i].n_n[2]].eta_y) + 150.0*(metric[i].eta_y) + 150.0*(metric[node[i].n_n[0]].eta_y) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_y) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_y));
			metric[i].eta_zjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_z) - 25.0*(metric[node[i].n_n[2]].eta_z) + 150.0*(metric[i].eta_z) + 150.0*(metric[node[i].n_n[0]].eta_z) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_z) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].eta_z));
			metric[i].xi_xjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_x) - 25.0*(metric[node[i].n_n[2]].xi_x) + 150.0*(metric[i].xi_x) + 150.0*(metric[node[i].n_n[0]].xi_x) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_x) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_x));
			metric[i].xi_yjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_y) - 25.0*(metric[node[i].n_n[2]].xi_y) + 150.0*(metric[i].xi_y) + 150.0*(metric[node[i].n_n[0]].xi_y) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_y) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_y));
			metric[i].xi_zjp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_z) - 25.0*(metric[node[i].n_n[2]].xi_z) + 150.0*(metric[i].xi_z) + 150.0*(metric[node[i].n_n[0]].xi_z) - 25.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_z) + 3.0*(metric[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_z));
			deter[i].jp = (1.0 / 256.0)*(3.0*(det[node[node[i].n_n[2]].n_n[2]]) - 25.0*(det[node[i].n_n[2]]) + 150.0*(det[i]) + 150.0*(det[node[i].n_n[0]]) - 25.0*(det[node[node[i].n_n[0]].n_n[0]]) + 3.0*(det[node[node[node[i].n_n[0]].n_n[0]].n_n[0]]));

			metric[i].zeta_xjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_x) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_x) + 150.0*(metric[node[i].n_n[2]].zeta_x) + 150.0*(metric[i].zeta_x) - 25.0*(metric[node[i].n_n[0]].zeta_x) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_x));
			metric[i].zeta_yjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_y) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_y) + 150.0*(metric[node[i].n_n[2]].zeta_y) + 150.0*(metric[i].zeta_y) - 25.0*(metric[node[i].n_n[0]].zeta_y) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_y));
			metric[i].zeta_zjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_z) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].zeta_z) + 150.0*(metric[node[i].n_n[2]].zeta_z) + 150.0*(metric[i].zeta_z) - 25.0*(metric[node[i].n_n[0]].zeta_z) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].zeta_z));
			metric[i].eta_xjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_x) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_x) + 150.0*(metric[node[i].n_n[2]].eta_x) + 150.0*(metric[i].eta_x) - 25.0*(metric[node[i].n_n[0]].eta_x) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_x));
			metric[i].eta_yjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_y) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_y) + 150.0*(metric[node[i].n_n[2]].eta_y) + 150.0*(metric[i].eta_y) - 25.0*(metric[node[i].n_n[0]].eta_y) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_y));
			metric[i].eta_zjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].eta_z) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].eta_z) + 150.0*(metric[node[i].n_n[2]].eta_z) + 150.0*(metric[i].eta_z) - 25.0*(metric[node[i].n_n[0]].eta_z) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].eta_z));
			metric[i].xi_xjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_x) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_x) + 150.0*(metric[node[i].n_n[2]].xi_x) + 150.0*(metric[i].xi_x) - 25.0*(metric[node[i].n_n[0]].xi_x) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_x));
			metric[i].xi_yjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_y) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_y) + 150.0*(metric[node[i].n_n[2]].xi_y) + 150.0*(metric[i].xi_y) - 25.0*(metric[node[i].n_n[0]].xi_y) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_y));
			metric[i].xi_zjm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_z) - 25.0*(metric[node[node[i].n_n[2]].n_n[2]].xi_z) + 150.0*(metric[node[i].n_n[2]].xi_z) + 150.0*(metric[i].xi_z) - 25.0*(metric[node[i].n_n[0]].xi_z) + 3.0*(metric[node[node[i].n_n[0]].n_n[0]].xi_z));
			deter[i].jm = (1.0 / 256.0)*(3.0*(det[node[node[node[i].n_n[2]].n_n[2]].n_n[2]]) - 25.0*(det[node[node[i].n_n[2]].n_n[2]]) + 150.0*(det[node[i].n_n[2]]) + 150.0*(det[i]) - 25.0*(det[node[i].n_n[0]]) + 3.0*(det[node[node[i].n_n[0]].n_n[0]]));

			metric[i].zeta_xkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_x) - 25.0*(metric[node[i].n_n[5]].zeta_x) + 150.0*(metric[i].zeta_x) + 150.0*(metric[node[i].n_n[4]].zeta_x) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_x) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_x));
			metric[i].zeta_ykp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_y) - 25.0*(metric[node[i].n_n[5]].zeta_y) + 150.0*(metric[i].zeta_y) + 150.0*(metric[node[i].n_n[4]].zeta_y) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_y) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_y));
			metric[i].zeta_zkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_z) - 25.0*(metric[node[i].n_n[5]].zeta_z) + 150.0*(metric[i].zeta_z) + 150.0*(metric[node[i].n_n[4]].zeta_z) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_z) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_z));
			metric[i].eta_xkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_x) - 25.0*(metric[node[i].n_n[5]].eta_x) + 150.0*(metric[i].eta_x) + 150.0*(metric[node[i].n_n[4]].eta_x) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_x) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_x));
			metric[i].eta_ykp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_y) - 25.0*(metric[node[i].n_n[5]].eta_y) + 150.0*(metric[i].eta_y) + 150.0*(metric[node[i].n_n[4]].eta_y) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_y) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_y));
			metric[i].eta_zkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_z) - 25.0*(metric[node[i].n_n[5]].eta_z) + 150.0*(metric[i].eta_z) + 150.0*(metric[node[i].n_n[4]].eta_z) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_z) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_z));
			metric[i].xi_xkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_x) - 25.0*(metric[node[i].n_n[5]].xi_x) + 150.0*(metric[i].xi_x) + 150.0*(metric[node[i].n_n[4]].xi_x) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_x) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_x));
			metric[i].xi_ykp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_y) - 25.0*(metric[node[i].n_n[5]].xi_y) + 150.0*(metric[i].xi_y) + 150.0*(metric[node[i].n_n[4]].xi_y) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_y) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_y));
			metric[i].xi_zkp = (1.0 / 256.0)*(3.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_z) - 25.0*(metric[node[i].n_n[5]].xi_z) + 150.0*(metric[i].xi_z) + 150.0*(metric[node[i].n_n[4]].xi_z) - 25.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_z) + 3.0*(metric[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].xi_z));
			deter[i].kp = (1.0 / 256.0)*(3.0*(det[node[node[i].n_n[5]].n_n[5]]) - 25.0*(det[node[i].n_n[5]]) + 150.0*(det[i]) + 150.0*(det[node[i].n_n[4]]) - 25.0*(det[node[node[i].n_n[4]].n_n[4]]) + 3.0*(det[node[node[node[i].n_n[4]].n_n[4]].n_n[4]]));

			metric[i].zeta_xkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_x) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_x) + 150.0*(metric[node[i].n_n[5]].zeta_x) + 150.0*(metric[i].zeta_x) - 25.0*(metric[node[i].n_n[4]].zeta_x) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_x));
			metric[i].zeta_ykm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_y) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_y) + 150.0*(metric[node[i].n_n[5]].zeta_y) + 150.0*(metric[i].zeta_y) - 25.0*(metric[node[i].n_n[4]].zeta_y) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_y));
			metric[i].zeta_zkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_z) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].zeta_z) + 150.0*(metric[node[i].n_n[5]].zeta_z) + 150.0*(metric[i].zeta_z) - 25.0*(metric[node[i].n_n[4]].zeta_z) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].zeta_z));
			metric[i].eta_xkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_x) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_x) + 150.0*(metric[node[i].n_n[5]].eta_x) + 150.0*(metric[i].eta_x) - 25.0*(metric[node[i].n_n[4]].eta_x) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_x));
			metric[i].eta_ykm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_y) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_y) + 150.0*(metric[node[i].n_n[5]].eta_y) + 150.0*(metric[i].eta_y) - 25.0*(metric[node[i].n_n[4]].eta_y) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_y));
			metric[i].eta_zkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_z) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].eta_z) + 150.0*(metric[node[i].n_n[5]].eta_z) + 150.0*(metric[i].eta_z) - 25.0*(metric[node[i].n_n[4]].eta_z) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].eta_z));
			metric[i].xi_xkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_x) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_x) + 150.0*(metric[node[i].n_n[5]].xi_x) + 150.0*(metric[i].xi_x) - 25.0*(metric[node[i].n_n[4]].xi_x) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_x));
			metric[i].xi_ykm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_y) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_y) + 150.0*(metric[node[i].n_n[5]].xi_y) + 150.0*(metric[i].xi_y) - 25.0*(metric[node[i].n_n[4]].xi_y) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_y));
			metric[i].xi_zkm = (1.0 / 256.0)*(3.0*(metric[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].xi_z) - 25.0*(metric[node[node[i].n_n[5]].n_n[5]].xi_z) + 150.0*(metric[node[i].n_n[5]].xi_z) + 150.0*(metric[i].xi_z) - 25.0*(metric[node[i].n_n[4]].xi_z) + 3.0*(metric[node[node[i].n_n[4]].n_n[4]].xi_z));
			deter[i].km = (1.0 / 256.0)*(3.0*(det[node[node[node[i].n_n[5]].n_n[5]].n_n[5]]) - 25.0*(det[node[node[i].n_n[5]].n_n[5]]) + 150.0*(det[node[i].n_n[5]]) + 150.0*(det[i]) - 25.0*(det[node[i].n_n[4]]) + 3.0*(det[node[node[i].n_n[4]].n_n[4]]));

		}

		/********* Repeating for evaluating metrics at half points for ghost points***********************************/
		//	g_node=sd_node+sd_gh+1;
		for (i = 0; i<all_bou_node; i++)
		{
			gar = 1;
			if (node[all_boundary_nodes[i]].loc > 14 && node[all_boundary_nodes[i]].loc < 27)
			{
				gar++;
			}

			if (node[all_boundary_nodes[i]].corner_ID == 0)
			{
				temp1 = 0;
				temp2 = 0;
				temp3 = 0;
				temp4 = 0;
				temp5 = 0;
				temp6 = 0;
				for (temp7 = 0; temp7 < gar; temp7++)
				{
					if (node[node[all_boundary_nodes[i]].n_n[3]].loc == 100 && temp1 == 0)
					{
						j = 3; /********no node on WEST**************/
						k = 1;
						m = 4;
						gar5 = 1;
						gar4 = 4;
						temp1++;
					}

					if (node[node[all_boundary_nodes[i]].n_n[2]].loc == 100 && temp2 == 0)
					{
						j = 2; /********no element on SOUTH**************/
						k = 0;
						m = 3;
						gar5 = 2;
						gar4 = 3;
						temp2++;
					}

					if (node[node[all_boundary_nodes[i]].n_n[1]].loc == 100 && temp3 == 0)
					{
						j = 1; /********no element on EAST**************/
						k = 3;
						m = 2;
						gar5 = 1;
						gar4 = 2;
						temp3++;
					}

					if (node[node[all_boundary_nodes[i]].n_n[0]].loc == 100 && temp4 == 0)
					{
						j = 0; /********no element on NORTH**************/
						k = 2;
						m = 1;
						gar5 = 2;
						gar4 = 1;
						temp4++;
					}

					if (node[node[all_boundary_nodes[i]].n_n[4]].loc == 100 && temp5 == 0)
					{
						j = 4; /********no element on NORTH**************/
						k = 5;
						m = 1;
						gar5 = 3;
						gar4 = 5;
						temp5++;
					}

					if (node[node[all_boundary_nodes[i]].n_n[5]].loc == 100 && temp6 == 0)
					{
						j = 5; /********no element on NORTH**************/
						k = 4;
						m = 1;
						gar5 = 3;
						gar4 = 6;
						temp6++;
					}

					if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
					{
						gar5 = 4;
						gar4 = 2;
					}
					if (node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 56)
					{
						gar5 = 5;
					}

					switch (gar5)
					{
					case 1:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];										
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;						
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;						
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
						break;

					case 2:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
						break;

					case 3:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].km;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;
						break;

					case 4:
						for (temp = 0; temp<3; temp++)
						{
							if (temp == 0 && node[i].loc == 28)
							{
								j = 0; /********no element on NORTH**************/
								k = 2;
								m = 1;
							}
							if (temp == 1 && node[i].loc == 28)
							{
								j = 2; /********no element on SOUTH**************/
								k = 0;
								m = 3;
							}
							if (temp == 0 && node[i].loc == 27)
							{
								j = 4; /********no element on NORTH**************/
								k = 5;
								m = 1;
							}
							if (temp == 1 && node[i].loc == 27)
							{
								j = 5; /********no element on SOUTH**************/
								k = 4;
								m = 3;
							}
							if (temp == 2)
							{
								j = 1; /********no element on EAST**************/
								k = 3;
								m = 2;
							}
							if (temp <= 1)
							{
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
								deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
								deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].im;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
								deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
								deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
								deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
								deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;




								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;



								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

							}
							if (temp == 2)
							{
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
								deter[singular[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].im;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
								deter[singular[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;


								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;


								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
							}
						}
						break;

					case 5:
						if ((node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 36) || (node[all_boundary_nodes[i]].loc >= 41 && node[all_boundary_nodes[i]].loc <= 44) || (node[all_boundary_nodes[i]].loc >= 49 && node[all_boundary_nodes[i]].loc <= 52))
						{
							temp6 = 2;
						}
						if ((node[all_boundary_nodes[i]].loc >= 37 && node[all_boundary_nodes[i]].loc <= 40) || (node[all_boundary_nodes[i]].loc >= 45 && node[all_boundary_nodes[i]].loc <= 48))
						{
							temp6 = 3;
						}
						if (node[all_boundary_nodes[i]].loc >= 53 && node[all_boundary_nodes[i]].loc <= 56)
						{
							temp6 = 1;
						}

						for (temp = 0; temp<temp6; temp++)
						{
							if (temp == 0 && ((node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 36) || (node[all_boundary_nodes[i]].loc >= 37 && node[all_boundary_nodes[i]].loc <= 40)))
							{
								j = 2;
								k = 0;
								m = 1;
							}
							if (temp == 0 && ((node[all_boundary_nodes[i]].loc >= 41 && node[all_boundary_nodes[i]].loc <= 44) || (node[all_boundary_nodes[i]].loc >= 45 && node[all_boundary_nodes[i]].loc <= 48)))
							{
								j = 0;
								k = 2;
								m = 1;
							}

							if (temp == 1 && (node[all_boundary_nodes[i]].loc == 33 || node[all_boundary_nodes[i]].loc == 43 || node[all_boundary_nodes[i]].loc == 51 || node[all_boundary_nodes[i]].loc == 52 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 47 || node[all_boundary_nodes[i]].loc == 48))
							{
								j = 3;
								k = 1;
								m = 2;
							}
							if (temp == 1 && (node[all_boundary_nodes[i]].loc == 35 || node[all_boundary_nodes[i]].loc == 41 || node[all_boundary_nodes[i]].loc == 49 || node[all_boundary_nodes[i]].loc == 50 || node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 45 || node[all_boundary_nodes[i]].loc == 46))
							{
								j = 1;
								k = 3;
								m = 2;
							}
							if (temp == 0 && (node[all_boundary_nodes[i]].loc == 49 || node[all_boundary_nodes[i]].loc == 52))
							{
								j = 5;
								k = 4;
								m = 2;
							}
							if (temp == 0 && (node[all_boundary_nodes[i]].loc == 50 || node[all_boundary_nodes[i]].loc == 51))
							{
								j = 4;
								k = 5;
								m = 2;
							}
							if (temp == 1 && (node[all_boundary_nodes[i]].loc == 34 || node[all_boundary_nodes[i]].loc == 44))
							{
								j = 5;
								k = 4;
								m = 2;
							}
							if (temp == 1 && (node[all_boundary_nodes[i]].loc == 36 || node[all_boundary_nodes[i]].loc == 42))
							{
								j = 4;
								k = 5;
								m = 2;
							}
							if (temp == 2 && (node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 48 || node[all_boundary_nodes[i]].loc == 45))
							{
								j = 5;
								k = 4;
								m = 2;
							}
							if (temp == 2 && (node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 46 || node[all_boundary_nodes[i]].loc == 47))
							{
								j = 4;
								k = 5;
								m = 2;
							}

							if (j == 0 || j == 2)
							{
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
								deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
								deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].im;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
								deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
								deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
								deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;

								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
								metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
								metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
								deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;




								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;



								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

							}

							if (j == 1 || j == 3)
							{
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
								deter[singular[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].im;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
								deter[singular[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].km;


								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;


								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;
							}

							if (j == 4 || j == 5)
							{
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zip;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zip;
								deter[singular[all_boundary_nodes[i]].n_n[j]].ip = deter[singular[all_boundary_nodes[i]].n_n[k]].ip;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zim;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zim;
								deter[singular[all_boundary_nodes[i]].n_n[j]].im = deter[singular[all_boundary_nodes[i]].n_n[k]].im;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jp = deter[singular[all_boundary_nodes[i]].n_n[k]].jp;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_yjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zjm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zjm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].jm = deter[singular[all_boundary_nodes[i]].n_n[k]].jm;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkm;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkm;
								deter[singular[all_boundary_nodes[i]].n_n[j]].kp = deter[singular[all_boundary_nodes[i]].n_n[k]].km;

								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_xkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_ykp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].eta_zkp;
								metric[singular[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[singular[all_boundary_nodes[i]].n_n[k]].xi_zkp;
								deter[singular[all_boundary_nodes[i]].n_n[j]].km = deter[singular[all_boundary_nodes[i]].n_n[k]].kp;


								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;

								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;


								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
								metric[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
								deter[node[node[singular[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[singular[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

							}

						}
						break;
					}

					switch (gar4)
					{
					case 1:
						metric[all_boundary_nodes[i]].zeta_xjp = metric[all_boundary_nodes[i]].zeta_xjm;
						metric[all_boundary_nodes[i]].eta_xjp = metric[all_boundary_nodes[i]].eta_xjm;
						metric[all_boundary_nodes[i]].xi_xjp = metric[all_boundary_nodes[i]].xi_xjm;
						metric[all_boundary_nodes[i]].zeta_yjp = metric[all_boundary_nodes[i]].zeta_yjm;
						metric[all_boundary_nodes[i]].eta_yjp = metric[all_boundary_nodes[i]].eta_yjm;
						metric[all_boundary_nodes[i]].xi_yjp = metric[all_boundary_nodes[i]].xi_yjm;
						metric[all_boundary_nodes[i]].zeta_zjp = metric[all_boundary_nodes[i]].zeta_zjm;
						metric[all_boundary_nodes[i]].eta_zjp = metric[all_boundary_nodes[i]].eta_zjm;
						metric[all_boundary_nodes[i]].xi_zjp = metric[all_boundary_nodes[i]].xi_zjm;
						deter[all_boundary_nodes[i]].jp = deter[all_boundary_nodes[i]].jm;
						break;
					case 2:
						metric[all_boundary_nodes[i]].zeta_xip = metric[all_boundary_nodes[i]].zeta_xim;
						metric[all_boundary_nodes[i]].eta_xip = metric[all_boundary_nodes[i]].eta_xim;
						metric[all_boundary_nodes[i]].xi_xip = metric[all_boundary_nodes[i]].xi_xim;
						metric[all_boundary_nodes[i]].zeta_yip = metric[all_boundary_nodes[i]].zeta_yim;
						metric[all_boundary_nodes[i]].eta_yip = metric[all_boundary_nodes[i]].eta_yim;
						metric[all_boundary_nodes[i]].xi_yip = metric[all_boundary_nodes[i]].xi_yim;
						metric[all_boundary_nodes[i]].zeta_zip = metric[all_boundary_nodes[i]].zeta_zim;
						metric[all_boundary_nodes[i]].eta_zip = metric[all_boundary_nodes[i]].eta_zim;
						metric[all_boundary_nodes[i]].xi_zip = metric[all_boundary_nodes[i]].xi_zim;
						deter[all_boundary_nodes[i]].ip = deter[all_boundary_nodes[i]].im;
						break;
					case 3:
						metric[all_boundary_nodes[i]].zeta_xjm = metric[all_boundary_nodes[i]].zeta_xjp;
						metric[all_boundary_nodes[i]].eta_xjm = metric[all_boundary_nodes[i]].eta_xjp;
						metric[all_boundary_nodes[i]].xi_xjm = metric[all_boundary_nodes[i]].xi_xjp;
						metric[all_boundary_nodes[i]].zeta_yjm = metric[all_boundary_nodes[i]].zeta_yjp;
						metric[all_boundary_nodes[i]].eta_yjm = metric[all_boundary_nodes[i]].eta_yjp;
						metric[all_boundary_nodes[i]].xi_yjm = metric[all_boundary_nodes[i]].xi_yjp;
						metric[all_boundary_nodes[i]].zeta_zjm = metric[all_boundary_nodes[i]].zeta_zjp;
						metric[all_boundary_nodes[i]].eta_zjm = metric[all_boundary_nodes[i]].eta_zjp;
						metric[all_boundary_nodes[i]].xi_zjm = metric[all_boundary_nodes[i]].xi_zjp;
						deter[all_boundary_nodes[i]].jm = deter[all_boundary_nodes[i]].jp;
						break;
					case 4:
						metric[all_boundary_nodes[i]].zeta_xim = metric[all_boundary_nodes[i]].zeta_xip;
						metric[all_boundary_nodes[i]].eta_xim = metric[all_boundary_nodes[i]].eta_xip;
						metric[all_boundary_nodes[i]].xi_xim = metric[all_boundary_nodes[i]].xi_xip;
						metric[all_boundary_nodes[i]].zeta_yim = metric[all_boundary_nodes[i]].zeta_yip;
						metric[all_boundary_nodes[i]].eta_yim = metric[all_boundary_nodes[i]].eta_yip;
						metric[all_boundary_nodes[i]].xi_yim = metric[all_boundary_nodes[i]].xi_yip;
						metric[all_boundary_nodes[i]].zeta_zim = metric[all_boundary_nodes[i]].zeta_zip;
						metric[all_boundary_nodes[i]].eta_zim = metric[all_boundary_nodes[i]].eta_zip;
						metric[all_boundary_nodes[i]].xi_zim = metric[all_boundary_nodes[i]].xi_zip;
						deter[all_boundary_nodes[i]].im = deter[all_boundary_nodes[i]].ip;
						break;
					case 5:
						metric[all_boundary_nodes[i]].zeta_xkp = metric[all_boundary_nodes[i]].zeta_xkm;
						metric[all_boundary_nodes[i]].eta_xkp = metric[all_boundary_nodes[i]].eta_xkm;
						metric[all_boundary_nodes[i]].xi_xkp = metric[all_boundary_nodes[i]].xi_xkm;
						metric[all_boundary_nodes[i]].zeta_ykp = metric[all_boundary_nodes[i]].zeta_ykm;
						metric[all_boundary_nodes[i]].eta_ykp = metric[all_boundary_nodes[i]].eta_ykm;
						metric[all_boundary_nodes[i]].xi_ykp = metric[all_boundary_nodes[i]].xi_ykm;
						metric[all_boundary_nodes[i]].zeta_zkp = metric[all_boundary_nodes[i]].zeta_zkm;
						metric[all_boundary_nodes[i]].eta_zkp = metric[all_boundary_nodes[i]].eta_zkm;
						metric[all_boundary_nodes[i]].xi_zkp = metric[all_boundary_nodes[i]].xi_zkm;
						deter[all_boundary_nodes[i]].kp = deter[all_boundary_nodes[i]].km;
						break;
					case 6:
						metric[all_boundary_nodes[i]].zeta_xkm = metric[all_boundary_nodes[i]].zeta_xkp;
						metric[all_boundary_nodes[i]].eta_xkm = metric[all_boundary_nodes[i]].eta_xkp;
						metric[all_boundary_nodes[i]].xi_xkm = metric[all_boundary_nodes[i]].xi_xkp;
						metric[all_boundary_nodes[i]].zeta_ykm = metric[all_boundary_nodes[i]].zeta_ykp;
						metric[all_boundary_nodes[i]].eta_ykm = metric[all_boundary_nodes[i]].eta_ykp;
						metric[all_boundary_nodes[i]].xi_ykm = metric[all_boundary_nodes[i]].xi_ykp;
						metric[all_boundary_nodes[i]].zeta_zkm = metric[all_boundary_nodes[i]].zeta_zkp;
						metric[all_boundary_nodes[i]].eta_zkm = metric[all_boundary_nodes[i]].eta_zkp;
						metric[all_boundary_nodes[i]].xi_zkm = metric[all_boundary_nodes[i]].xi_zkp;
						deter[all_boundary_nodes[i]].km = deter[all_boundary_nodes[i]].kp;
						break;
					}


				}
			}

			else if (node[all_boundary_nodes[i]].corner_ID != 0)
			{
				for (gar = 0; gar<3; gar++)
				{
					if (node[all_boundary_nodes[i]].corner_ID == 1)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 2)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 3)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 4)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 2;
							j = 0;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;
						}
					}

					if (node[all_boundary_nodes[i]].corner_ID == 5)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 6)
					{
						if (gar == 0)
						{
							k = 1;
							j = 3;
							m = 8;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 7)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 10;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 4;
							j = 5;
							m = 10;
							gar5 = 3;
						}
					}
					else if (node[all_boundary_nodes[i]].corner_ID == 8)
					{
						if (gar == 0)
						{
							k = 3;
							j = 1;
							m = 9;
							gar5 = 1;
						}
						else if (gar == 1)
						{
							k = 0;
							j = 2;
							m = 11;
							gar5 = 2;
						}
						else if (gar == 2)
						{
							k = 5;
							j = 4;
							m = 11;
							gar5 = 3;
						}
					}
					switch (gar5)
					{
					case 1:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];										
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;						
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;						
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

						break;

					case 2:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

						break;

					case 3:
						//node[all_boundary_nodes[i]].n_n[j] = g_node;
						//node[g_node].n_n[k] = all_boundary_nodes[i];
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yip;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zip = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zip;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zip = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zip;
						deter[node[all_boundary_nodes[i]].n_n[j]].ip = deter[node[all_boundary_nodes[i]].n_n[k]].ip;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yim;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zim = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zim;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zim = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zim;
						deter[node[all_boundary_nodes[i]].n_n[j]].im = deter[node[all_boundary_nodes[i]].n_n[k]].im;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjm;
						deter[node[all_boundary_nodes[i]].n_n[j]].jp = deter[node[all_boundary_nodes[i]].n_n[k]].jm;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_yjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_yjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zjp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zjm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zjp;
						deter[node[all_boundary_nodes[i]].n_n[j]].jm = deter[node[all_boundary_nodes[i]].n_n[k]].jp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykp;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkp;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkp = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkp;
						deter[node[all_boundary_nodes[i]].n_n[j]].kp = deter[node[all_boundary_nodes[i]].n_n[k]].kp;

						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_xkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_xkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_ykm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_ykm;
						metric[node[all_boundary_nodes[i]].n_n[j]].zeta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].zeta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].eta_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].eta_zkm;
						metric[node[all_boundary_nodes[i]].n_n[j]].xi_zkm = metric[node[all_boundary_nodes[i]].n_n[k]].xi_zkm;
						deter[node[all_boundary_nodes[i]].n_n[j]].km = deter[node[all_boundary_nodes[i]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zip = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].ip = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].ip;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zim = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].im = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].im;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jm;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].jm = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].jp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].kp = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].kp;

						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].km = deter[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].km;


						//node[g_node-1].n_n[j] = g_node;
						//node[g_node].n_n[k] = g_node-1;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zip;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zip = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zip;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].ip = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].ip;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zim;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zim = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zim;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].im = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].im;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jm;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_yjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_yjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zjp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zjm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zjp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].jm = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].jp;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkm;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkp = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkm;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].kp = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].km;

						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_xkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_xkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_ykm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_ykp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].zeta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].zeta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].eta_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].eta_zkp;
						metric[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].xi_zkm = metric[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].xi_zkp;
						deter[node[node[node[all_boundary_nodes[i]].n_n[j]].n_n[j]].n_n[j]].km = deter[node[node[node[all_boundary_nodes[i]].n_n[k]].n_n[k]].n_n[k]].kp;

						break;
					}
				}
			}
		}

		for (i = 0; i<all_bou_node; i++)
		{
			if (node[all_boundary_nodes[i]].loc >= 27 && node[all_boundary_nodes[i]].loc <= 32)
			{
				metric[all_boundary_nodes[i]].zeta_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_xjp;
				metric[all_boundary_nodes[i]].eta_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_xjp;
				metric[all_boundary_nodes[i]].xi_xjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_xjp;
				metric[all_boundary_nodes[i]].zeta_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_yjp;
				metric[all_boundary_nodes[i]].eta_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_yjp;
				metric[all_boundary_nodes[i]].xi_yjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_yjp;
				metric[all_boundary_nodes[i]].zeta_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].zeta_zjp;
				metric[all_boundary_nodes[i]].eta_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].eta_zjp;
				metric[all_boundary_nodes[i]].xi_zjp = metric[singular[all_boundary_nodes[i]].n_n[2]].xi_zjp;
				deter[all_boundary_nodes[i]].jp = 1.0 / ((metric[all_boundary_nodes[i]].zeta_xjp)*(metric[all_boundary_nodes[i]].eta_yjp) - (metric[all_boundary_nodes[i]].eta_xjp)*(metric[all_boundary_nodes[i]].zeta_yjp));

				metric[all_boundary_nodes[i]].zeta_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_xip;
				metric[all_boundary_nodes[i]].eta_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_xip;
				metric[all_boundary_nodes[i]].xi_xip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_xip;
				metric[all_boundary_nodes[i]].zeta_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_yip;
				metric[all_boundary_nodes[i]].eta_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_yip;
				metric[all_boundary_nodes[i]].xi_yip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_yip;
				metric[all_boundary_nodes[i]].zeta_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].zeta_zip;
				metric[all_boundary_nodes[i]].eta_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].eta_zip;
				metric[all_boundary_nodes[i]].xi_zip = metric[singular[all_boundary_nodes[i]].n_n[3]].xi_zip;
				deter[all_boundary_nodes[i]].ip = 1.0 / ((metric[all_boundary_nodes[i]].zeta_xip)*(metric[all_boundary_nodes[i]].eta_yip) - (metric[all_boundary_nodes[i]].eta_xip)*(metric[all_boundary_nodes[i]].zeta_yip));

				metric[all_boundary_nodes[i]].zeta_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_xjm;
				metric[all_boundary_nodes[i]].eta_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_xjm;
				metric[all_boundary_nodes[i]].xi_xjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_xjm;
				metric[all_boundary_nodes[i]].zeta_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_yjm;
				metric[all_boundary_nodes[i]].eta_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_yjm;
				metric[all_boundary_nodes[i]].xi_yjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_yjm;
				metric[all_boundary_nodes[i]].zeta_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].zeta_zjm;
				metric[all_boundary_nodes[i]].eta_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].eta_zjm;
				metric[all_boundary_nodes[i]].xi_zjm = metric[singular[all_boundary_nodes[i]].n_n[0]].xi_zjm;
				deter[all_boundary_nodes[i]].jm = 1.0 / ((metric[all_boundary_nodes[i]].zeta_xjm)*(metric[all_boundary_nodes[i]].eta_yjm) - (metric[all_boundary_nodes[i]].eta_xjm)*(metric[all_boundary_nodes[i]].zeta_yjm));

			}
			if (node[all_boundary_nodes[i]].loc >= 33 && node[all_boundary_nodes[i]].loc <= 40)
			{
				if (node[all_boundary_nodes[i]].loc == 33)
				{
					gar4 = 2;
				}
				if (node[all_boundary_nodes[i]].loc == 34)
				{
					gar4 = 5;
				}
				if (node[all_boundary_nodes[i]].loc == 35)
				{
					gar4 = 4;
				}
				if (node[all_boundary_nodes[i]].loc == 36)
				{
					gar4 = 6;
				}
				metric[all_boundary_nodes[i]].zeta_xjm = metric[all_boundary_nodes[i]].zeta_xjp;
				metric[all_boundary_nodes[i]].eta_xjm = metric[all_boundary_nodes[i]].eta_xjp;
				metric[all_boundary_nodes[i]].xi_xjm = metric[all_boundary_nodes[i]].xi_xjp;
				metric[all_boundary_nodes[i]].zeta_yjm = metric[all_boundary_nodes[i]].zeta_yjp;
				metric[all_boundary_nodes[i]].eta_yjm = metric[all_boundary_nodes[i]].eta_yjp;
				metric[all_boundary_nodes[i]].xi_yjm = metric[all_boundary_nodes[i]].xi_yjp;
				metric[all_boundary_nodes[i]].zeta_zjm = metric[all_boundary_nodes[i]].zeta_zjp;
				metric[all_boundary_nodes[i]].eta_zjm = metric[all_boundary_nodes[i]].eta_zjp;
				metric[all_boundary_nodes[i]].xi_zjm = metric[all_boundary_nodes[i]].xi_zjp;
				deter[all_boundary_nodes[i]].jm = deter[all_boundary_nodes[i]].jp;
			}
			if (node[all_boundary_nodes[i]].loc >= 41 && node[all_boundary_nodes[i]].loc <= 48)
			{
				metric[all_boundary_nodes[i]].zeta_xjp = metric[all_boundary_nodes[i]].zeta_xjm;
				metric[all_boundary_nodes[i]].eta_xjp = metric[all_boundary_nodes[i]].eta_xjm;
				metric[all_boundary_nodes[i]].xi_xjp = metric[all_boundary_nodes[i]].xi_xjm;
				metric[all_boundary_nodes[i]].zeta_yjp = metric[all_boundary_nodes[i]].zeta_yjm;
				metric[all_boundary_nodes[i]].eta_yjp = metric[all_boundary_nodes[i]].eta_yjm;
				metric[all_boundary_nodes[i]].xi_yjp = metric[all_boundary_nodes[i]].xi_yjm;
				metric[all_boundary_nodes[i]].zeta_zjp = metric[all_boundary_nodes[i]].zeta_zjm;
				metric[all_boundary_nodes[i]].eta_zjp = metric[all_boundary_nodes[i]].eta_zjm;
				metric[all_boundary_nodes[i]].xi_zjp = metric[all_boundary_nodes[i]].xi_zjm;
				deter[all_boundary_nodes[i]].jp = deter[all_boundary_nodes[i]].jm;
			}
			if (node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 35 || node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 50 || node[all_boundary_nodes[i]].loc == 46 || node[all_boundary_nodes[i]].loc == 41 || node[all_boundary_nodes[i]].loc == 45 || node[all_boundary_nodes[i]].loc == 49)
			{
				metric[all_boundary_nodes[i]].zeta_xip = metric[all_boundary_nodes[i]].zeta_xim;
				metric[all_boundary_nodes[i]].eta_xip = metric[all_boundary_nodes[i]].eta_xim;
				metric[all_boundary_nodes[i]].xi_xip = metric[all_boundary_nodes[i]].xi_xim;
				metric[all_boundary_nodes[i]].zeta_yip = metric[all_boundary_nodes[i]].zeta_yim;
				metric[all_boundary_nodes[i]].eta_yip = metric[all_boundary_nodes[i]].eta_yim;
				metric[all_boundary_nodes[i]].xi_yip = metric[all_boundary_nodes[i]].xi_yim;
				metric[all_boundary_nodes[i]].zeta_zip = metric[all_boundary_nodes[i]].zeta_zim;
				metric[all_boundary_nodes[i]].eta_zip = metric[all_boundary_nodes[i]].eta_zim;
				metric[all_boundary_nodes[i]].xi_zip = metric[all_boundary_nodes[i]].xi_zim;
				deter[all_boundary_nodes[i]].ip = deter[all_boundary_nodes[i]].im;
			}
			if (node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 33 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 51 || node[all_boundary_nodes[i]].loc == 47 || node[all_boundary_nodes[i]].loc == 43 || node[all_boundary_nodes[i]].loc == 48 || node[all_boundary_nodes[i]].loc == 52)
			{
				metric[all_boundary_nodes[i]].zeta_xim = metric[all_boundary_nodes[i]].zeta_xip;
				metric[all_boundary_nodes[i]].eta_xim = metric[all_boundary_nodes[i]].eta_xip;
				metric[all_boundary_nodes[i]].xi_xim = metric[all_boundary_nodes[i]].xi_xip;
				metric[all_boundary_nodes[i]].zeta_yim = metric[all_boundary_nodes[i]].zeta_yip;
				metric[all_boundary_nodes[i]].eta_yim = metric[all_boundary_nodes[i]].eta_yip;
				metric[all_boundary_nodes[i]].xi_yim = metric[all_boundary_nodes[i]].xi_yip;
				metric[all_boundary_nodes[i]].zeta_zim = metric[all_boundary_nodes[i]].zeta_zip;
				metric[all_boundary_nodes[i]].eta_zim = metric[all_boundary_nodes[i]].eta_zip;
				metric[all_boundary_nodes[i]].xi_zim = metric[all_boundary_nodes[i]].xi_zip;
				deter[all_boundary_nodes[i]].im = deter[all_boundary_nodes[i]].ip;
			}
			if (node[all_boundary_nodes[i]].loc == 37 || node[all_boundary_nodes[i]].loc == 34 || node[all_boundary_nodes[i]].loc == 40 || node[all_boundary_nodes[i]].loc == 52 || node[all_boundary_nodes[i]].loc == 48 || node[all_boundary_nodes[i]].loc == 44 || node[all_boundary_nodes[i]].loc == 45 || node[all_boundary_nodes[i]].loc == 49)
			{
				metric[all_boundary_nodes[i]].zeta_xkm = metric[all_boundary_nodes[i]].zeta_xkp;
				metric[all_boundary_nodes[i]].eta_xkm = metric[all_boundary_nodes[i]].eta_xkp;
				metric[all_boundary_nodes[i]].xi_xkm = metric[all_boundary_nodes[i]].xi_xkp;
				metric[all_boundary_nodes[i]].zeta_ykm = metric[all_boundary_nodes[i]].zeta_ykp;
				metric[all_boundary_nodes[i]].eta_ykm = metric[all_boundary_nodes[i]].eta_ykp;
				metric[all_boundary_nodes[i]].xi_ykm = metric[all_boundary_nodes[i]].xi_ykp;
				metric[all_boundary_nodes[i]].zeta_zkm = metric[all_boundary_nodes[i]].zeta_zkp;
				metric[all_boundary_nodes[i]].eta_zkm = metric[all_boundary_nodes[i]].eta_zkp;
				metric[all_boundary_nodes[i]].xi_zkm = metric[all_boundary_nodes[i]].xi_zkp;
				deter[all_boundary_nodes[i]].km = deter[all_boundary_nodes[i]].kp;
			}
			if (node[all_boundary_nodes[i]].loc == 38 || node[all_boundary_nodes[i]].loc == 36 || node[all_boundary_nodes[i]].loc == 39 || node[all_boundary_nodes[i]].loc == 51 || node[all_boundary_nodes[i]].loc == 47 || node[all_boundary_nodes[i]].loc == 42 || node[all_boundary_nodes[i]].loc == 46 || node[all_boundary_nodes[i]].loc == 50)
			{
				metric[all_boundary_nodes[i]].zeta_xkp = metric[all_boundary_nodes[i]].zeta_xkm;
				metric[all_boundary_nodes[i]].eta_xkp = metric[all_boundary_nodes[i]].eta_xkm;
				metric[all_boundary_nodes[i]].xi_xkp = metric[all_boundary_nodes[i]].xi_xkm;
				metric[all_boundary_nodes[i]].zeta_ykp = metric[all_boundary_nodes[i]].zeta_ykm;
				metric[all_boundary_nodes[i]].eta_ykp = metric[all_boundary_nodes[i]].eta_ykm;
				metric[all_boundary_nodes[i]].xi_ykp = metric[all_boundary_nodes[i]].xi_ykm;
				metric[all_boundary_nodes[i]].zeta_zkp = metric[all_boundary_nodes[i]].zeta_zkm;
				metric[all_boundary_nodes[i]].eta_zkp = metric[all_boundary_nodes[i]].eta_zkm;
				metric[all_boundary_nodes[i]].xi_zkp = metric[all_boundary_nodes[i]].xi_zkm;
				deter[all_boundary_nodes[i]].kp = deter[all_boundary_nodes[i]].km;
			}

		}


		/*******************************************************malloc for variables********************************************/
		u.resize(5);
		v.resize(5);
		w.resize(5);
		rho.resize(5);
		p.resize(5);
		t.resize(5);
		mu.resize(5);
		a.resize(5);
		e.resize(5);

		for (i = 0; i < 5; i++)
		{
			u[i].resize(nodes);
			v[i].resize(nodes);
			w[i].resize(nodes);
			rho[i].resize(nodes);
			p[i].resize(nodes);
			t[i].resize(nodes);
			mu[i].resize(nodes);
			a[i].resize(nodes);
			e[i].resize(nodes);
		}

		tauzz.resize(nodes);
		tauee.resize(nodes);
		tauxx.resize(nodes);
		tauze.resize(nodes);
		tauzx.resize(nodes);
		tauex.resize(nodes);

		TAU_SGS_XX.resize(nodes);
		TAU_SGS_YY.resize(nodes);
		TAU_SGS_ZZ.resize(nodes);
		TAU_SGS_XY.resize(nodes);
		TAU_SGS_XZ.resize(nodes);
		TAU_SGS_YZ.resize(nodes);
		H_SGS_X.resize(nodes);
		H_SGS_Y.resize(nodes);
		H_SGS_Z.resize(nodes);
		D_SGS_X.resize(nodes);
		D_SGS_Y.resize(nodes);
		D_SGS_Z.resize(nodes);
		DUCROS.resize(nodes);

		qz.resize(nodes);
		qe.resize(nodes);
		qx.resize(nodes);

		del_cfl.resize(nodes);
		v_dash1.resize(nodes);


		/*******************************malloc for solver variables********************************************/
		roe_rho_ip.resize(nodes);
		roe_u_ip.resize(nodes);
		roe_v_ip.resize(nodes);
		roe_w_ip.resize(nodes);
		roe_h_ip.resize(nodes);

		final_U.resize(5);

		roe_rho_im.resize(nodes);
		roe_u_im.resize(nodes);
		roe_v_im.resize(nodes);
		roe_w_im.resize(nodes);
		roe_h_im.resize(nodes);

		roe_rho_jp.resize(nodes);
		roe_u_jp.resize(nodes);
		roe_v_jp.resize(nodes);
		roe_w_jp.resize(nodes);
		roe_h_jp.resize(nodes);

		roe_rho_jm.resize(nodes);
		roe_u_jm.resize(nodes);
		roe_v_jm.resize(nodes);
		roe_w_jm.resize(nodes);
		roe_h_jm.resize(nodes);

		roe_rho_kp.resize(nodes);
		roe_u_kp.resize(nodes);
		roe_v_kp.resize(nodes);
		roe_w_kp.resize(nodes);
		roe_h_kp.resize(nodes);

		roe_rho_km.resize(nodes);
		roe_u_km.resize(nodes);
		roe_v_km.resize(nodes);
		roe_w_km.resize(nodes);
		roe_h_km.resize(nodes);
		
		roe_R.resize(nodes);
		roe_a.resize(nodes);

		diver.resize(nodes);

		Qi_iminus.resize(nodes);
		Qi_iplus.resize(nodes);
		Qj_iminus.resize(nodes);
		Qj_iplus.resize(nodes);
		Qk_iminus.resize(nodes);
		Qk_iplus.resize(nodes);

		Qi_iminus_n.resize(nodes);
		Qi_iplus_n.resize(nodes);
		Qj_iminus_n.resize(nodes);
		Qj_iplus_n.resize(nodes);
		Qk_iminus_n.resize(nodes);
		Qk_iplus_n.resize(nodes);

		w_Qip.resize(5);
		w_Qinp.resize(5);
		w_Qim.resize(5);
		w_Qinm.resize(5);
		w_Qjp.resize(5);
		w_Qjnp.resize(5);
		w_Qjm.resize(5);
		w_Qjnm.resize(5);
		w_Qkp.resize(5);
		w_Qknp.resize(5);
		w_Qkm.resize(5);
		w_Qknm.resize(5);


		W_Qip.resize(5);
		W_Qinp.resize(5);
		W_Qim.resize(5);
		W_Qinm.resize(5);
		W_Qjp.resize(5);
		W_Qjnp.resize(5);
		W_Qjm.resize(5);
		W_Qjnm.resize(5);
		W_Qkp.resize(5);
		W_Qknp.resize(5);
		W_Qkm.resize(5);
		W_Qknm.resize(5);

		L.resize(nodes);
		U.resize(nodes);
		E.resize(nodes);
		F.resize(nodes);
		G.resize(nodes);
		Ev.resize(nodes);
		Fv.resize(nodes);
		Gv.resize(nodes);
		E1.resize(nodes);
		F1.resize(nodes);
		G1.resize(nodes);
		Ev1.resize(nodes);
		Fv1.resize(nodes);
		Gv1.resize(nodes);
		dF.resize(nodes);
		dE.resize(nodes);
		dG.resize(nodes);
		dFv.resize(nodes);
		dEv.resize(nodes);
		dGv.resize(nodes);


		Qi_half_m.resize(5);
		Qi_halfn_m.resize(5);
		Qi_half_p.resize(5);
		Qi_half_np.resize(5);
		Qj_half_m.resize(5);
		Qj_halfn_m.resize(5);
		Qj_half_p.resize(5);
		Qj_half_np.resize(5);
		Qk_half_m.resize(5);
		Qk_halfn_m.resize(5);
		Qk_half_p.resize(5);
		Qk_half_np.resize(5);

		IS_Qim.resize(5);
		IS_Qinm.resize(5);
		IS_Qip.resize(5);
		IS_Qinp.resize(5);
		IS_Qjm.resize(5);
		IS_Qjnm.resize(5);
		IS_Qjp.resize(5);
		IS_Qjnp.resize(5);
		IS_Qkm.resize(5);
		IS_Qknm.resize(5);
		IS_Qkp.resize(5);
		IS_Qknp.resize(5);

		r_eigen_Qip.resize(5);
		r_eigen_Qjp.resize(5);
		r_eigen_Qkp.resize(5);
		l_eigen_Qip.resize(5);
		l_eigen_Qjp.resize(5);
		l_eigen_Qkp.resize(5);

		r_eigen_Qim.resize(5);
		r_eigen_Qjm.resize(5);
		r_eigen_Qkm.resize(5);
		l_eigen_Qim.resize(5);
		l_eigen_Qjm.resize(5);
		l_eigen_Qkm.resize(5);

		for (i = 0; i<nodes; i++)
		{

			Qi_iminus[i].resize(5);
			Qi_iplus[i].resize(5);
			Qj_iminus[i].resize(5);
			Qj_iplus[i].resize(5);
			Qk_iminus[i].resize(5);
			Qk_iplus[i].resize(5);

			diver[i].resize(10);

			Qi_iminus_n[i].resize(5);
			Qi_iplus_n[i].resize(5);
			Qj_iminus_n[i].resize(5);
			Qj_iplus_n[i].resize(5);
			Qk_iminus_n[i].resize(5);
			Qk_iplus_n[i].resize(5);

			L[i].resize(5);
			U[i].resize(5);
			E[i].resize(5);
			F[i].resize(5);
			G[i].resize(5);
			Ev[i].resize(5);
			Fv[i].resize(5);
			Gv[i].resize(5);
			E1[i].resize(5);
			F1[i].resize(5);
			G1[i].resize(5);
			Ev1[i].resize(5);
			Fv1[i].resize(5);
			Gv1[i].resize(5);
		}

		for (i = 0; i<5; i++)
		{
			w_Qip[i].resize(5);
			w_Qinp[i].resize(5);
			w_Qim[i].resize(5);
			w_Qinm[i].resize(5);
			w_Qjp[i].resize(5);
			w_Qjnp[i].resize(5);
			w_Qjm[i].resize(5);
			w_Qjnm[i].resize(5);
			w_Qkp[i].resize(5);
			w_Qknp[i].resize(5);
			w_Qkm[i].resize(5);
			w_Qknm[i].resize(5);

			W_Qip[i].resize(5);
			W_Qinp[i].resize(5);
			W_Qim[i].resize(5);
			W_Qinm[i].resize(5);
			W_Qjp[i].resize(5);
			W_Qjnp[i].resize(5);
			W_Qjm[i].resize(5);
			W_Qjnm[i].resize(5);
			W_Qkp[i].resize(5);
			W_Qknp[i].resize(5);
			W_Qkm[i].resize(5);
			W_Qknm[i].resize(5);

			IS_Qim[i].resize(5);
			IS_Qinm[i].resize(5);
			IS_Qip[i].resize(5);
			IS_Qinp[i].resize(5);
			IS_Qjm[i].resize(5);
			IS_Qjnm[i].resize(5);
			IS_Qjp[i].resize(5);
			IS_Qjnp[i].resize(5);
			IS_Qkm[i].resize(5);
			IS_Qknm[i].resize(5);
			IS_Qkp[i].resize(5);
			IS_Qknp[i].resize(5);
			Qi_half_m[i].resize(5);
			Qi_halfn_m[i].resize(5);
			Qi_half_p[i].resize(5);
			Qi_half_np[i].resize(5);
			Qj_half_m[i].resize(5);
			Qj_halfn_m[i].resize(5);
			Qj_half_p[i].resize(5);
			Qj_half_np[i].resize(5);
			Qk_half_m[i].resize(5);
			Qk_halfn_m[i].resize(5);
			Qk_half_p[i].resize(5);
			Qk_half_np[i].resize(5);

			r_eigen_Qip[i].resize(5);
			r_eigen_Qjp[i].resize(5);
			r_eigen_Qkp[i].resize(5);
			l_eigen_Qip[i].resize(5);
			l_eigen_Qjp[i].resize(5);
			l_eigen_Qkp[i].resize(5);
			l_eigen_Qim[i].resize(5);
			l_eigen_Qjm[i].resize(5);
			l_eigen_Qkm[i].resize(5);

			r_eigen_Qim[i].resize(5);
			r_eigen_Qjm[i].resize(5);
			r_eigen_Qkm[i].resize(5);
		}

		for (i = 0; i<nodes; i++)
		{
			for (j = 0; j<5; j++)
			{
				L[i][j].resize(5);
				U[i][j].resize(5);
			}
		}
		
		temp_a.resize(sd_node + 10);

		/******************************Reading restart file**********************************************/
		sprintf(filename2, "restart_file_%d.neu", restart_num);
		inFile.open(filename2);
		if (inFile.good())
		{
			i = 1;
			while (inFile.good())
			{
				getline(inFile, line, '\n');
				istringstream streamA(line);
				for (k = 1; k <= 9; k++)
				{
					if (k == 1)
						streamA >> u[0][d_node[i].local];
					else if (k == 2)
						streamA >> v[0][d_node[i].local];
					else if (k == 3)
						streamA >> w[0][d_node[i].local];
					else if (k == 4)
						streamA >> rho[0][d_node[i].local];
					else if (k == 5)
						streamA >> p[0][d_node[i].local];
					else if (k == 6)
						streamA >> t[0][d_node[i].local];
					else if (k == 7)
						streamA >> e[0][d_node[i].local];
					else if (k == 8)
						streamA >> mu[0][d_node[i].local];
					else if (k == 9)
						streamA >> itera;
				}
				i++;
			}
			itera++;
			itera++;
			inFile.close();
		}
		else
		{
			itera = 1;
		}

		d_node.clear();
		/***************************SENDING LOCAL NODE DATA TO MASTER************************************/
		MPI_Send(&sd_node, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD);

		position = 0;
		MPI_Pack_size(sd_node, MPI_INT, MPI_COMM_WORLD, &memsize);
		buffer = new char[memsize];
		for (i = 1; i <= sd_node; i++)
		{
			temp_a[i] = node[i].global;
			MPI_Pack(&temp_a[i], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
		}
		MPI_Send(buffer, memsize, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		delete[] buffer;

		/*************************************************************************************************
		Flow variables initializing
		*************************************************************************************************/
		int back_p;
		j = 0;
		initial = 0;

		if (itera != 1)
		{
			initial = 1;
		}
		intialise(j);

		temp1 = sd_node;
		position = 0;
		for (i = 0; i<neigh_pro; i++)
		{
			MPI_Pack_size(recv_c[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize1);
			buffer2 = new char[memsize1];                                          /***********carefull with buffer1 and buffer2******************/
			MPI_Pack_size(proc_node[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			position = 0;
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&u[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&v[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&w[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&p[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&t[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&rho[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&e[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}
			for (j = 0; j<recv_c[c[i]]; j++)
			{
				MPI_Pack(&mu[0][loc_dat[i][j]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
			}

			MPI_Sendrecv(buffer2, memsize1, MPI_PACKED, c[i], c[i], buffer, memsize, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
			position = 0;
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &u[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &v[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &w[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &p[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &t[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &rho[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &e[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			for (j = 0; j<proc_node[c[i]]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &mu[0][recv_b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
			}
			delete[] buffer;
			delete[] buffer2;
		}

		j = 0;
		initial = 1;
		intialise(j);

		sd_node = temp1;
		for (k = 1; k <= no_of_tip_send; k++)
		{
			u[0][tip[k].node] = 0.0;
			v[0][tip[k].node] = 0.0;
			w[0][tip[k].node] = 0.0;
			p[0][tip[k].node] = p[0][node[tip[k].node].n_n[3]];
			t[0][tip[k].node] = t[0][node[tip[k].node].n_n[3]];
			rho[0][tip[k].node] = (1.4*Mach*Mach)*(p[0][tip[k].node] / t[0][tip[k].node]);
			e[0][tip[k].node] = p[0][tip[k].node] / (0.4*rho[0][tip[k].node]);
		}

		for (k = 1; k <= no_of_tip_recv; k++)
		{
			u[0][tip_recv[k].node] = 0.0;
			v[0][tip_recv[k].node] = 0.0;
			w[0][tip_recv[k].node] = 0.0;
			p[0][tip_recv[k].node] = p[0][node[tip_recv[k].node].n_n[3]];
			t[0][tip_recv[k].node] = t[0][node[tip_recv[k].node].n_n[3]];
			rho[0][tip_recv[k].node] = (1.4*Mach*Mach)*(p[0][tip_recv[k].node] / t[0][tip_recv[k].node]);
			e[0][tip_recv[k].node] = p[0][tip_recv[k].node] / (0.4*rho[0][tip_recv[k].node]);
		}

		j = 0;
		sd_node = temp1;

		if (itera == 1)
		{
			position = 0;
			MPI_Pack_size(sd_node * 10, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Pack(&u[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&v[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&w[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&rho[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&p[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&t[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&e[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&mu[0][1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
			MPI_Pack(&DUCROS[1], sd_node, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);

			MPI_Send(buffer, memsize, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
			delete[] buffer;
		}

		MPI_Barrier(comm_slaves);

		/***********************************************************************************************************************************
		SOLVER LOOP
		************************************************************************************************************************************/
		solver();
	}

	/**********************DATA GATHERING AND FILE WRITING AND PRINTING AFTER ITERATIONS************************************************/
	if (myrank == 0)
	{
		printf("\nfreeing memory of temp arrays\n");
		
		max_div_u.resize(size);
		max_div_v.resize(size);
		max_div_e.resize(size);
		sd_count.resize(size);
		b.resize(size);
		sd_gh_list.resize(size);

		u.resize(5);
		v.resize(5);
		w.resize(5);
		rho.resize(5);
		p.resize(5);
		t.resize(5);
		mu.resize(5);
		e.resize(5);
		DUCROS.resize(nodes);

		for (i = 0; i<5; i++)
		{
			u[i].resize(nodes);
			v[i].resize(nodes);
			w[i].resize(nodes);
			rho[i].resize(nodes);
			p[i].resize(nodes);
			t[i].resize(nodes);
			mu[i].resize(nodes);
			e[i].resize(nodes);
		}

		/*****************************************************************************************************************/
		
		neigh_pro_list.resize(size + 1);
		all_c_list.resize(size + 1);
		
		for (i = 1; i <= size - 1; i++)
		{
			MPI_Recv(&neigh_pro_list[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);

			all_c_list[i].resize(neigh_pro_list[i] + 1);
			position = 0;
			MPI_Pack_size(neigh_pro_list[i], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);
			for (j = 1; j <= neigh_pro_list[i]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &glob, 1, MPI_INT, MPI_COMM_WORLD);
				all_c_list[i][j] = glob;
			}
			delete[] buffer;
		}



		for (i = 1; i <= size - 1; i++)
		{
			for (k = 1; k <= size - 1; k++)
			{
				MPI_Send(&neigh_pro_list[k], 1, MPI_INT, i, i, MPI_COMM_WORLD);
				position = 0;
				MPI_Pack_size(neigh_pro_list[k], MPI_INT, MPI_COMM_WORLD, &memsize);
				buffer = new char[memsize];
				for (j = 1; j <= neigh_pro_list[k]; j++)
				{
					MPI_Pack(&all_c_list[k][j], 1, MPI_INT, buffer, memsize, &position, MPI_COMM_WORLD);
				}

				MPI_Send(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD);
				delete[] buffer;
			}
		}

		/********************************************************************************************************************/

		for (i = 1; i <= size - 1; i++)
		{
			MPI_Recv(&sd_gh_list[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
			position = 0;
			MPI_Pack_size(sd_gh_list[i], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);
			for (j = 1; j< sd_gh_list[i]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &glob, 1, MPI_INT, MPI_COMM_WORLD);
				node[glob].req = i;
			}
			delete[] buffer;
		}

		
		jacobian.clear();
		metric.clear();
		det.clear();
		deter.clear();

		
		/*********************************************************/
		proc_request();


		for (i = 1; i <= size - 1; i++)
		{
			MPI_Recv(&sd_count[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
			position = 0;
			MPI_Pack_size(sd_count[i], MPI_INT, MPI_COMM_WORLD, &memsize);
			buffer = new char[memsize];
			MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);
			b[i].resize(sd_count[i] + 10);
			for (j = 1; j <= sd_count[i]; j++)
			{
				MPI_Unpack(buffer, memsize, &position, &glob, 1, MPI_INT, MPI_COMM_WORLD);
				b[i][j] = glob;
			}
			delete[] buffer;
		}

		ifstream inFile;
		sprintf(filename2, "restart_file_%d.neu", restart_num);
		inFile.open(filename2);
		if (inFile.good())
		{
			getline(inFile, line2, '\n');
			istringstream streamA(line2);
			for (k = 1; k <= 9; k++)
			{
				if (k<9)
				{
					streamA >> temp_u;
				}
				if (k == 9)
				{
					streamA >> itera;
				}
			}
			inFile.close();
			itera++;
		}
		else
		{
			itera = 0;
		}

		printf("deleting all temperory files\n");
		system("rm -rf *.txt *.o node_neighbour.neu elem_neighbour.neu");
		printf("calculations started\n");

		restart_num = 0;
		for (iter = itera; iter<iterations; iter++)
		{
			if (iter > itera)
			{
				for (i = 1; i <= size - 1; i++)
				{
					position = 0;
					MPI_Pack_size(3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
					buffer = new char[memsize];
					MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);
					max_div_u[i] = temp_d;
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);
					max_div_v[i] = temp_d;
					MPI_Unpack(buffer, memsize, &position, &temp_d, 1, MPI_DOUBLE, MPI_COMM_WORLD);
					max_div_e[i] = temp_d;

					delete[] buffer;
				}

				for (i = 1; i <= size - 1; i++)
				{
					for (i = 1; i <= size - 1; i++)
					{
						if (max_div_u[i] < max_div_u[i + 1])
						{
							temp_d = max_div_u[i + 1];
							max_div_u[i + 1] = max_div_u[i];
							max_div_u[i] = temp_d;
						}
						if (max_div_v[i] < max_div_v[i + 1])
						{
							temp_d = max_div_v[i + 1];
							max_div_v[i + 1] = max_div_v[i];
							max_div_v[i] = temp_d;
						}
						if (max_div_e[i] < max_div_e[i + 1])
						{
							temp_d = max_div_e[i + 1];
							max_div_e[i + 1] = max_div_e[i];
							max_div_e[i] = temp_d;
						}
					}
				}
				cout << "iteration = " << iter - 1 << " max_div_u = " << max_div_u[1] << " max_div_v = " << max_div_v[1] << " max_div_e = " << max_div_e[1] << "\n";
			}

			if (iter % 100 == 0)
			{
				restart_num++;
				for (i = 1; i <= size - 1; i++)
				{
					position = 0;
					MPI_Pack_size(sd_count[i] * 10, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
					buffer = new char[memsize];
					temp = sd_count[i] * 8;

					MPI_Recv(buffer, memsize, MPI_PACKED, i, i, MPI_COMM_WORLD, &status);
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &u[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &v[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &w[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &rho[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &p[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &t[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &e[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &mu[0][b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					for (j = 1; j <= sd_count[i]; j++)
					{
						MPI_Unpack(buffer, memsize, &position, &DUCROS[b[i][j]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
					}
					delete[] buffer;
				}
				if (iter % 100 == 0)
				{
					writepltfile();

					sprintf(filename2, "nohup ./preplot nodefile_%d.dat > preplot.out", iter);
					system(filename2);
					system("rm -rf preplot.out *.dat");
				}

				if (restart_num > 2)
				{
					restart_num = 1;
				}

				write_restart_file();
			}
		}
	}
}

