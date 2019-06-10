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

void solver()
{
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	int k, m, weno_exec, gar, lm, value, u_loc, v_loc, e_loc;
	double epsi, rc;
	int position, memsize, memsize1;
	char *buffer, *buffer2;

	
	e_inf = e[0][10650];
	epsi = 0.003;

	rc = 0.3;
	epsi = 0.003;

	max_div_um = -1e25;
	max_div_vm = -1e25;
	max_div_wm = -1e25;
	max_div_em = -1e25;

	alpha_u.resize(5);
	alpha_v.resize(5);
	alpha_w.resize(5);
	alpha_u_ip.resize(5);
	alpha_u_im.resize(5);
	alpha_v_ip.resize(5);
	alpha_v_im.resize(5);
	alpha_w_ip.resize(5);
	alpha_w_im.resize(5);
	alpha_u_jp.resize(5);
	alpha_u_jm.resize(5);
	alpha_v_jp.resize(5);
	alpha_v_jm.resize(5);
	alpha_w_jp.resize(5);
	alpha_w_jm.resize(5);
	alpha_u_kp.resize(5);
	alpha_u_km.resize(5);
	alpha_v_kp.resize(5);
	alpha_v_km.resize(5);
	alpha_w_kp.resize(5);
	alpha_w_km.resize(5);


	lm = 0;
	for (iter = itera; iter<iterations; iter++)
	{
		max_div_um = -1e25;
		max_div_vm = -1e25;
		max_div_em = -1e25;
		gar = 0;

		weno_exec = 1;

		for (j = 0; j<3; j++)
		{
			lm++;
			gar++;
			for (i = 1; i<g_node; i++)
			{
				roe_average();

				U[i][j][0] = (det[i])*rho[j][i];
				U[i][j][1] = (det[i])*rho[j][i] * u[j][i];
				U[i][j][2] = (det[i])*rho[j][i] * v[j][i];
				U[i][j][3] = (det[i])*rho[j][i] * w[j][i];
				U[i][j][4] = (det[i])*((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]));

				Q[i][0].ip = rho[j][i];
				Q[i][1].ip = rho[j][i] * u[j][i];
				Q[i][2].ip = rho[j][i] * v[j][i];
				Q[i][3].ip = rho[j][i] * w[j][i];
				Q[i][4].ip = ((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]));
			}

			for (i = 1; i<g_node; i++)
			{
				viscousflux_variables();

				F1[i][0] = rho[j][i] * u[j][i];
				F1[i][1] = rho[j][i] * u[j][i] * u[j][i] + p[j][i];
				F1[i][2] = rho[j][i] * v[j][i] * u[j][i];
				F1[i][3] = rho[j][i] * u[j][i] * w[j][i];
				F1[i][4] = (((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i])) + p[j][i])*u[j][i];

				E1[i][0] = rho[j][i] * v[j][i];
				E1[i][1] = rho[j][i] * v[j][i] * u[j][i];
				E1[i][2] = rho[j][i] * v[j][i] * v[j][i] + p[j][i];
				E1[i][3] = rho[j][i] * v[j][i] * w[j][i];
				E1[i][4] = (((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i])) + p[j][i])*v[j][i];

				G1[i][0] = rho[j][i] * w[j][i];
				G1[i][1] = rho[j][i] * u[j][i] * w[j][i];
				G1[i][2] = rho[j][i] * v[j][i] * w[j][i];
				G1[i][3] = rho[j][i] * w[j][i] * w[j][i] + p[j][i];
				G1[i][4] = (((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i])) + p[j][i])*w[j][i];

				Fv1[i][0] = 0.0;
				Fv1[i][1] = (-1.0)*tauzz[i] + TAU_SGS_XX[i];
				Fv1[i][2] = (-1.0)*tauze[i] + TAU_SGS_XY[i];
				Fv1[i][3] = (-1.0)*tauzx[i] + TAU_SGS_XZ[i];
				Fv1[i][4] = (((-1.0)*u[j][i] * tauzz[i] - v[j][i] * tauze[i] - w[j][i] * tauzx[i]) + qz[i] + H_SGS_X[i] + D_SGS_X[i]);

				Ev1[i][0] = 0.0;
				Ev1[i][1] = (-1.0)*tauze[i] + TAU_SGS_XY[i];
				Ev1[i][2] = (-1.0)*tauee[i] + TAU_SGS_YY[i];
				Ev1[i][3] = (-1.0)*tauex[i] + TAU_SGS_YZ[i];
				Ev1[i][4] = (((-1.0)*u[j][i] * tauze[i] - v[j][i] * tauee[i] - w[j][i] * tauex[i]) + qe[i] + H_SGS_Y[i] + D_SGS_Y[i]);

				Gv1[i][0] = 0.0;
				Gv1[i][1] = (-1.0)*tauzx[i] + TAU_SGS_XZ[i];
				Gv1[i][2] = (-1.0)*tauex[i] + TAU_SGS_YZ[i];
				Gv1[i][3] = (-1.0)*tauxx[i] + TAU_SGS_ZZ[i];
				Gv1[i][4] = (((-1.0)*u[j][i] * tauzx[i] - v[j][i] * tauex[i] - w[j][i] * tauxx[i]) + qx[i] + H_SGS_Z[i] + D_SGS_Z[i]);

				for (k = 0; k<5; k++)
				{
					F[i][k] = (det[i])*(F1[i][k] * metric[i].zeta_x + E1[i][k] * metric[i].zeta_y + G1[i][k] * metric[i].zeta_z);
					E[i][k] = (det[i])*(F1[i][k] * metric[i].eta_x + E1[i][k] * metric[i].eta_y + G1[i][k] * metric[i].eta_z);
					G[i][k] = (det[i])*(F1[i][k] * metric[i].xi_x + E1[i][k] * metric[i].xi_y + G1[i][k] * metric[i].xi_z);

					Fv[i][k] = (det[i])*(Fv1[i][k] * metric[i].zeta_x + Ev1[i][k] * metric[i].zeta_y + Gv1[i][k] * metric[i].zeta_z);
					Ev[i][k] = (det[i])*(Fv1[i][k] * metric[i].eta_x + Ev1[i][k] * metric[i].eta_y + Gv1[i][k] * metric[i].eta_z);
					Gv[i][k] = (det[i])*(Fv1[i][k] * metric[i].xi_x + Ev1[i][k] * metric[i].xi_y + Gv1[i][k] * metric[i].xi_z);
				}
				node[i].val = 0;
			}
			/********************************************************************************************************************************************************************************************/
			for (i = 1; i <= sd_node; i++)
			{
				//	if (node[i].loc == 0 && node[i].corner_ID == 0 )
				{
					value = node[i].val;

					weno_solver_1();
					weno_solver_3();
					weno_solver_0();
					weno_solver_2();
					weno_solver_4();
					weno_solver_5();

					for (k = 0; k<5; k++)
					{
						/** Qi(i+1/2)_plus**/
						Qi_half_p[0][k] = (15.0 / 8.0)*Qi_iplus[node[i].n_n[1]][k] - (5.0 / 4.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k] + (3.0 / 8.0)*Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k];
						Qi_half_p[1][k] = (3.0 / 8.0)*Qi_iplus[i][k] + (3.0 / 4.0)*Qi_iplus[node[i].n_n[1]][k] - (1.0 / 8.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k];
						Qi_half_p[2][k] = (-1.0 / 8.0)*Qi_iplus[node[i].n_n[3]][k] + (3.0 / 4.0)*Qi_iplus[i][k] + (3.0 / 8.0)*Qi_iplus[node[i].n_n[1]][k];

						/** Qi(i+1/2)_minus**/
						Qi_half_m[0][k] = (3.0 / 8.0)*Qi_iplus[i][k] + (3.0 / 4.0)*Qi_iplus[node[i].n_n[1]][k] - (1.0 / 8.0)*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k];
						Qi_half_m[1][k] = (-1.0 / 8.0)*Qi_iplus[node[i].n_n[3]][k] + (3.0 / 4.0)*Qi_iplus[i][k] + (3.0 / 8.0)*Qi_iplus[node[i].n_n[1]][k];
						Qi_half_m[2][k] = (3.0 / 8.0)*Qi_iplus[node[node[i].n_n[3]].n_n[3]][k] - (5.0 / 4.0)*Qi_iplus[node[i].n_n[3]][k] + (15.0 / 8.0)*Qi_iplus[i][k];

						/** Qi(i-1/2)_minus**/
						Qi_halfn_m[0][k] = (3.0 / 8.0)*Qi_iminus_n[node[i].n_n[3]][k] + (3.0 / 4.0)*Qi_iminus_n[i][k] - (1.0 / 8.0)*Qi_iminus_n[node[i].n_n[1]][k];
						Qi_halfn_m[1][k] = (-1.0 / 8.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] + (3.0 / 4.0)*Qi_iminus_n[node[i].n_n[3]][k] + (3.0 / 8.0)*Qi_iminus_n[i][k];
						Qi_halfn_m[2][k] = (3.0 / 8.0)*Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] - (5.0 / 4.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] + (15.0 / 8.0)*Qi_iminus_n[node[i].n_n[3]][k];

						/** Qi(i-1/2)_plus**/
						Qi_half_np[0][k] = (15.0 / 8.0)*Qi_iminus_n[i][k] - (5.0 / 4.0)*Qi_iminus_n[node[i].n_n[1]][k] + (3.0 / 8.0)*Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k];
						Qi_half_np[1][k] = (3.0 / 8.0)*Qi_iminus_n[node[i].n_n[3]][k] + (3.0 / 4.0)*Qi_iminus_n[i][k] - (1.0 / 8.0)*Qi_iminus_n[node[i].n_n[1]][k];
						Qi_half_np[2][k] = (-1.0 / 8.0)*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] + (3.0 / 4.0)*Qi_iminus_n[node[i].n_n[3]][k] + (3.0 / 8.0)*Qi_iminus_n[i][k];

						/** Qj(i+1/2)_plus**/
						Qj_half_p[0][k] = (15.0 / 8.0)*Qj_iplus[node[i].n_n[0]][k] - (5.0 / 4.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k] + (3.0 / 8.0)*Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k];
						Qj_half_p[1][k] = (3.0 / 8.0)*Qj_iplus[i][k] + (3.0 / 4.0)*Qj_iplus[node[i].n_n[0]][k] - (1.0 / 8.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k];
						Qj_half_p[2][k] = (-1.0 / 8.0)*Qj_iplus[node[i].n_n[2]][k] + (3.0 / 4.0)*Qj_iplus[i][k] + (3.0 / 8.0)*Qj_iplus[node[i].n_n[0]][k];

						/** Qj(i+1/2)_minus**/
						Qj_half_m[0][k] = (3.0 / 8.0)*Qj_iplus[i][k] + (3.0 / 4.0)*Qj_iplus[node[i].n_n[0]][k] - (1.0 / 8.0)*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k];
						Qj_half_m[1][k] = (-1.0 / 8.0)*Qj_iplus[node[i].n_n[2]][k] + (3.0 / 4.0)*Qj_iplus[i][k] + (3.0 / 8.0)*Qj_iplus[node[i].n_n[0]][k];
						Qj_half_m[2][k] = (3.0 / 8.0)*Qj_iplus[node[node[i].n_n[2]].n_n[2]][k] - (5.0 / 4.0)*Qj_iplus[node[i].n_n[2]][k] + (15.0 / 8.0)*Qj_iplus[i][k];

						/** Qj(i-1/2)_minus**/
						Qj_halfn_m[0][k] = (3.0 / 8.0)*Qj_iminus_n[node[i].n_n[2]][k] + (3.0 / 4.0)*Qj_iminus_n[i][k] - (1.0 / 8.0)*Qj_iminus_n[node[i].n_n[0]][k];
						Qj_halfn_m[1][k] = (-1.0 / 8.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] + (3.0 / 4.0)*Qj_iminus_n[node[i].n_n[2]][k] + (3.0 / 8.0)*Qj_iminus_n[i][k];
						Qj_halfn_m[2][k] = (3.0 / 8.0)*Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] - (5.0 / 4.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] + (15.0 / 8.0)*Qj_iminus_n[node[i].n_n[2]][k];

						/** Qj(i-1/2)_plus**/
						Qj_half_np[0][k] = (15.0 / 8.0)*Qj_iminus_n[i][k] - (5.0 / 4.0)*Qj_iminus_n[node[i].n_n[0]][k] + (3.0 / 8.0)*Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k];
						Qj_half_np[1][k] = (3.0 / 8.0)*Qj_iminus_n[node[i].n_n[2]][k] + (3.0 / 4.0)*Qj_iminus_n[i][k] - (1.0 / 8.0)*Qj_iminus_n[node[i].n_n[0]][k];
						Qj_half_np[2][k] = (-1.0 / 8.0)*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] + (3.0 / 4.0)*Qj_iminus_n[node[i].n_n[2]][k] + (3.0 / 8.0)*Qj_iminus_n[i][k];

						/** Qk(i+1/2)_plus**/
						Qk_half_p[0][k] = (15.0 / 8.0)*Qk_iplus[node[i].n_n[4]][k] - (5.0 / 4.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k] + (3.0 / 8.0)*Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k];
						Qk_half_p[1][k] = (3.0 / 8.0)*Qk_iplus[i][k] + (3.0 / 4.0)*Qk_iplus[node[i].n_n[4]][k] - (1.0 / 8.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k];
						Qk_half_p[2][k] = (-1.0 / 8.0)*Qk_iplus[node[i].n_n[5]][k] + (3.0 / 4.0)*Qk_iplus[i][k] + (3.0 / 8.0)*Qk_iplus[node[i].n_n[4]][k];

						/** Qk(i+1/2)_minus**/
						Qk_half_m[0][k] = (3.0 / 8.0)*Qk_iplus[i][k] + (3.0 / 4.0)*Qk_iplus[node[i].n_n[4]][k] - (1.0 / 8.0)*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k];
						Qk_half_m[1][k] = (-1.0 / 8.0)*Qk_iplus[node[i].n_n[5]][k] + (3.0 / 4.0)*Qk_iplus[i][k] + (3.0 / 8.0)*Qk_iplus[node[i].n_n[4]][k];
						Qk_half_m[2][k] = (3.0 / 8.0)*Qk_iplus[node[node[i].n_n[5]].n_n[5]][k] - (5.0 / 4.0)*Qk_iplus[node[i].n_n[5]][k] + (15.0 / 8.0)*Qk_iplus[i][k];

						/** Qk(i-1/2)_minus**/
						Qk_halfn_m[0][k] = (3.0 / 8.0)*Qk_iminus_n[node[i].n_n[5]][k] + (3.0 / 4.0)*Qk_iminus_n[i][k] - (1.0 / 8.0)*Qk_iminus_n[node[i].n_n[4]][k];
						Qk_halfn_m[1][k] = (-1.0 / 8.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] + (3.0 / 4.0)*Qk_iminus_n[node[i].n_n[5]][k] + (3.0 / 8.0)*Qk_iminus_n[i][k];
						Qk_halfn_m[2][k] = (3.0 / 8.0)*Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] - (5.0 / 4.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] + (15.0 / 8.0)*Qk_iminus_n[node[i].n_n[5]][k];

						/** Qk(i-1/2)_plus**/
						Qk_half_np[0][k] = (15.0 / 8.0)*Qk_iminus_n[i][k] - (5.0 / 4.0)*Qk_iminus_n[node[i].n_n[4]][k] + (3.0 / 8.0)*Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k];
						Qk_half_np[1][k] = (3.0 / 8.0)*Qk_iminus_n[node[i].n_n[5]][k] + (3.0 / 4.0)*Qk_iminus_n[i][k] - (1.0 / 8.0)*Qk_iminus_n[node[i].n_n[4]][k];
						Qk_half_np[2][k] = (-1.0 / 8.0)*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] + (3.0 / 4.0)*Qk_iminus_n[node[i].n_n[5]][k] + (3.0 / 8.0)*Qk_iminus_n[i][k];

						/**********************************SMOOTHNESS INDICATOR OR SMOOTHNESS FUNCTION********************************/

						IS_Qim[2][k] = (13.0 / 12.0)*(pow(Qi_iplus[node[node[i].n_n[3]].n_n[3]][k] - 2.0*Qi_iplus[node[i].n_n[3]][k] + Qi_iplus[i][k], 2)) + (1.0 / 4.0)*(pow(Qi_iplus[node[node[i].n_n[3]].n_n[3]][k] - 4.0*Qi_iplus[node[i].n_n[3]][k] + 3.0*Qi_iplus[i][k], 2));
						IS_Qim[1][k] = (13.0 / 12.0)*(pow(Qi_iplus[node[i].n_n[3]][k] - 2.0*Qi_iplus[i][k] + Qi_iplus[node[i].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(Qi_iplus[node[i].n_n[3]][k] - Qi_iplus[node[i].n_n[1]][k], 2));
						IS_Qim[0][k] = (13.0 / 12.0)*(pow(Qi_iplus[i][k] - 2.0*Qi_iplus[node[i].n_n[1]][k] + Qi_iplus[node[node[i].n_n[1]].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus[i][k] - 4.0*Qi_iplus[node[i].n_n[1]][k] + Qi_iplus[node[node[i].n_n[1]].n_n[1]][k], 2));

						IS_Qip[2][k] = (13.0 / 12.0)*(pow(Qi_iplus[node[i].n_n[3]][k] - 2.0*Qi_iplus[i][k] + Qi_iplus[node[i].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(Qi_iplus[node[i].n_n[3]][k] - 4.0*Qi_iplus[i][k] + 3.0*Qi_iplus[node[i].n_n[1]][k], 2));
						IS_Qip[1][k] = (13.0 / 12.0)*(pow(Qi_iplus[i][k] - 2.0*Qi_iplus[node[i].n_n[1]][k] + Qi_iplus[node[node[i].n_n[1]].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(Qi_iplus[i][k] - Qi_iplus[node[node[i].n_n[1]].n_n[1]][k], 2));
						IS_Qip[0][k] = (13.0 / 12.0)*(pow(Qi_iplus[node[i].n_n[1]][k] - 2.0*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k] + Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iplus[node[i].n_n[1]][k] - 4.0*Qi_iplus[node[node[i].n_n[1]].n_n[1]][k] + Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k], 2));

						IS_Qinm[2][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] - 2.0*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] + Qi_iminus_n[node[i].n_n[3]][k], 2)) + (1.0 / 4.0)*(pow(Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] - 4.0*Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] + 3.0*Qi_iminus_n[node[i].n_n[3]][k], 2));
						IS_Qinm[1][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] - 2.0*Qi_iminus_n[node[i].n_n[3]][k] + Qi_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] - Qi_iminus_n[i][k], 2));
						IS_Qinm[0][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k] - 2.0*Qi_iminus_n[i][k] + Qi_iminus_n[node[i].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_n[node[i].n_n[3]][k] - 4.0*Qi_iminus_n[i][k] + Qi_iminus_n[node[i].n_n[1]][k], 2));

						IS_Qinp[2][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] - 2.0*Qi_iminus_n[node[i].n_n[3]][k] + Qi_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] - 4.0*Qi_iminus_n[node[i].n_n[3]][k] + 3.0*Qi_iminus_n[i][k], 2));
						IS_Qinp[1][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k] - 2.0*Qi_iminus_n[i][k] + Qi_iminus_n[node[i].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(Qi_iminus_n[node[i].n_n[3]][k] - Qi_iminus_n[node[i].n_n[1]][k], 2));
						IS_Qinp[0][k] = (13.0 / 12.0)*(pow(Qi_iminus_n[i][k] - 2.0*Qi_iminus_n[node[i].n_n[1]][k] + Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qi_iminus_n[i][k] - 4.0*Qi_iminus_n[node[i].n_n[1]][k] + Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k], 2));

						IS_Qjm[2][k] = (13.0 / 12.0)*(pow(Qj_iplus[node[node[i].n_n[2]].n_n[2]][k] - 2.0*Qj_iplus[node[i].n_n[2]][k] + Qj_iplus[i][k], 2)) + (1.0 / 4.0)*(pow(Qj_iplus[node[node[i].n_n[2]].n_n[2]][k] - 4.0*Qj_iplus[node[i].n_n[2]][k] + 3.0*Qj_iplus[i][k], 2));
						IS_Qjm[1][k] = (13.0 / 12.0)*(pow(Qj_iplus[node[i].n_n[2]][k] - 2.0*Qj_iplus[i][k] + Qj_iplus[node[i].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(Qj_iplus[node[i].n_n[2]][k] - Qj_iplus[node[i].n_n[0]][k], 2));
						IS_Qjm[0][k] = (13.0 / 12.0)*(pow(Qj_iplus[i][k] - 2.0*Qj_iplus[node[i].n_n[0]][k] + Qj_iplus[node[node[i].n_n[0]].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus[i][k] - 4.0*Qj_iplus[node[i].n_n[0]][k] + Qj_iplus[node[node[i].n_n[0]].n_n[0]][k], 2));

						IS_Qjp[2][k] = (13.0 / 12.0)*(pow(Qj_iplus[node[i].n_n[2]][k] - 2.0*Qj_iplus[i][k] + Qj_iplus[node[i].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(Qj_iplus[node[i].n_n[2]][k] - 4.0*Qj_iplus[i][k] + 3.0*Qj_iplus[node[i].n_n[0]][k], 2));
						IS_Qjp[1][k] = (13.0 / 12.0)*(pow(Qj_iplus[i][k] - 2.0*Qj_iplus[node[i].n_n[0]][k] + Qj_iplus[node[node[i].n_n[0]].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(Qj_iplus[i][k] - Qj_iplus[node[node[i].n_n[0]].n_n[0]][k], 2));
						IS_Qjp[0][k] = (13.0 / 12.0)*(pow(Qj_iplus[node[i].n_n[0]][k] - 2.0*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k] + Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iplus[node[i].n_n[0]][k] - 4.0*Qj_iplus[node[node[i].n_n[0]].n_n[0]][k] + Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k], 2));

						IS_Qjnm[2][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] - 2.0*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] + Qj_iminus_n[node[i].n_n[2]][k], 2)) + (1.0 / 4.0)*(pow(Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] - 4.0*Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] + 3.0*Qj_iminus_n[node[i].n_n[2]][k], 2));
						IS_Qjnm[1][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] - 2.0*Qj_iminus_n[node[i].n_n[2]][k] + Qj_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] - Qj_iminus_n[i][k], 2));
						IS_Qjnm[0][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k] - 2.0*Qj_iminus_n[i][k] + Qj_iminus_n[node[i].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_n[node[i].n_n[2]][k] - 4.0*Qj_iminus_n[i][k] + Qj_iminus_n[node[i].n_n[0]][k], 2));

						IS_Qjnp[2][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] - 2.0*Qj_iminus_n[node[i].n_n[2]][k] + Qj_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] - 4.0*Qj_iminus_n[node[i].n_n[2]][k] + 3.0*Qj_iminus_n[i][k], 2));
						IS_Qjnp[1][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k] - 2.0*Qj_iminus_n[i][k] + Qj_iminus_n[node[i].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(Qj_iminus_n[node[i].n_n[2]][k] - Qj_iminus_n[node[i].n_n[0]][k], 2));
						IS_Qjnp[0][k] = (13.0 / 12.0)*(pow(Qj_iminus_n[i][k] - 2.0*Qj_iminus_n[node[i].n_n[0]][k] + Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qj_iminus_n[i][k] - 4.0*Qj_iminus_n[node[i].n_n[0]][k] + Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k], 2));

						IS_Qkm[2][k] = (13.0 / 12.0)*(pow(Qk_iplus[node[node[i].n_n[5]].n_n[5]][k] - 2.0*Qk_iplus[node[i].n_n[5]][k] + Qk_iplus[i][k], 2)) + (1.0 / 4.0)*(pow(Qk_iplus[node[node[i].n_n[5]].n_n[5]][k] - 4.0*Qk_iplus[node[i].n_n[5]][k] + 3.0*Qk_iplus[i][k], 2));
						IS_Qkm[1][k] = (13.0 / 12.0)*(pow(Qk_iplus[node[i].n_n[5]][k] - 2.0*Qk_iplus[i][k] + Qk_iplus[node[i].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(Qk_iplus[node[i].n_n[5]][k] - Qk_iplus[node[i].n_n[4]][k], 2));
						IS_Qkm[0][k] = (13.0 / 12.0)*(pow(Qk_iplus[i][k] - 2.0*Qk_iplus[node[i].n_n[4]][k] + Qk_iplus[node[node[i].n_n[4]].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus[i][k] - 4.0*Qk_iplus[node[i].n_n[4]][k] + Qk_iplus[node[node[i].n_n[4]].n_n[4]][k], 2));

						IS_Qkp[2][k] = (13.0 / 12.0)*(pow(Qk_iplus[node[i].n_n[5]][k] - 2.0*Qk_iplus[i][k] + Qk_iplus[node[i].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(Qk_iplus[node[i].n_n[5]][k] - 4.0*Qk_iplus[i][k] + 3.0*Qk_iplus[node[i].n_n[4]][k], 2));
						IS_Qkp[1][k] = (13.0 / 12.0)*(pow(Qk_iplus[i][k] - 2.0*Qk_iplus[node[i].n_n[4]][k] + Qk_iplus[node[node[i].n_n[4]].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(Qk_iplus[i][k] - Qk_iplus[node[node[i].n_n[4]].n_n[4]][k], 2));
						IS_Qkp[0][k] = (13.0 / 12.0)*(pow(Qk_iplus[node[i].n_n[4]][k] - 2.0*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k] + Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iplus[node[i].n_n[4]][k] - 4.0*Qk_iplus[node[node[i].n_n[4]].n_n[4]][k] + Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k], 2));

						IS_Qknm[2][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] - 2.0*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] + Qk_iminus_n[node[i].n_n[5]][k], 2)) + (1.0 / 4.0)*(pow(Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] - 4.0*Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] + 3.0*Qk_iminus_n[node[i].n_n[5]][k], 2));
						IS_Qknm[1][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] - 2.0*Qk_iminus_n[node[i].n_n[5]][k] + Qk_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] - Qk_iminus_n[i][k], 2));
						IS_Qknm[0][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k] - 2.0*Qk_iminus_n[i][k] + Qk_iminus_n[node[i].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_n[node[i].n_n[5]][k] - 4.0*Qk_iminus_n[i][k] + Qk_iminus_n[node[i].n_n[4]][k], 2));

						IS_Qknp[2][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] - 2.0*Qk_iminus_n[node[i].n_n[5]][k] + Qk_iminus_n[i][k], 2)) + (1.0 / 4.0)*(pow(Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] - 4.0*Qk_iminus_n[node[i].n_n[5]][k] + 3.0*Qk_iminus_n[i][k], 2));
						IS_Qknp[1][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k] - 2.0*Qk_iminus_n[i][k] + Qk_iminus_n[node[i].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(Qk_iminus_n[node[i].n_n[5]][k] - Qk_iminus_n[node[i].n_n[4]][k], 2));
						IS_Qknp[0][k] = (13.0 / 12.0)*(pow(Qk_iminus_n[i][k] - 2.0*Qk_iminus_n[node[i].n_n[4]][k] + Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k], 2)) + (1.0 / 4.0)*(pow(3.0*Qk_iminus_n[i][k] - 4.0*Qk_iminus_n[node[i].n_n[4]][k] + Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k], 2));

						w_Qip[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qip[0][k]), 2.0));
						w_Qip[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qip[1][k]), 2.0));
						w_Qip[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qip[2][k]), 2.0));

						w_Qim[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qim[0][k]), 2.0));
						w_Qim[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qim[1][k]), 2.0));
						w_Qim[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qim[2][k]), 2.0));

						w_Qinp[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinp[0][k]), 2.0));
						w_Qinp[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinp[1][k]), 2.0));
						w_Qinp[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinp[2][k]), 2.0));

						w_Qinm[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qinm[0][k]), 2.0));
						w_Qinm[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qinm[1][k]), 2.0));
						w_Qinm[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qinm[2][k]), 2.0));

						w_Qjp[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjp[0][k]), 2.0));
						w_Qjp[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjp[1][k]), 2.0));
						w_Qjp[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjp[2][k]), 2.0));

						w_Qjm[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjm[0][k]), 2.0));
						w_Qjm[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjm[1][k]), 2.0));
						w_Qjm[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjm[2][k]), 2.0));

						w_Qjnp[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[0][k]), 2.0));
						w_Qjnp[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnp[1][k]), 2.0));
						w_Qjnp[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnp[2][k]), 2.0));

						w_Qjnm[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[0][k]), 2.0));
						w_Qjnm[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qjnm[1][k]), 2.0));
						w_Qjnm[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qjnm[2][k]), 2.0));

						w_Qkp[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkp[0][k]), 2.0));
						w_Qkp[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkp[1][k]), 2.0));
						w_Qkp[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkp[2][k]), 2.0));

						w_Qkm[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qkm[0][k]), 2.0));
						w_Qkm[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qkm[1][k]), 2.0));
						w_Qkm[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qkm[2][k]), 2.0));

						w_Qknp[0][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknp[0][k]), 2.0));
						w_Qknp[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknp[1][k]), 2.0));
						w_Qknp[2][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknp[2][k]), 2.0));

						w_Qknm[0][k] = (5.0 / 16.0) / (pow(Epsilon + (IS_Qknm[0][k]), 2.0));
						w_Qknm[1][k] = (5.0 / 8.0) / (pow(Epsilon + (IS_Qknm[1][k]), 2.0));
						w_Qknm[2][k] = (1.0 / 16.0) / (pow(Epsilon + (IS_Qknm[2][k]), 2.0));

						W_Qip[0][k] = w_Qip[0][k] / (w_Qip[0][k] + w_Qip[1][k] + w_Qip[2][k]);
						W_Qip[1][k] = w_Qip[1][k] / (w_Qip[0][k] + w_Qip[1][k] + w_Qip[2][k]);
						W_Qip[2][k] = w_Qip[2][k] / (w_Qip[0][k] + w_Qip[1][k] + w_Qip[2][k]);

						W_Qinp[0][k] = w_Qinp[0][k] / (w_Qinp[0][k] + w_Qinp[1][k] + w_Qinp[2][k]);
						W_Qinp[1][k] = w_Qinp[1][k] / (w_Qinp[0][k] + w_Qinp[1][k] + w_Qinp[2][k]);
						W_Qinp[2][k] = w_Qinp[2][k] / (w_Qinp[0][k] + w_Qinp[1][k] + w_Qinp[2][k]);

						W_Qim[0][k] = w_Qim[0][k] / (w_Qim[0][k] + w_Qim[1][k] + w_Qim[2][k]);
						W_Qim[1][k] = w_Qim[1][k] / (w_Qim[0][k] + w_Qim[1][k] + w_Qim[2][k]);
						W_Qim[2][k] = w_Qim[2][k] / (w_Qim[0][k] + w_Qim[1][k] + w_Qim[2][k]);

						W_Qinm[0][k] = w_Qinm[0][k] / (w_Qinm[0][k] + w_Qinm[1][k] + w_Qinm[2][k]);
						W_Qinm[1][k] = w_Qinm[1][k] / (w_Qinm[0][k] + w_Qinm[1][k] + w_Qinm[2][k]);
						W_Qinm[2][k] = w_Qinm[2][k] / (w_Qinm[0][k] + w_Qinm[1][k] + w_Qinm[2][k]);

						W_Qjp[0][k] = w_Qjp[0][k] / (w_Qjp[0][k] + w_Qjp[1][k] + w_Qjp[2][k]);
						W_Qjp[1][k] = w_Qjp[1][k] / (w_Qjp[0][k] + w_Qjp[1][k] + w_Qjp[2][k]);
						W_Qjp[2][k] = w_Qjp[2][k] / (w_Qjp[0][k] + w_Qjp[1][k] + w_Qjp[2][k]);

						W_Qjnp[0][k] = w_Qjnp[0][k] / (w_Qjnp[0][k] + w_Qjnp[1][k] + w_Qjnp[2][k]);
						W_Qjnp[1][k] = w_Qjnp[1][k] / (w_Qjnp[0][k] + w_Qjnp[1][k] + w_Qjnp[2][k]);
						W_Qjnp[2][k] = w_Qjnp[2][k] / (w_Qjnp[0][k] + w_Qjnp[1][k] + w_Qjnp[2][k]);

						W_Qjm[0][k] = w_Qjm[0][k] / (w_Qjm[0][k] + w_Qjm[1][k] + w_Qjm[2][k]);
						W_Qjm[1][k] = w_Qjm[1][k] / (w_Qjm[0][k] + w_Qjm[1][k] + w_Qjm[2][k]);
						W_Qjm[2][k] = w_Qjm[2][k] / (w_Qjm[0][k] + w_Qjm[1][k] + w_Qjm[2][k]);

						W_Qjnm[0][k] = w_Qjnm[0][k] / (w_Qjnm[0][k] + w_Qjnm[1][k] + w_Qjnm[2][k]);
						W_Qjnm[1][k] = w_Qjnm[1][k] / (w_Qjnm[0][k] + w_Qjnm[1][k] + w_Qjnm[2][k]);
						W_Qjnm[2][k] = w_Qjnm[2][k] / (w_Qjnm[0][k] + w_Qjnm[1][k] + w_Qjnm[2][k]);

						W_Qkp[0][k] = w_Qkp[0][k] / (w_Qkp[0][k] + w_Qkp[1][k] + w_Qkp[2][k]);
						W_Qkp[1][k] = w_Qkp[1][k] / (w_Qkp[0][k] + w_Qkp[1][k] + w_Qkp[2][k]);
						W_Qkp[2][k] = w_Qkp[2][k] / (w_Qkp[0][k] + w_Qkp[1][k] + w_Qkp[2][k]);

						W_Qknp[0][k] = w_Qknp[0][k] / (w_Qknp[0][k] + w_Qknp[1][k] + w_Qknp[2][k]);
						W_Qknp[1][k] = w_Qknp[1][k] / (w_Qknp[0][k] + w_Qknp[1][k] + w_Qknp[2][k]);
						W_Qknp[2][k] = w_Qknp[2][k] / (w_Qknp[0][k] + w_Qknp[1][k] + w_Qknp[2][k]);

						W_Qkm[0][k] = w_Qkm[0][k] / (w_Qkm[0][k] + w_Qkm[1][k] + w_Qkm[2][k]);
						W_Qkm[1][k] = w_Qkm[1][k] / (w_Qkm[0][k] + w_Qkm[1][k] + w_Qkm[2][k]);
						W_Qkm[2][k] = w_Qkm[2][k] / (w_Qkm[0][k] + w_Qkm[1][k] + w_Qkm[2][k]);

						W_Qknm[0][k] = w_Qknm[0][k] / (w_Qknm[0][k] + w_Qknm[1][k] + w_Qknm[2][k]);
						W_Qknm[1][k] = w_Qknm[1][k] / (w_Qknm[0][k] + w_Qknm[1][k] + w_Qknm[2][k]);
						W_Qknm[2][k] = w_Qknm[2][k] / (w_Qknm[0][k] + w_Qknm[1][k] + w_Qknm[2][k]);

						Qi_iplus_half_pos_char[k] = W_Qip[0][k] * Qi_half_p[0][k] + W_Qip[1][k] * Qi_half_p[1][k] + W_Qip[2][k] * Qi_half_p[2][k];
						Qi_iminus_half_pos_char[k] = W_Qinp[0][k] * Qi_half_np[0][k] + W_Qinp[1][k] * Qi_half_np[1][k] + W_Qinp[2][k] * Qi_half_np[2][k];

						Qi_iplus_half_neg_char[k] = W_Qim[0][k] * Qi_half_m[0][k] + W_Qim[1][k] * Qi_half_m[1][k] + W_Qim[2][k] * Qi_half_m[2][k];
						Qi_iminus_half_neg_char[k] = W_Qinm[0][k] * Qi_halfn_m[0][k] + W_Qinm[1][k] * Qi_halfn_m[1][k] + W_Qinm[2][k] * Qi_halfn_m[2][k];

						Qj_iplus_half_pos_char[k] = W_Qjp[0][k] * Qj_half_p[0][k] + W_Qjp[1][k] * Qj_half_p[1][k] + W_Qjp[2][k] * Qj_half_p[2][k];
						Qj_iminus_half_pos_char[k] = W_Qjnp[0][k] * Qj_half_np[0][k] + W_Qjnp[1][k] * Qj_half_np[1][k] + W_Qjnp[2][k] * Qj_half_np[2][k];

						Qj_iplus_half_neg_char[k] = W_Qjm[0][k] * Qj_half_m[0][k] + W_Qjm[1][k] * Qj_half_m[1][k] + W_Qjm[2][k] * Qj_half_m[2][k];
						Qj_iminus_half_neg_char[k] = W_Qjnm[0][k] * Qj_halfn_m[0][k] + W_Qjnm[1][k] * Qj_halfn_m[1][k] + W_Qjnm[2][k] * Qj_halfn_m[2][k];

						Qk_iplus_half_pos_char[k] = W_Qkp[0][k] * Qk_half_p[0][k] + W_Qkp[1][k] * Qk_half_p[1][k] + W_Qkp[2][k] * Qk_half_p[2][k];
						Qk_iminus_half_pos_char[k] = W_Qknp[0][k] * Qk_half_np[0][k] + W_Qknp[1][k] * Qk_half_np[1][k] + W_Qknp[2][k] * Qk_half_np[2][k];

						Qk_iplus_half_neg_char[k] = W_Qkm[0][k] * Qk_half_m[0][k] + W_Qkm[1][k] * Qk_half_m[1][k] + W_Qkm[2][k] * Qk_half_m[2][k];
						Qk_iminus_half_neg_char[k] = W_Qknm[0][k] * Qk_halfn_m[0][k] + W_Qknm[1][k] * Qk_halfn_m[1][k] + W_Qknm[2][k] * Qk_halfn_m[2][k];


					}

					m = 0;
					for (k = 0; k<5; k++)
					{
						Qi_iplus_half_pos[k] = r_eigen_Qip[m][0] * (Qi_iplus_half_pos_char[0]) + r_eigen_Qip[m][1] * (Qi_iplus_half_pos_char[1]) + r_eigen_Qip[m][2] * (Qi_iplus_half_pos_char[2]) + r_eigen_Qip[m][3] * (Qi_iplus_half_pos_char[3]) + r_eigen_Qip[m][4] * (Qi_iplus_half_pos_char[4]);
						Qi_iplus_half_neg[k] = r_eigen_Qip[m][0] * (Qi_iplus_half_neg_char[0]) + r_eigen_Qip[m][1] * (Qi_iplus_half_neg_char[1]) + r_eigen_Qip[m][2] * (Qi_iplus_half_neg_char[2]) + r_eigen_Qip[m][3] * (Qi_iplus_half_neg_char[3]) + r_eigen_Qip[m][4] * (Qi_iplus_half_neg_char[4]);

						Qi_iminus_half_pos[k] = r_eigen_Qim[m][0] * (Qi_iminus_half_pos_char[0]) + r_eigen_Qim[m][1] * (Qi_iminus_half_pos_char[1]) + r_eigen_Qim[m][2] * (Qi_iminus_half_pos_char[2]) + r_eigen_Qim[m][3] * (Qi_iminus_half_pos_char[3]) + r_eigen_Qim[m][4] * (Qi_iminus_half_pos_char[4]);
						Qi_iminus_half_neg[k] = r_eigen_Qim[m][0] * (Qi_iminus_half_neg_char[0]) + r_eigen_Qim[m][1] * (Qi_iminus_half_neg_char[1]) + r_eigen_Qim[m][2] * (Qi_iminus_half_neg_char[2]) + r_eigen_Qim[m][3] * (Qi_iminus_half_neg_char[3]) + r_eigen_Qim[m][4] * (Qi_iminus_half_neg_char[4]);

						Qj_iplus_half_pos[k] = r_eigen_Qjp[m][0] * (Qj_iplus_half_pos_char[0]) + r_eigen_Qjp[m][1] * (Qj_iplus_half_pos_char[1]) + r_eigen_Qjp[m][2] * (Qj_iplus_half_pos_char[2]) + r_eigen_Qjp[m][3] * (Qj_iplus_half_pos_char[3]) + r_eigen_Qjp[m][4] * (Qj_iplus_half_pos_char[4]);
						Qj_iplus_half_neg[k] = r_eigen_Qjp[m][0] * (Qj_iplus_half_neg_char[0]) + r_eigen_Qjp[m][1] * (Qj_iplus_half_neg_char[1]) + r_eigen_Qjp[m][2] * (Qj_iplus_half_neg_char[2]) + r_eigen_Qjp[m][3] * (Qj_iplus_half_neg_char[3]) + r_eigen_Qjp[m][4] * (Qj_iplus_half_neg_char[4]);

						Qj_iminus_half_pos[k] = r_eigen_Qjm[m][0] * (Qj_iminus_half_pos_char[0]) + r_eigen_Qjm[m][1] * (Qj_iminus_half_pos_char[1]) + r_eigen_Qjm[m][2] * (Qj_iminus_half_pos_char[2]) + r_eigen_Qjm[m][3] * (Qj_iminus_half_pos_char[3]) + r_eigen_Qjm[m][4] * (Qj_iminus_half_pos_char[4]);
						Qj_iminus_half_neg[k] = r_eigen_Qjm[m][0] * (Qj_iminus_half_neg_char[0]) + r_eigen_Qjm[m][1] * (Qj_iminus_half_neg_char[1]) + r_eigen_Qjm[m][2] * (Qj_iminus_half_neg_char[2]) + r_eigen_Qjm[m][3] * (Qj_iminus_half_neg_char[3]) + r_eigen_Qjm[m][4] * (Qj_iminus_half_neg_char[4]);

						Qk_iplus_half_pos[k] = r_eigen_Qkp[m][0] * (Qk_iplus_half_pos_char[0]) + r_eigen_Qkp[m][1] * (Qk_iplus_half_pos_char[1]) + r_eigen_Qkp[m][2] * (Qk_iplus_half_pos_char[2]) + r_eigen_Qkp[m][3] * (Qk_iplus_half_pos_char[3]) + r_eigen_Qkp[m][4] * (Qk_iplus_half_pos_char[4]);
						Qk_iplus_half_neg[k] = r_eigen_Qkp[m][0] * (Qk_iplus_half_neg_char[0]) + r_eigen_Qkp[m][1] * (Qk_iplus_half_neg_char[1]) + r_eigen_Qkp[m][2] * (Qk_iplus_half_neg_char[2]) + r_eigen_Qkp[m][3] * (Qk_iplus_half_neg_char[3]) + r_eigen_Qkp[m][4] * (Qk_iplus_half_neg_char[4]);

						Qk_iminus_half_pos[k] = r_eigen_Qkm[m][0] * (Qk_iminus_half_pos_char[0]) + r_eigen_Qkm[m][1] * (Qk_iminus_half_pos_char[1]) + r_eigen_Qkm[m][2] * (Qk_iminus_half_pos_char[2]) + r_eigen_Qkm[m][3] * (Qk_iminus_half_pos_char[3]) + r_eigen_Qkm[m][4] * (Qk_iminus_half_pos_char[4]);
						Qk_iminus_half_neg[k] = r_eigen_Qkm[m][0] * (Qk_iminus_half_neg_char[0]) + r_eigen_Qkm[m][1] * (Qk_iminus_half_neg_char[1]) + r_eigen_Qkm[m][2] * (Qk_iminus_half_neg_char[2]) + r_eigen_Qkm[m][3] * (Qk_iminus_half_neg_char[3]) + r_eigen_Qkm[m][4] * (Qk_iminus_half_neg_char[4]);

						m++;
					}

					F_ip_pos[0] = Qi_iplus_half_pos[1];
					F_ip_pos[1] = (pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
					F_ip_pos[2] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[2]) / Qi_iplus_half_pos[0];
					F_ip_pos[3] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
					F_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]);

					F_ip_neg[0] = Qi_iplus_half_neg[1];
					F_ip_neg[1] = (pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
					F_ip_neg[2] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[2]) / Qi_iplus_half_neg[0];
					F_ip_neg[3] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
					F_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]);

					F_im_pos[0] = Qi_iminus_half_pos[1];
					F_im_pos[1] = (pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
					F_im_pos[2] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[2]) / Qi_iminus_half_pos[0];
					F_im_pos[3] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
					F_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]);

					F_im_neg[0] = Qi_iminus_half_neg[1];
					F_im_neg[1] = (pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
					F_im_neg[2] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[2]) / Qi_iminus_half_neg[0];
					F_im_neg[3] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
					F_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]);

					E_ip_pos[0] = Qi_iplus_half_pos[2];
					E_ip_pos[1] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[2]) / Qi_iplus_half_pos[0];
					E_ip_pos[2] = (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
					E_ip_pos[3] = (Qi_iplus_half_pos[2] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
					E_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]);

					E_ip_neg[0] = Qi_iplus_half_neg[2];
					E_ip_neg[1] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[2]) / Qi_iplus_half_neg[0];
					E_ip_neg[2] = (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
					E_ip_neg[3] = (Qi_iplus_half_neg[2] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
					E_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]);

					E_im_pos[0] = Qi_iminus_half_pos[2];
					E_im_pos[1] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[2]) / Qi_iminus_half_pos[0];
					E_im_pos[2] = (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
					E_im_pos[3] = (Qi_iminus_half_pos[2] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
					E_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]);

					E_im_neg[0] = Qi_iminus_half_neg[2];
					E_im_neg[1] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[2]) / Qi_iminus_half_neg[0];
					E_im_neg[2] = (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
					E_im_neg[3] = (Qi_iminus_half_neg[2] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
					E_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]);

					G_ip_pos[0] = Qi_iplus_half_pos[3];
					G_ip_pos[1] = (Qi_iplus_half_pos[1] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
					G_ip_pos[2] = (Qi_iplus_half_pos[2] * Qi_iplus_half_pos[3]) / Qi_iplus_half_pos[0];
					G_ip_pos[3] = (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0]) + 0.4*(Qi_iplus_half_pos[4] - 0.5*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])));
					G_ip_pos[4] = (1.4*Qi_iplus_half_pos[4] - 0.2*((pow(Qi_iplus_half_pos[1], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[2], 2) / Qi_iplus_half_pos[0]) + (pow(Qi_iplus_half_pos[3], 2) / Qi_iplus_half_pos[0])))*(Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]);

					G_ip_neg[0] = Qi_iplus_half_neg[3];
					G_ip_neg[1] = (Qi_iplus_half_neg[1] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
					G_ip_neg[2] = (Qi_iplus_half_neg[2] * Qi_iplus_half_neg[3]) / Qi_iplus_half_neg[0];
					G_ip_neg[3] = (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0]) + 0.4*(Qi_iplus_half_neg[4] - 0.5*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])));
					G_ip_neg[4] = (1.4*Qi_iplus_half_neg[4] - 0.2*((pow(Qi_iplus_half_neg[1], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[2], 2) / Qi_iplus_half_neg[0]) + (pow(Qi_iplus_half_neg[3], 2) / Qi_iplus_half_neg[0])))*(Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]);

					G_im_pos[0] = Qi_iminus_half_pos[3];
					G_im_pos[1] = (Qi_iminus_half_pos[1] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
					G_im_pos[2] = (Qi_iminus_half_pos[2] * Qi_iminus_half_pos[3]) / Qi_iminus_half_pos[0];
					G_im_pos[3] = (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0]) + 0.4*(Qi_iminus_half_pos[4] - 0.5*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])));
					G_im_pos[4] = (1.4*Qi_iminus_half_pos[4] - 0.2*((pow(Qi_iminus_half_pos[1], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[2], 2) / Qi_iminus_half_pos[0]) + (pow(Qi_iminus_half_pos[3], 2) / Qi_iminus_half_pos[0])))*(Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]);

					G_im_neg[0] = Qi_iminus_half_neg[3];
					G_im_neg[1] = (Qi_iminus_half_neg[1] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
					G_im_neg[2] = (Qi_iminus_half_neg[2] * Qi_iminus_half_neg[3]) / Qi_iminus_half_neg[0];
					G_im_neg[3] = (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0]) + 0.4*(Qi_iminus_half_neg[4] - 0.5*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])));
					G_im_neg[4] = (1.4*Qi_iminus_half_neg[4] - 0.2*((pow(Qi_iminus_half_neg[1], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[2], 2) / Qi_iminus_half_neg[0]) + (pow(Qi_iminus_half_neg[3], 2) / Qi_iminus_half_neg[0])))*(Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]);


					F_jp_pos[0] = Qj_iplus_half_pos[1];
					F_jp_pos[1] = (pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
					F_jp_pos[2] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[2]) / Qj_iplus_half_pos[0];
					F_jp_pos[3] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
					F_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]);

					F_jp_neg[0] = Qj_iplus_half_neg[1];
					F_jp_neg[1] = (pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
					F_jp_neg[2] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[2]) / Qj_iplus_half_neg[0];
					F_jp_neg[3] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
					F_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]);

					F_jm_pos[0] = Qj_iminus_half_pos[1];
					F_jm_pos[1] = (pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
					F_jm_pos[2] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[2]) / Qj_iminus_half_pos[0];
					F_jm_pos[3] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
					F_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]);

					F_jm_neg[0] = Qj_iminus_half_neg[1];
					F_jm_neg[1] = (pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
					F_jm_neg[2] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[2]) / Qj_iminus_half_neg[0];
					F_jm_neg[3] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
					F_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]);

					E_jp_pos[0] = Qj_iplus_half_pos[2];
					E_jp_pos[1] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[2]) / Qj_iplus_half_pos[0];
					E_jp_pos[2] = (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
					E_jp_pos[3] = (Qj_iplus_half_pos[2] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
					E_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]);

					E_jp_neg[0] = Qj_iplus_half_neg[2];
					E_jp_neg[1] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[2]) / Qj_iplus_half_neg[0];
					E_jp_neg[2] = (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
					E_jp_neg[3] = (Qj_iplus_half_neg[2] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
					E_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]);

					E_jm_pos[0] = Qj_iminus_half_pos[2];
					E_jm_pos[1] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[2]) / Qj_iminus_half_pos[0];
					E_jm_pos[2] = (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
					E_jm_pos[3] = (Qj_iminus_half_pos[2] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
					E_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]);

					E_jm_neg[0] = Qj_iminus_half_neg[2];
					E_jm_neg[1] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[2]) / Qj_iminus_half_neg[0];
					E_jm_neg[2] = (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
					E_jm_neg[3] = (Qj_iminus_half_neg[2] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
					E_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]);

					G_jp_pos[0] = Qj_iplus_half_pos[3];
					G_jp_pos[1] = (Qj_iplus_half_pos[1] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
					G_jp_pos[2] = (Qj_iplus_half_pos[2] * Qj_iplus_half_pos[3]) / Qj_iplus_half_pos[0];
					G_jp_pos[3] = (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0]) + 0.4*(Qj_iplus_half_pos[4] - 0.5*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])));
					G_jp_pos[4] = (1.4*Qj_iplus_half_pos[4] - 0.2*((pow(Qj_iplus_half_pos[1], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[2], 2) / Qj_iplus_half_pos[0]) + (pow(Qj_iplus_half_pos[3], 2) / Qj_iplus_half_pos[0])))*(Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]);

					G_jp_neg[0] = Qj_iplus_half_neg[3];
					G_jp_neg[1] = (Qj_iplus_half_neg[1] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
					G_jp_neg[2] = (Qj_iplus_half_neg[2] * Qj_iplus_half_neg[3]) / Qj_iplus_half_neg[0];
					G_jp_neg[3] = (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0]) + 0.4*(Qj_iplus_half_neg[4] - 0.5*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])));
					G_jp_neg[4] = (1.4*Qj_iplus_half_neg[4] - 0.2*((pow(Qj_iplus_half_neg[1], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[2], 2) / Qj_iplus_half_neg[0]) + (pow(Qj_iplus_half_neg[3], 2) / Qj_iplus_half_neg[0])))*(Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]);

					G_jm_pos[0] = Qj_iminus_half_pos[3];
					G_jm_pos[1] = (Qj_iminus_half_pos[1] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
					G_jm_pos[2] = (Qj_iminus_half_pos[2] * Qj_iminus_half_pos[3]) / Qj_iminus_half_pos[0];
					G_jm_pos[3] = (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0]) + 0.4*(Qj_iminus_half_pos[4] - 0.5*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])));
					G_jm_pos[4] = (1.4*Qj_iminus_half_pos[4] - 0.2*((pow(Qj_iminus_half_pos[1], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[2], 2) / Qj_iminus_half_pos[0]) + (pow(Qj_iminus_half_pos[3], 2) / Qj_iminus_half_pos[0])))*(Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]);

					G_jm_neg[0] = Qj_iminus_half_neg[3];
					G_jm_neg[1] = (Qj_iminus_half_neg[1] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
					G_jm_neg[2] = (Qj_iminus_half_neg[2] * Qj_iminus_half_neg[3]) / Qj_iminus_half_neg[0];
					G_jm_neg[3] = (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0]) + 0.4*(Qj_iminus_half_neg[4] - 0.5*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])));
					G_jm_neg[4] = (1.4*Qj_iminus_half_neg[4] - 0.2*((pow(Qj_iminus_half_neg[1], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[2], 2) / Qj_iminus_half_neg[0]) + (pow(Qj_iminus_half_neg[3], 2) / Qj_iminus_half_neg[0])))*(Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]);


					F_kp_pos[0] = Qk_iplus_half_pos[1];
					F_kp_pos[1] = (pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
					F_kp_pos[2] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[2]) / Qk_iplus_half_pos[0];
					F_kp_pos[3] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
					F_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]);

					F_kp_neg[0] = Qk_iplus_half_neg[1];
					F_kp_neg[1] = (pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
					F_kp_neg[2] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[2]) / Qk_iplus_half_neg[0];
					F_kp_neg[3] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
					F_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]);

					F_km_pos[0] = Qk_iminus_half_pos[1];
					F_km_pos[1] = (pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
					F_km_pos[2] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[2]) / Qk_iminus_half_pos[0];
					F_km_pos[3] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
					F_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]);

					F_km_neg[0] = Qk_iminus_half_neg[1];
					F_km_neg[1] = (pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
					F_km_neg[2] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[2]) / Qk_iminus_half_neg[0];
					F_km_neg[3] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
					F_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]);

					E_kp_pos[0] = Qk_iplus_half_pos[2];
					E_kp_pos[1] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[2]) / Qk_iplus_half_pos[0];
					E_kp_pos[2] = (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
					E_kp_pos[3] = (Qk_iplus_half_pos[2] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
					E_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]);

					E_kp_neg[0] = Qk_iplus_half_neg[2];
					E_kp_neg[1] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[2]) / Qk_iplus_half_neg[0];
					E_kp_neg[2] = (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
					E_kp_neg[3] = (Qk_iplus_half_neg[2] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
					E_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]);

					E_km_pos[0] = Qk_iminus_half_pos[2];
					E_km_pos[1] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[2]) / Qk_iminus_half_pos[0];
					E_km_pos[2] = (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
					E_km_pos[3] = (Qk_iminus_half_pos[2] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
					E_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]);

					E_km_neg[0] = Qk_iminus_half_neg[2];
					E_km_neg[1] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[2]) / Qk_iminus_half_neg[0];
					E_km_neg[2] = (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
					E_km_neg[3] = (Qk_iminus_half_neg[2] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
					E_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]);

					G_kp_pos[0] = Qk_iplus_half_pos[3];
					G_kp_pos[1] = (Qk_iplus_half_pos[1] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
					G_kp_pos[2] = (Qk_iplus_half_pos[2] * Qk_iplus_half_pos[3]) / Qk_iplus_half_pos[0];
					G_kp_pos[3] = (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0]) + 0.4*(Qk_iplus_half_pos[4] - 0.5*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])));
					G_kp_pos[4] = (1.4*Qk_iplus_half_pos[4] - 0.2*((pow(Qk_iplus_half_pos[1], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[2], 2) / Qk_iplus_half_pos[0]) + (pow(Qk_iplus_half_pos[3], 2) / Qk_iplus_half_pos[0])))*(Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]);

					G_kp_neg[0] = Qk_iplus_half_neg[3];
					G_kp_neg[1] = (Qk_iplus_half_neg[1] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
					G_kp_neg[2] = (Qk_iplus_half_neg[2] * Qk_iplus_half_neg[3]) / Qk_iplus_half_neg[0];
					G_kp_neg[3] = (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0]) + 0.4*(Qk_iplus_half_neg[4] - 0.5*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])));
					G_kp_neg[4] = (1.4*Qk_iplus_half_neg[4] - 0.2*((pow(Qk_iplus_half_neg[1], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[2], 2) / Qk_iplus_half_neg[0]) + (pow(Qk_iplus_half_neg[3], 2) / Qk_iplus_half_neg[0])))*(Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]);

					G_km_pos[0] = Qk_iminus_half_pos[3];
					G_km_pos[1] = (Qk_iminus_half_pos[1] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
					G_km_pos[2] = (Qk_iminus_half_pos[2] * Qk_iminus_half_pos[3]) / Qk_iminus_half_pos[0];
					G_km_pos[3] = (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0]) + 0.4*(Qk_iminus_half_pos[4] - 0.5*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])));
					G_km_pos[4] = (1.4*Qk_iminus_half_pos[4] - 0.2*((pow(Qk_iminus_half_pos[1], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[2], 2) / Qk_iminus_half_pos[0]) + (pow(Qk_iminus_half_pos[3], 2) / Qk_iminus_half_pos[0])))*(Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]);

					G_km_neg[0] = Qk_iminus_half_neg[3];
					G_km_neg[1] = (Qk_iminus_half_neg[1] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
					G_km_neg[2] = (Qk_iminus_half_neg[2] * Qk_iminus_half_neg[3]) / Qk_iminus_half_neg[0];
					G_km_neg[3] = (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0]) + 0.4*(Qk_iminus_half_neg[4] - 0.5*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])));
					G_km_neg[4] = (1.4*Qk_iminus_half_neg[4] - 0.2*((pow(Qk_iminus_half_neg[1], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[2], 2) / Qk_iminus_half_neg[0]) + (pow(Qk_iminus_half_neg[3], 2) / Qk_iminus_half_neg[0])))*(Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]);

					zeta_xip = metric[i].zeta_xip;
					zeta_yip = metric[i].zeta_yip;
					zeta_zip = metric[i].zeta_zip;
					eta_xjp = metric[i].eta_xjp;
					eta_yjp = metric[i].eta_yjp;
					eta_zjp = metric[i].eta_zjp;
					xi_xkp = metric[i].xi_xkp;
					xi_ykp = metric[i].xi_ykp;
					xi_zkp = metric[i].xi_zkp;

					zeta_xim = metric[i].zeta_xim;
					zeta_yim = metric[i].zeta_yim;
					zeta_zim = metric[i].zeta_zim;
					eta_xjm = metric[i].eta_xjm;
					eta_yjm = metric[i].eta_yjm;
					eta_zjm = metric[i].eta_zjm;
					xi_xkm = metric[i].xi_xkm;
					xi_ykm = metric[i].xi_ykm;
					xi_zkm = metric[i].xi_zkm;

					eigen_Qip[0] = (zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qi_iplus_half_pos[4] / Qi_iplus_half_pos[0]) - ((pow(Qi_iplus_half_pos[1], 2) + pow(Qi_iplus_half_pos[2], 2) + pow(Qi_iplus_half_pos[3], 2)) / (2.0*pow(Qi_iplus_half_pos[0], 2))))*(zeta_xip*zeta_xip + zeta_yip * zeta_yip + zeta_zip * zeta_zip));
					eigen_Qip[1] = (zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
					eigen_Qip[2] = (zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
					eigen_Qip[3] = (zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0]));
					eigen_Qip[4] = (zeta_xip*(Qi_iplus_half_pos[1] / Qi_iplus_half_pos[0]) + zeta_yip * (Qi_iplus_half_pos[2] / Qi_iplus_half_pos[0]) + zeta_zip * (Qi_iplus_half_pos[3] / Qi_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qi_iplus_half_pos[4] / Qi_iplus_half_pos[0]) - ((pow(Qi_iplus_half_pos[1], 2) + pow(Qi_iplus_half_pos[2], 2) + pow(Qi_iplus_half_pos[3], 2)) / (2.0*pow(Qi_iplus_half_pos[0], 2))))*(zeta_xip*zeta_xip + zeta_yip * zeta_yip + zeta_zip * zeta_zip));

					eigen_Qim[0] = (zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qi_iplus_half_neg[4] / Qi_iplus_half_neg[0]) - ((pow(Qi_iplus_half_neg[1], 2) + pow(Qi_iplus_half_neg[2], 2) + pow(Qi_iplus_half_neg[3], 2)) / (2.0*pow(Qi_iplus_half_neg[0], 2))))*(zeta_xip*zeta_xip + zeta_yip * zeta_yip + zeta_zip * zeta_zip));
					eigen_Qim[1] = (zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
					eigen_Qim[2] = (zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
					eigen_Qim[3] = (zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0]));
					eigen_Qim[4] = (zeta_xip*(Qi_iplus_half_neg[1] / Qi_iplus_half_neg[0]) + zeta_yip * (Qi_iplus_half_neg[2] / Qi_iplus_half_neg[0]) + zeta_zip * (Qi_iplus_half_neg[3] / Qi_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qi_iplus_half_neg[4] / Qi_iplus_half_neg[0]) - ((pow(Qi_iplus_half_neg[1], 2) + pow(Qi_iplus_half_neg[2], 2) + pow(Qi_iplus_half_neg[3], 2)) / (2.0*pow(Qi_iplus_half_neg[0], 2))))*(zeta_xip*zeta_xip + zeta_yip * zeta_yip + zeta_zip * zeta_zip));

					eigen_Qinp[0] = (zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qi_iminus_half_pos[4] / Qi_iminus_half_pos[0]) - ((pow(Qi_iminus_half_pos[1], 2) + pow(Qi_iminus_half_pos[2], 2) + pow(Qi_iminus_half_pos[3], 2)) / (2.0*pow(Qi_iminus_half_pos[0], 2))))*(zeta_xim*zeta_xim + zeta_yim * zeta_yim + zeta_zim * zeta_zim));
					eigen_Qinp[1] = (zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
					eigen_Qinp[2] = (zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
					eigen_Qinp[3] = (zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0]));
					eigen_Qinp[4] = (zeta_xim*(Qi_iminus_half_pos[1] / Qi_iminus_half_pos[0]) + zeta_yim * (Qi_iminus_half_pos[2] / Qi_iminus_half_pos[0]) + zeta_zim * (Qi_iminus_half_pos[3] / Qi_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qi_iminus_half_pos[4] / Qi_iminus_half_pos[0]) - ((pow(Qi_iminus_half_pos[1], 2) + pow(Qi_iminus_half_pos[2], 2) + pow(Qi_iminus_half_pos[3], 2)) / (2.0*pow(Qi_iminus_half_pos[0], 2))))*(zeta_xim*zeta_xim + zeta_yim * zeta_yim + zeta_zim * zeta_zim));

					eigen_Qinm[0] = (zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qi_iminus_half_neg[4] / Qi_iminus_half_neg[0]) - ((pow(Qi_iminus_half_neg[1], 2) + pow(Qi_iminus_half_neg[2], 2) + pow(Qi_iminus_half_neg[3], 2)) / (2.0*pow(Qi_iminus_half_neg[0], 2))))*(zeta_xim*zeta_xim + zeta_yim * zeta_yim + zeta_zim * zeta_zim));
					eigen_Qinm[1] = (zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
					eigen_Qinm[2] = (zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
					eigen_Qinm[3] = (zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0]));
					eigen_Qinm[4] = (zeta_xim*(Qi_iminus_half_neg[1] / Qi_iminus_half_neg[0]) + zeta_yim * (Qi_iminus_half_neg[2] / Qi_iminus_half_neg[0]) + zeta_zim * (Qi_iminus_half_neg[3] / Qi_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qi_iminus_half_neg[4] / Qi_iminus_half_neg[0]) - ((pow(Qi_iminus_half_neg[1], 2) + pow(Qi_iminus_half_neg[2], 2) + pow(Qi_iminus_half_neg[3], 2)) / (2.0*pow(Qi_iminus_half_neg[0], 2))))*(zeta_xim*zeta_xim + zeta_yim * zeta_yim + zeta_zim * zeta_zim));

					eigen_Qjp[0] = (eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qj_iplus_half_pos[4] / Qj_iplus_half_pos[0]) - ((pow(Qj_iplus_half_pos[1], 2) + pow(Qj_iplus_half_pos[2], 2) + pow(Qj_iplus_half_pos[3], 2)) / (2.0*pow(Qj_iplus_half_pos[0], 2))))*(eta_xjp*eta_xjp + eta_yjp * eta_yjp + eta_zjp * eta_zjp));
					eigen_Qjp[1] = (eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
					eigen_Qjp[2] = (eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
					eigen_Qjp[3] = (eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0]));
					eigen_Qjp[4] = (eta_xjp*(Qj_iplus_half_pos[1] / Qj_iplus_half_pos[0]) + eta_yjp * (Qj_iplus_half_pos[2] / Qj_iplus_half_pos[0]) + eta_zjp * (Qj_iplus_half_pos[3] / Qj_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qj_iplus_half_pos[4] / Qj_iplus_half_pos[0]) - ((pow(Qj_iplus_half_pos[1], 2) + pow(Qj_iplus_half_pos[2], 2) + pow(Qj_iplus_half_pos[3], 2)) / (2.0*pow(Qj_iplus_half_pos[0], 2))))*(eta_xjp*eta_xjp + eta_yjp * eta_yjp + eta_zjp * eta_zjp));

					eigen_Qjm[0] = (eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qj_iplus_half_neg[4] / Qj_iplus_half_neg[0]) - ((pow(Qj_iplus_half_neg[1], 2) + pow(Qj_iplus_half_neg[2], 2) + pow(Qj_iplus_half_neg[3], 2)) / (2.0*pow(Qj_iplus_half_neg[0], 2))))*(eta_xjp*eta_xjp + eta_yjp * eta_yjp + eta_zjp * eta_zjp));
					eigen_Qjm[1] = (eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
					eigen_Qjm[2] = (eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
					eigen_Qjm[3] = (eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0]));
					eigen_Qjm[4] = (eta_xjp*(Qj_iplus_half_neg[1] / Qj_iplus_half_neg[0]) + eta_yjp * (Qj_iplus_half_neg[2] / Qj_iplus_half_neg[0]) + eta_zjp * (Qj_iplus_half_neg[3] / Qj_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qj_iplus_half_neg[4] / Qj_iplus_half_neg[0]) - ((pow(Qj_iplus_half_neg[1], 2) + pow(Qj_iplus_half_neg[2], 2) + pow(Qj_iplus_half_neg[3], 2)) / (2.0*pow(Qj_iplus_half_neg[0], 2))))*(eta_xjp*eta_xjp + eta_yjp * eta_yjp + eta_zjp * eta_zjp));

					eigen_Qjnp[0] = (eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qj_iminus_half_pos[4] / Qj_iminus_half_pos[0]) - ((pow(Qj_iminus_half_pos[1], 2) + pow(Qj_iminus_half_pos[2], 2) + pow(Qj_iminus_half_pos[3], 2)) / (2.0*pow(Qj_iminus_half_pos[0], 2))))*(eta_xjm*eta_xjm + eta_yjm * eta_yjm + eta_zjm * eta_zjm));
					eigen_Qjnp[1] = (eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
					eigen_Qjnp[2] = (eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
					eigen_Qjnp[3] = (eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0]));
					eigen_Qjnp[4] = (eta_xjm*(Qj_iminus_half_pos[1] / Qj_iminus_half_pos[0]) + eta_yjm * (Qj_iminus_half_pos[2] / Qj_iminus_half_pos[0]) + eta_zjm * (Qj_iminus_half_pos[3] / Qj_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qj_iminus_half_pos[4] / Qj_iminus_half_pos[0]) - ((pow(Qj_iminus_half_pos[1], 2) + pow(Qj_iminus_half_pos[2], 2) + pow(Qj_iminus_half_pos[3], 2)) / (2.0*pow(Qj_iminus_half_pos[0], 2))))*(eta_xjm*eta_xjm + eta_yjm * eta_yjm + eta_zjm * eta_zjm));

					eigen_Qjnm[0] = (eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qj_iminus_half_neg[4] / Qj_iminus_half_neg[0]) - ((pow(Qj_iminus_half_neg[1], 2) + pow(Qj_iminus_half_neg[2], 2) + pow(Qj_iminus_half_neg[3], 2)) / (2.0*pow(Qj_iminus_half_neg[0], 2))))*(eta_xjm*eta_xjm + eta_yjm * eta_yjm + eta_zjm * eta_zjm));
					eigen_Qjnm[1] = (eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
					eigen_Qjnm[2] = (eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
					eigen_Qjnm[3] = (eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0]));
					eigen_Qjnm[4] = (eta_xjm*(Qj_iminus_half_neg[1] / Qj_iminus_half_neg[0]) + eta_yjm * (Qj_iminus_half_neg[2] / Qj_iminus_half_neg[0]) + eta_zjm * (Qj_iminus_half_neg[3] / Qj_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qj_iminus_half_neg[4] / Qj_iminus_half_neg[0]) - ((pow(Qj_iminus_half_neg[1], 2) + pow(Qj_iminus_half_neg[2], 2) + pow(Qj_iminus_half_neg[3], 2)) / (2.0*pow(Qj_iminus_half_neg[0], 2))))*(eta_xjm*eta_xjm + eta_yjm * eta_yjm + eta_zjm * eta_zjm));

					eigen_Qkp[0] = (xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0])) - sqrt(1.4*0.4*((Qk_iplus_half_pos[4] / Qk_iplus_half_pos[0]) - ((pow(Qk_iplus_half_pos[1], 2) + pow(Qk_iplus_half_pos[2], 2) + pow(Qk_iplus_half_pos[3], 2)) / (2.0*pow(Qk_iplus_half_pos[0], 2))))*(xi_xkp*xi_xkp + xi_ykp * xi_ykp + xi_zkp * xi_zkp));
					eigen_Qkp[1] = (xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
					eigen_Qkp[2] = (xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
					eigen_Qkp[3] = (xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0]));
					eigen_Qkp[4] = (xi_xkp*(Qk_iplus_half_pos[1] / Qk_iplus_half_pos[0]) + xi_ykp * (Qk_iplus_half_pos[2] / Qk_iplus_half_pos[0]) + xi_zkp * (Qk_iplus_half_pos[3] / Qk_iplus_half_pos[0])) + sqrt(1.4*0.4*((Qk_iplus_half_pos[4] / Qk_iplus_half_pos[0]) - ((pow(Qk_iplus_half_pos[1], 2) + pow(Qk_iplus_half_pos[2], 2) + pow(Qk_iplus_half_pos[3], 2)) / (2.0*pow(Qk_iplus_half_pos[0], 2))))*(xi_xkp*xi_xkp + xi_ykp * xi_ykp + xi_zkp * xi_zkp));

					eigen_Qkm[0] = (xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0])) - sqrt(1.4*0.4*((Qk_iplus_half_neg[4] / Qk_iplus_half_neg[0]) - ((pow(Qk_iplus_half_neg[1], 2) + pow(Qk_iplus_half_neg[2], 2) + pow(Qk_iplus_half_neg[3], 2)) / (2.0*pow(Qk_iplus_half_neg[0], 2))))*(xi_xkp*xi_xkp + xi_ykp * xi_ykp + xi_zkp * xi_zkp));
					eigen_Qkm[1] = (xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
					eigen_Qkm[2] = (xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
					eigen_Qkm[3] = (xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0]));
					eigen_Qkm[4] = (xi_xkp*(Qk_iplus_half_neg[1] / Qk_iplus_half_neg[0]) + xi_ykp * (Qk_iplus_half_neg[2] / Qk_iplus_half_neg[0]) + xi_zkp * (Qk_iplus_half_neg[3] / Qk_iplus_half_neg[0])) + sqrt(1.4*0.4*((Qk_iplus_half_neg[4] / Qk_iplus_half_neg[0]) - ((pow(Qk_iplus_half_neg[1], 2) + pow(Qk_iplus_half_neg[2], 2) + pow(Qk_iplus_half_neg[3], 2)) / (2.0*pow(Qk_iplus_half_neg[0], 2))))*(xi_xkp*xi_xkp + xi_ykp * xi_ykp + xi_zkp * xi_zkp));

					eigen_Qknp[0] = (xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0])) - sqrt(1.4*0.4*((Qk_iminus_half_pos[4] / Qk_iminus_half_pos[0]) - ((pow(Qk_iminus_half_pos[1], 2) + pow(Qk_iminus_half_pos[2], 2) + pow(Qk_iminus_half_pos[3], 2)) / (2.0*pow(Qk_iminus_half_pos[0], 2))))*(xi_xkm*xi_xkm + xi_ykm * xi_ykm + xi_zkm * xi_zkm));
					eigen_Qknp[1] = (xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
					eigen_Qknp[2] = (xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
					eigen_Qknp[3] = (xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0]));
					eigen_Qknp[4] = (xi_xkm*(Qk_iminus_half_pos[1] / Qk_iminus_half_pos[0]) + xi_ykm * (Qk_iminus_half_pos[2] / Qk_iminus_half_pos[0]) + xi_zkm * (Qk_iminus_half_pos[3] / Qk_iminus_half_pos[0])) + sqrt(1.4*0.4*((Qk_iminus_half_pos[4] / Qk_iminus_half_pos[0]) - ((pow(Qk_iminus_half_pos[1], 2) + pow(Qk_iminus_half_pos[2], 2) + pow(Qk_iminus_half_pos[3], 2)) / (2.0*pow(Qk_iminus_half_pos[0], 2))))*(xi_xkm*xi_xkm + xi_ykm * xi_ykm + xi_zkm * xi_zkm));

					eigen_Qknm[0] = (xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0])) - sqrt(1.4*0.4*((Qk_iminus_half_neg[4] / Qk_iminus_half_neg[0]) - ((pow(Qk_iminus_half_neg[1], 2) + pow(Qk_iminus_half_neg[2], 2) + pow(Qk_iminus_half_neg[3], 2)) / (2.0*pow(Qk_iminus_half_neg[0], 2))))*(xi_xkm*xi_xkm + xi_ykm * xi_ykm + xi_zkm * xi_zkm));
					eigen_Qknm[1] = (xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
					eigen_Qknm[2] = (xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
					eigen_Qknm[3] = (xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0]));
					eigen_Qknm[4] = (xi_xkm*(Qk_iminus_half_neg[1] / Qk_iminus_half_neg[0]) + xi_ykm * (Qk_iminus_half_neg[2] / Qk_iminus_half_neg[0]) + xi_zkm * (Qk_iminus_half_neg[3] / Qk_iminus_half_neg[0])) + sqrt(1.4*0.4*((Qk_iminus_half_neg[4] / Qk_iminus_half_neg[0]) - ((pow(Qk_iminus_half_neg[1], 2) + pow(Qk_iminus_half_neg[2], 2) + pow(Qk_iminus_half_neg[3], 2)) / (2.0*pow(Qk_iminus_half_neg[0], 2))))*(xi_xkm*xi_xkm + xi_ykm * xi_ykm + xi_zkm * xi_zkm));

					for (k = 0; k<5; k++)
					{
						F_ip_pos_comp[k] = deter[i].ip*(zeta_xip*F_ip_pos[k] + zeta_yip * E_ip_pos[k] + zeta_zip * G_ip_pos[k]);
						F_ip_neg_comp[k] = deter[i].ip*(zeta_xip*F_ip_neg[k] + zeta_yip * E_ip_neg[k] + zeta_zip * G_ip_neg[k]);

						F_im_pos_comp[k] = deter[i].im*(zeta_xim*F_im_pos[k] + zeta_yim * E_im_pos[k] + zeta_zim * G_im_pos[k]);
						F_im_neg_comp[k] = deter[i].im*(zeta_xim*F_im_neg[k] + zeta_yim * E_im_neg[k] + zeta_zim * G_im_neg[k]);

						E_ip_pos_comp[k] = deter[i].jp*(eta_xjp*F_jp_pos[k] + eta_yjp * E_jp_pos[k] + eta_zjp * G_jp_pos[k]);
						E_ip_neg_comp[k] = deter[i].jp*(eta_xjp*F_jp_neg[k] + eta_yjp * E_jp_neg[k] + eta_zjp * G_jp_neg[k]);

						E_im_pos_comp[k] = deter[i].jm*(eta_xjm*F_jm_pos[k] + eta_yjm * E_jm_pos[k] + eta_zjm * G_jm_pos[k]);
						E_im_neg_comp[k] = deter[i].jm*(eta_xjm*F_jm_neg[k] + eta_yjm * E_jm_neg[k] + eta_zjm * G_jm_neg[k]);

						G_ip_pos_comp[k] = deter[i].kp*(xi_xkp*F_kp_pos[k] + xi_ykp * E_kp_pos[k] + xi_zkp * G_kp_pos[k]);
						G_ip_neg_comp[k] = deter[i].kp*(xi_xkp*F_kp_neg[k] + xi_ykp * E_kp_neg[k] + xi_zkp * G_kp_neg[k]);

						G_im_pos_comp[k] = deter[i].km*(xi_xkm*F_km_pos[k] + xi_ykm * E_km_pos[k] + xi_zkm * G_km_pos[k]);
						G_im_neg_comp[k] = deter[i].km*(xi_xkm*F_km_neg[k] + xi_ykm * E_km_neg[k] + xi_zkm * G_km_neg[k]);

						Qi_iplus_half_pos[k] = deter[i].ip*Qi_iplus_half_pos[k];
						Qi_iplus_half_neg[k] = deter[i].ip*Qi_iplus_half_neg[k];

						Qi_iminus_half_pos[k] = deter[i].im*Qi_iminus_half_pos[k];
						Qi_iminus_half_neg[k] = deter[i].im*Qi_iminus_half_neg[k];

						Qj_iplus_half_pos[k] = deter[i].jp*Qj_iplus_half_pos[k];
						Qj_iplus_half_neg[k] = deter[i].jp*Qj_iplus_half_neg[k];

						Qj_iminus_half_pos[k] = deter[i].jm*Qj_iminus_half_pos[k];
						Qj_iminus_half_neg[k] = deter[i].jm*Qj_iminus_half_neg[k];

						Qk_iplus_half_pos[k] = deter[i].kp*Qk_iplus_half_pos[k];
						Qk_iplus_half_neg[k] = deter[i].kp*Qk_iplus_half_neg[k];

						Qk_iminus_half_pos[k] = deter[i].km*Qk_iminus_half_pos[k];
						Qk_iminus_half_neg[k] = deter[i].km*Qk_iminus_half_neg[k];

						/******************************************************************************************************************************************************/
						/******************************************************CENTRAL 4th ORDER VISCOUS TERMS*****************************************************************/

						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[1]].loc == 100))
						{
							dFv[k] = (1.0 / 6.0)*(11.0*Fv[i][k] - 18.0*Fv[node[i].n_n[3]][k] + 9.0*Fv[node[node[i].n_n[3]].n_n[3]][k] - 2.0*Fv[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]);
							dF[k] = (1.0 / 6.0)*(11.0*F[i][k] - 18.0*F[node[i].n_n[3]][k] + 9.0*F[node[node[i].n_n[3]].n_n[3]][k] - 2.0*F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k]);
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[3]].loc == 100))
						{
							dFv[k] = (1.0 / 6.0)*(-11.0*Fv[i][k] + 18.0*Fv[node[i].n_n[1]][k] - 9.0*Fv[node[node[i].n_n[1]].n_n[1]][k] + 2.0*Fv[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
							dF[k] = (1.0 / 6.0)*(-11.0*F[i][k] + 18.0*F[node[i].n_n[1]][k] - 9.0*F[node[node[i].n_n[1]].n_n[1]][k] + 2.0*F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[1]].loc > 0 && node[node[i].n_n[1]].loc <= 52 && (node[node[node[i].n_n[1]].n_n[1]].loc == 100))
						{
							dFv[k] = (1.0 / 6.0)*((2.0)*Fv[node[i].n_n[1]][k] + 3.0*Fv[i][k] - 6.0*Fv[node[i].n_n[3]][k] + Fv[node[node[i].n_n[3]].n_n[3]][k]);
							dF[k] = (1.0 / 6.0)*((2.0)*F[node[i].n_n[1]][k] + 3.0*F[i][k] - 6.0*F[node[i].n_n[3]][k] + F[node[node[i].n_n[3]].n_n[3]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[3]].loc > 0 && node[node[i].n_n[3]].loc <= 52 && (node[node[node[i].n_n[3]].n_n[3]].loc == 100))
						{
							dFv[k] = (1.0 / 6.0)*((-2.0)*Fv[node[i].n_n[3]][k] - 3.0*Fv[i][k] + 6.0*Fv[node[i].n_n[1]][k] - Fv[node[node[i].n_n[1]].n_n[1]][k]);
							dF[k] = (1.0 / 6.0)*((-2.0)*F[node[i].n_n[3]][k] - 3.0*F[i][k] + 6.0*F[node[i].n_n[1]][k] - F[node[node[i].n_n[1]].n_n[1]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[3]].n_n[3]].loc != 100 && node[node[i].n_n[3]].loc != 100 && node[node[node[i].n_n[1]].n_n[1]].loc != 100 && node[node[i].n_n[1]].loc != 100)
						{
							dFv[k] = (1.0 / 12.0)*(8.0*(Fv[node[i].n_n[1]][k] - Fv[node[i].n_n[3]][k]) - (Fv[node[node[i].n_n[1]].n_n[1]][k] - Fv[node[node[i].n_n[3]].n_n[3]][k]));
							dF[k] = (1.0 / 12.0)*(8.0*(F[node[i].n_n[1]][k] - F[node[i].n_n[3]][k]) - (F[node[node[i].n_n[1]].n_n[1]][k] - F[node[node[i].n_n[3]].n_n[3]][k]));
						}


						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[0]].loc == 100))
						{
							dEv[k] = (1.0 / 6.0)*(11.0*Ev[i][k] - 18.0*Ev[node[i].n_n[2]][k] + 9.0*Ev[node[node[i].n_n[2]].n_n[2]][k] - 2.0*Ev[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]);
							dE[k] = (1.0 / 6.0)*(11.0*E[i][k] - 18.0*E[node[i].n_n[2]][k] + 9.0*E[node[node[i].n_n[2]].n_n[2]][k] - 2.0*E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k]);
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[2]].loc == 100))
						{
							dEv[k] = (1.0 / 6.0)*(-11.0*Ev[i][k] + 18.0*Ev[node[i].n_n[0]][k] - 9.0*Ev[node[node[i].n_n[0]].n_n[0]][k] + 2.0*Ev[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
							dE[k] = (1.0 / 6.0)*(-11.0*E[i][k] + 18.0*E[node[i].n_n[0]][k] - 9.0*E[node[node[i].n_n[0]].n_n[0]][k] + 2.0*E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[0]].loc > 0 && node[node[i].n_n[0]].loc <= 52 && (node[node[node[i].n_n[0]].n_n[0]].loc == 100))
						{
							dEv[k] = (1.0 / 6.0)*(2.0*Ev[node[i].n_n[0]][k] + 3.0*Ev[i][k] - 6.0*Ev[node[i].n_n[2]][k] + Ev[node[node[i].n_n[2]].n_n[2]][k]);
							dE[k] = (1.0 / 6.0)*(2.0*E[node[i].n_n[0]][k] + 3.0*E[i][k] - 6.0*E[node[i].n_n[2]][k] + E[node[node[i].n_n[2]].n_n[2]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[2]].loc > 0 && node[node[i].n_n[2]].loc <= 52 && (node[node[node[i].n_n[2]].n_n[2]].loc == 100))
						{
							dEv[k] = (1.0 / 6.0)*(-2.0*Ev[node[i].n_n[2]][k] - 3.0*Ev[i][k] + 6.0*Ev[node[i].n_n[0]][k] - Ev[node[node[i].n_n[0]].n_n[0]][k]);
							dE[k] = (1.0 / 6.0)*(-2.0*E[node[i].n_n[2]][k] - 3.0*E[i][k] + 6.0*E[node[i].n_n[0]][k] - E[node[node[i].n_n[0]].n_n[0]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[2]].n_n[2]].loc != 100 && node[node[i].n_n[2]].loc != 100 && node[node[node[i].n_n[0]].n_n[0]].loc != 100 && node[node[i].n_n[0]].loc != 100)
						{
							dEv[k] = (1.0 / 12.0)*(8.0*(Ev[node[i].n_n[0]][k] - Ev[node[i].n_n[2]][k]) - (Ev[node[node[i].n_n[0]].n_n[0]][k] - Ev[node[node[i].n_n[2]].n_n[2]][k]));
							dE[k] = (1.0 / 12.0)*(8.0*(E[node[i].n_n[0]][k] - E[node[i].n_n[2]][k]) - (E[node[node[i].n_n[0]].n_n[0]][k] - E[node[node[i].n_n[2]].n_n[2]][k]));
						}


						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[4]].loc == 100))
						{
							dGv[k] = (1.0 / 6.0)*(11.0*Gv[i][k] - 18.0*Gv[node[i].n_n[5]][k] + 9.0*Gv[node[node[i].n_n[5]].n_n[5]][k] - 2.0*Gv[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]);
							dG[k] = (1.0 / 6.0)*(11.0*G[i][k] - 18.0*G[node[i].n_n[5]][k] + 9.0*G[node[node[i].n_n[5]].n_n[5]][k] - 2.0*G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k]);
						}
						if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[5]].loc == 100))
						{
							dGv[k] = (1.0 / 6.0)*(-11.0*Gv[i][k] + 18.0*Gv[node[i].n_n[4]][k] - 9.0*Gv[node[node[i].n_n[4]].n_n[4]][k] + 2.0*Gv[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
							dG[k] = (1.0 / 6.0)*(-11.0*G[i][k] + 18.0*G[node[i].n_n[4]][k] - 9.0*G[node[node[i].n_n[4]].n_n[4]][k] + 2.0*G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[4]].loc > 0 && node[node[i].n_n[4]].loc <= 52 && (node[node[node[i].n_n[4]].n_n[4]].loc == 100))
						{
							dGv[k] = (1.0 / 6.0)*(2.0*Gv[node[i].n_n[4]][k] + 3.0*Gv[i][k] - 6.0*Gv[node[i].n_n[5]][k] + Gv[node[node[i].n_n[5]].n_n[5]][k]);
							dG[k] = (1.0 / 6.0)*(2.0*G[node[i].n_n[4]][k] + 3.0*G[i][k] - 6.0*G[node[i].n_n[5]][k] + G[node[node[i].n_n[5]].n_n[5]][k]);
						}
						if (node[i].loc == 0 && node[node[i].n_n[5]].loc > 0 && node[node[i].n_n[5]].loc <= 52 && (node[node[node[i].n_n[5]].n_n[5]].loc == 100))
						{
							dGv[k] = (1.0 / 6.0)*(-2.0*Gv[node[i].n_n[5]][k] - 3.0*Gv[i][k] + 6.0*Gv[node[i].n_n[4]][k] - Gv[node[node[i].n_n[4]].n_n[4]][k]);
							dG[k] = (1.0 / 6.0)*(-2.0*G[node[i].n_n[5]][k] - 3.0*G[i][k] + 6.0*G[node[i].n_n[4]][k] - G[node[node[i].n_n[4]].n_n[4]][k]);
						}
						if (node[i].loc == 0 && node[node[node[i].n_n[5]].n_n[5]].loc != 100 && node[node[i].n_n[5]].loc != 100 && node[node[node[i].n_n[4]].n_n[4]].loc != 100 && node[node[i].n_n[4]].loc != 100)
						{
							dGv[k] = (1.0 / 12.0)*(8.0*(Gv[node[i].n_n[4]][k] - Gv[node[i].n_n[5]][k]) - (Gv[node[node[i].n_n[4]].n_n[4]][k] - Gv[node[node[i].n_n[5]].n_n[5]][k]));
							dG[k] = (1.0 / 12.0)*(8.0*(G[node[i].n_n[4]][k] - G[node[i].n_n[5]][k]) - (G[node[node[i].n_n[4]].n_n[4]][k] - G[node[node[i].n_n[5]].n_n[5]][k]));
						}
						

						if (fabs(eigen_Qip[4]) >= fabs(eigen_Qim[4]))
						{
							alpha_u_ip[k] = fabs(eigen_Qip[4]);
						}
						else if (fabs(eigen_Qim[4]) >= fabs(eigen_Qip[4]))
						{
							alpha_u_ip[k] = fabs(eigen_Qim[4]);
						}

						if (fabs(eigen_Qjp[4]) >= fabs(eigen_Qjm[4]))
						{
							alpha_v_jp[k] = fabs(eigen_Qjp[4]);
						}
						else if (fabs(eigen_Qjm[4]) >= fabs(eigen_Qjp[4]))
						{
							alpha_v_jp[k] = fabs(eigen_Qjm[4]);
						}

						if (fabs(eigen_Qkp[4]) >= fabs(eigen_Qkm[4]))
						{
							alpha_w_kp[k] = fabs(eigen_Qkp[4]);
						}
						else if (fabs(eigen_Qkm[4]) >= fabs(eigen_Qkp[4]))
						{
							alpha_w_kp[k] = fabs(eigen_Qkm[4]);
						}

						if (fabs(eigen_Qinp[4]) >= fabs(eigen_Qinm[4]))
						{
							alpha_u_im[k] = fabs(eigen_Qinp[4]);
						}
						else if (fabs(eigen_Qinm[4]) >= fabs(eigen_Qinp[4]))
						{
							alpha_u_im[k] = fabs(eigen_Qinm[4]);
						}

						if (fabs(eigen_Qjnp[4]) >= fabs(eigen_Qjnm[4]))
						{
							alpha_v_jm[k] = fabs(eigen_Qjnp[4]);
						}
						else if (fabs(eigen_Qjnm[4]) >= fabs(eigen_Qjnp[4]))
						{
							alpha_v_jm[k] = fabs(eigen_Qjnm[4]);
						}

						if (fabs(eigen_Qknp[4]) >= fabs(eigen_Qknm[4]))
						{
							alpha_w_km[k] = fabs(eigen_Qknp[4]);
						}
						else if (fabs(eigen_Qknm[4]) >= fabs(eigen_Qknp[4]))
						{
							alpha_w_km[k] = fabs(eigen_Qknm[4]);
						}

						F_ip[k] = 0.5*(F_ip_pos_comp[k] + F_ip_neg_comp[k] - alpha_u_ip[k] * (Qi_iplus_half_pos[k] - Qi_iplus_half_neg[k]));
						F_im[k] = 0.5*(F_im_pos_comp[k] + F_im_neg_comp[k] - alpha_u_im[k] * (Qi_iminus_half_pos[k] - Qi_iminus_half_neg[k]));

						E_jp[k] = 0.5*(E_ip_pos_comp[k] + E_ip_neg_comp[k] - alpha_v_jp[k] * (Qj_iplus_half_pos[k] - Qj_iplus_half_neg[k]));
						E_jm[k] = 0.5*(E_im_pos_comp[k] + E_im_neg_comp[k] - alpha_v_jm[k] * (Qj_iminus_half_pos[k] - Qj_iminus_half_neg[k]));

						G_kp[k] = 0.5*(G_ip_pos_comp[k] + G_ip_neg_comp[k] - alpha_w_kp[k] * (Qk_iplus_half_pos[k] - Qk_iplus_half_neg[k]));
						G_km[k] = 0.5*(G_im_pos_comp[k] + G_im_neg_comp[k] - alpha_w_km[k] * (Qk_iminus_half_pos[k] - Qk_iminus_half_neg[k]));

						h_F_ip[k] = F_ip[k];
						h_F_im[k] = F_im[k];

						h_E_ip[k] = E_jp[k];
						h_E_im[k] = E_jm[k];

						h_G_ip[k] = G_kp[k];
						h_G_im[k] = G_km[k];

						d2F_d2z_ip[k] = (1.0 / 48.0)*(-5.0*F[node[node[i].n_n[3]].n_n[3]][k] + 39.0*F[node[i].n_n[3]][k] - 34.0*F[i][k] - 34.0*F[node[i].n_n[1]][k] + 39.0*F[node[node[i].n_n[1]].n_n[1]][k] - 5.0*F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
						d2F_d2z_im[k] = (1.0 / 48.0)*(-5.0*F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] + 39.0*F[node[node[i].n_n[3]].n_n[3]][k] - 34.0*F[node[i].n_n[3]][k] - 34.0*F[i][k] + 39.0*F[node[i].n_n[1]][k] - 5.0*F[node[node[i].n_n[1]].n_n[1]][k]);

						d2E_d2e_ip[k] = (1.0 / 48.0)*(-5.0*E[node[node[i].n_n[2]].n_n[2]][k] + 39.0*E[node[i].n_n[2]][k] - 34.0*E[i][k] - 34.0*E[node[i].n_n[0]][k] + 39.0*E[node[node[i].n_n[0]].n_n[0]][k] - 5.0*E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
						d2E_d2e_im[k] = (1.0 / 48.0)*(-5.0*E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] + 39.0*E[node[node[i].n_n[2]].n_n[2]][k] - 34.0*E[node[i].n_n[2]][k] - 34.0*E[i][k] + 39.0*E[node[i].n_n[0]][k] - 5.0*E[node[node[i].n_n[0]].n_n[0]][k]);

						d2G_d2x_ip[k] = (1.0 / 48.0)*(-5.0*G[node[node[i].n_n[5]].n_n[5]][k] + 39.0*G[node[i].n_n[5]][k] - 34.0*G[i][k] - 34.0*G[node[i].n_n[4]][k] + 39.0*G[node[node[i].n_n[4]].n_n[4]][k] - 5.0*G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
						d2G_d2x_im[k] = (1.0 / 48.0)*(-5.0*G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] + 39.0*G[node[node[i].n_n[5]].n_n[5]][k] - 34.0*G[node[i].n_n[5]][k] - 34.0*G[i][k] + 39.0*G[node[i].n_n[4]][k] - 5.0*G[node[node[i].n_n[4]].n_n[4]][k]);

						d4F_d4z_ip[k] = (1.0 / 2.0)*(F[node[node[i].n_n[3]].n_n[3]][k] - 3.0*F[node[i].n_n[3]][k] + 2.0*F[i][k] + 2.0*F[node[i].n_n[1]][k] - 3.0*F[node[node[i].n_n[1]].n_n[1]][k] + F[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k]);
						d4F_d4z_im[k] = (1.0 / 2.0)*(F[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] - 3.0*F[node[node[i].n_n[3]].n_n[3]][k] + 2.0*F[node[i].n_n[3]][k] + 2.0*F[i][k] - 3.0*F[node[i].n_n[1]][k] + F[node[node[i].n_n[1]].n_n[1]][k]);

						d4E_d4e_ip[k] = (1.0 / 2.0)*(E[node[node[i].n_n[2]].n_n[2]][k] - 3.0*E[node[i].n_n[2]][k] + 2.0*E[i][k] + 2.0*E[node[i].n_n[0]][k] - 3.0*E[node[node[i].n_n[0]].n_n[0]][k] + E[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k]);
						d4E_d4e_im[k] = (1.0 / 2.0)*(E[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] - 3.0*E[node[node[i].n_n[2]].n_n[2]][k] + 2.0*E[node[i].n_n[2]][k] + 2.0*E[i][k] - 3.0*E[node[i].n_n[0]][k] + E[node[node[i].n_n[0]].n_n[0]][k]);

						d4G_d4x_ip[k] = (1.0 / 2.0)*(G[node[node[i].n_n[5]].n_n[5]][k] - 3.0*G[node[i].n_n[5]][k] + 2.0*G[i][k] + 2.0*G[node[i].n_n[4]][k] - 3.0*G[node[node[i].n_n[4]].n_n[4]][k] + G[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k]);
						d4G_d4x_im[k] = (1.0 / 2.0)*(G[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] - 3.0*G[node[node[i].n_n[5]].n_n[5]][k] + 2.0*G[node[i].n_n[5]][k] + 2.0*G[i][k] - 3.0*G[node[i].n_n[4]][k] + G[node[node[i].n_n[4]].n_n[4]][k]);

						F_ip_h[k] = h_F_ip[k] - (1.0 / 24.0)*d2F_d2z_ip[k] + (7.0 / 5760.0)*d4F_d4z_ip[k];
						F_im_h[k] = h_F_im[k] - (1.0 / 24.0)*d2F_d2z_im[k] + (7.0 / 5760.0)*d4F_d4z_im[k];

						E_ip_h[k] = h_E_ip[k] - (1.0 / 24.0)*d2E_d2e_ip[k] + (7.0 / 5760.0)*d4E_d4e_ip[k];
						E_im_h[k] = h_E_im[k] - (1.0 / 24.0)*d2E_d2e_im[k] + (7.0 / 5760.0)*d4E_d4e_im[k];

						G_ip_h[k] = h_G_ip[k] - (1.0 / 24.0)*d2G_d2x_ip[k] + (7.0 / 5760.0)*d4G_d4x_ip[k];
						G_im_h[k] = h_G_im[k] - (1.0 / 24.0)*d2G_d2x_im[k] + (7.0 / 5760.0)*d4G_d4x_im[k];

						dF_W[k] = F_ip_h[k] - F_im_h[k];
						dE_W[k] = E_ip_h[k] - E_im_h[k];
						dG_W[k] = G_ip_h[k] - G_im_h[k];

						dF_C[k] = dF[k];
						dE_C[k] = dE[k];
						dG_C[k] = dG[k];

						dF[k] = DUCROS[i] * dF_W[k] + (1.0 - DUCROS[i])*dF_C[k];
						dE[k] = DUCROS[i] * dE_W[k] + (1.0 - DUCROS[i])*dE_C[k];
						dG[k] = DUCROS[i] * dG_W[k] + (1.0 - DUCROS[i])*dG_C[k];

						dF[k] = F_ip_h[k] - F_im_h[k];
						dE[k] = E_ip_h[k] - E_im_h[k];
						dG[k] = G_ip_h[k] - G_im_h[k];

						L[i][j][k] = (-1.0)*((1.0 / del_zeta)*dF[k] + (1.0 / del_eta)*dE[k] + (1.0 / del_eta)*dG[k] + (1.0 / del_zeta)*dFv[k] + (1.0 / del_eta)*dEv[k] + (1.0 / del_eta)*dGv[k]);
						/***********************3rd order RK Scheme***************************************************************************/
						if (j == 0)
						{
							U[i][gar][k] = U[i][0][k] + del_t * L[i][0][k];
						}
						else if (j == 1)
						{
							U[i][gar][k] = (3.0 / 4.0)*U[i][0][k] + (1.0 / 4.0)*U[i][1][k] + 0.25*del_t*L[i][1][k];
						}
						else if (j == 2)
						{
							U[i][gar][k] = (U[i][0][k] / 3.0) + ((2.0 / 3.0)*U[i][2][k]) + (2.0 / 3.0)*(del_t*L[i][2][k]);
							lm = 0;
						}
						final_U[k] = (1.0 / det[i])*U[i][gar][k];
					}
					
					rho[lm][i] = final_U[0];
					u[lm][i] = final_U[1] / final_U[0];
					v[lm][i] = final_U[2] / final_U[0];
					w[lm][i] = final_U[3] / final_U[0];
					e[lm][i] = (final_U[4] / final_U[0]) - 0.5*(pow(u[lm][i], 2.0) + pow(v[lm][i], 2.0) + pow(w[lm][i], 2.0));
					p[lm][i] = e[lm][i] * rho[lm][i] * 0.4;
					t[lm][i] = (1.4*Mach*Mach)*(p[lm][i] / rho[lm][i]);
					a[lm][i] = sqrt(1.4*p[lm][i] / rho[lm][i]);

					if (j == 2)
					{
						diver[i][1].u = u[0][i];
						diver[i][1].v = v[0][i];
						diver[i][1].e = e[0][i];

						if (max_div_um < fabs(diver[i][1].u - diver[i][0].u))
						{
							max_div_um = fabs(diver[i][1].u - diver[i][0].u);
							u_loc = i;
						}
						if (max_div_vm < fabs(diver[i][1].v - diver[i][0].v))
						{
							max_div_vm = fabs(diver[i][1].v - diver[i][0].v);
							v_loc = i;
						}
						if (max_div_em < fabs(diver[i][1].e - diver[i][0].e))
						{
							max_div_em = fabs(diver[i][1].e - diver[i][0].e);
							e_loc = i;
						}

						diver[i][0].u = diver[i][1].u;
						diver[i][0].v = diver[i][1].v;
						diver[i][0].e = diver[i][1].e;
					}
				}
			}

			position = 0;
			for (i = 0; i<neigh_pro; i++)
			{
				MPI_Pack_size(recv_c[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize1);
				buffer2 = new char[memsize1];                                          /***********carefull with buffer1 and buffer2******************/
				MPI_Pack_size(proc_node[c[i]] * 9, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
				buffer = new char[memsize];
				position = 0;
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&u[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&v[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&w[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&p[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&t[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&rho[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&e[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}
				for (k = 0; k<recv_c[c[i]]; k++)
				{
					MPI_Pack(&mu[lm][loc_dat[i][k]], 1, MPI_DOUBLE, buffer2, memsize1, &position, MPI_COMM_WORLD);
				}

				MPI_Sendrecv(buffer2, memsize1, MPI_PACKED, c[i], c[i], buffer, memsize, MPI_PACKED, c[i], myrank, MPI_COMM_WORLD, &status);
				position = 0;

				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &u[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &v[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &w[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &p[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &t[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &rho[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &e[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				for (k = 0; k<proc_node[c[i]]; k++)
				{
					MPI_Unpack(buffer, memsize, &position, &mu[lm][recv_b[i][k]], 1, MPI_DOUBLE, MPI_COMM_WORLD);
				}
				delete[] buffer;
				delete[] buffer2;
				k = 0;
			}

			initial = 1;
			intialise(lm);

			for (k = 1; k <= no_of_tip_send; k++)
			{
				u[lm][tip[k].node] = 0.0;
				v[lm][tip[k].node] = 0.0;
				w[lm][tip[k].node] = 0.0;
				p[lm][tip[k].node] = p[lm][node[tip[k].node].n_n[3]];
				t[lm][tip[k].node] = t[lm][node[tip[k].node].n_n[3]];
				rho[lm][tip[k].node] = (1.4*Mach*Mach)*(p[lm][tip[k].node] / t[lm][tip[k].node]);
				e[lm][tip[k].node] = p[lm][tip[k].node] / (0.4*rho[lm][tip[k].node]);
			}

			for (k = 1; k <= no_of_tip_recv; k++)
			{
				u[lm][tip_recv[k].node] = 0.0;
				v[lm][tip_recv[k].node] = 0.0;
				w[lm][tip_recv[k].node] = 0.0;
				p[lm][tip_recv[k].node] = p[lm][node[tip_recv[k].node].n_n[3]];
				t[lm][tip_recv[k].node] = t[lm][node[tip_recv[k].node].n_n[3]];
				rho[lm][tip_recv[k].node] = (1.4*Mach*Mach)*(p[lm][tip_recv[k].node] / t[lm][tip_recv[k].node]);
				e[lm][tip_recv[k].node] = p[lm][tip_recv[k].node] / (0.4*rho[lm][tip_recv[k].node]);
			}
		}
		
		position = 0;
		MPI_Pack_size(3, MPI_DOUBLE, MPI_COMM_WORLD, &memsize);
		buffer = new char[memsize];
		MPI_Pack(&max_div_um, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_vm, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Pack(&max_div_em, 1, MPI_DOUBLE, buffer, memsize, &position, MPI_COMM_WORLD);
		MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
		delete[] buffer;

		if ((iter % 100) == 0)
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

			MPI_Send(buffer, position, MPI_PACKED, 0, myrank, MPI_COMM_WORLD);
			delete[] buffer;
		}
	}
}
