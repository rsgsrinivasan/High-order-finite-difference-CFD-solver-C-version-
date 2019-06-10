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

void weno_solver_0()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_jp[i] - 0.5*(roe_u_jp[i] * roe_u_jp[i] + roe_v_jp[i] * roe_v_jp[i] + roe_w_jp[i] * roe_w_jp[i])));

	//	if (hk == 0)
	{
		zeta_xi = metric[i].eta_xjp;
		zeta_yi = metric[i].eta_yjp;
		zeta_zi = metric[i].eta_zjp;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_jp[i] + ky_bar * roe_v_jp[i] + kz_bar * roe_w_jp[i];
	phi_sq = 0.5*0.4*(roe_u_jp[i] * roe_u_jp[i] + roe_v_jp[i] * roe_v_jp[i] + roe_w_jp[i] * roe_w_jp[i]);
	alpha = roe_rho_jp[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_jp[i] * roe_a[i]);

	r_eigen_Qjp[0][0] = kx_bar;
	r_eigen_Qjp[0][1] = ky_bar;
	r_eigen_Qjp[0][2] = kz_bar;
	r_eigen_Qjp[0][3] = alpha;
	r_eigen_Qjp[0][4] = alpha;

	r_eigen_Qjp[1][0] = kx_bar * roe_u_jp[i];
	r_eigen_Qjp[1][1] = ky_bar * roe_u_jp[i] - kz_bar * roe_rho_jp[i];
	r_eigen_Qjp[1][2] = kz_bar * roe_u_jp[i] + ky_bar * roe_rho_jp[i];
	r_eigen_Qjp[1][3] = alpha * (roe_u_jp[i] + kx_bar * roe_a[i]);
	r_eigen_Qjp[1][4] = alpha * (roe_u_jp[i] - kx_bar * roe_a[i]);

	r_eigen_Qjp[2][0] = kx_bar * roe_v_jp[i] + kz_bar * roe_rho_jp[i];
	r_eigen_Qjp[2][1] = ky_bar * roe_v_jp[i];
	r_eigen_Qjp[2][2] = kz_bar * roe_v_jp[i] - kx_bar * roe_rho_jp[i];
	r_eigen_Qjp[2][3] = alpha * (roe_v_jp[i] + ky_bar * roe_a[i]);
	r_eigen_Qjp[2][4] = alpha * (roe_v_jp[i] - ky_bar * roe_a[i]);

	r_eigen_Qjp[3][0] = kx_bar * roe_w_jp[i] - ky_bar * roe_rho_jp[i];
	r_eigen_Qjp[3][1] = ky_bar * roe_w_jp[i] + kx_bar * roe_rho_jp[i];
	r_eigen_Qjp[3][2] = kz_bar * roe_w_jp[i];
	r_eigen_Qjp[3][3] = alpha * (roe_w_jp[i] + kz_bar * roe_a[i]);
	r_eigen_Qjp[3][4] = alpha * (roe_w_jp[i] - kz_bar * roe_a[i]);

	r_eigen_Qjp[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_jp[i] * (kz_bar*roe_v_jp[i] - ky_bar * roe_w_jp[i]);
	r_eigen_Qjp[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_jp[i] * (kx_bar*roe_w_jp[i] - kz_bar * roe_u_jp[i]);
	r_eigen_Qjp[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_jp[i] * (ky_bar*roe_u_jp[i] - kx_bar * roe_v_jp[i]);
	r_eigen_Qjp[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qjp[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qjp[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_jp[i] - ky_bar * roe_w_jp[i]) / roe_rho_jp[i]);
	l_eigen_Qjp[0][1] = kx_bar * 0.4*roe_u_jp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[0][2] = (kx_bar*0.4*roe_v_jp[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_jp[i]);
	l_eigen_Qjp[0][3] = (kx_bar*0.4*roe_w_jp[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_jp[i]);
	l_eigen_Qjp[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_jp[i] - kz_bar * roe_u_jp[i]) / roe_rho_jp[i]);
	l_eigen_Qjp[1][1] = (ky_bar*0.4*roe_u_jp[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_jp[i];
	l_eigen_Qjp[1][2] = ky_bar * 0.4*roe_v_jp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[1][3] = (ky_bar*0.4*roe_w_jp[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_jp[i]);
	l_eigen_Qjp[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_jp[i] - kx_bar * roe_v_jp[i]) / roe_rho_jp[i]);
	l_eigen_Qjp[2][1] = kz_bar * 0.4*roe_u_jp[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_jp[i]);
	l_eigen_Qjp[2][2] = kz_bar * 0.4*roe_v_jp[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_jp[i]);
	l_eigen_Qjp[2][3] = kz_bar * 0.4*roe_w_jp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjp[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjp[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qjp[3][1] = -beta * (0.4*roe_u_jp[i] - kx_bar * roe_a[i]);
	l_eigen_Qjp[3][2] = -beta * (0.4*roe_v_jp[i] - ky_bar * roe_a[i]);
	l_eigen_Qjp[3][3] = -beta * (0.4*roe_w_jp[i] - kz_bar * roe_a[i]);
	l_eigen_Qjp[3][4] = beta * 0.4;

	l_eigen_Qjp[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qjp[4][1] = -beta * (0.4*roe_u_jp[i] + kx_bar * roe_a[i]);
	l_eigen_Qjp[4][2] = -beta * (0.4*roe_v_jp[i] + ky_bar * roe_a[i]);
	l_eigen_Qjp[4][3] = -beta * (0.4*roe_w_jp[i] + kz_bar * roe_a[i]);
	l_eigen_Qjp[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 0)
		{
			Qj_iplus[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;
			Qj_iplus[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[node[i].n_n[0]].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[i].n_n[0]].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[i].n_n[0]].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[i].n_n[0]].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[i].n_n[0]].n_n[0]][4].ip;
			Qj_iplus[node[i].n_n[0]][k] = l_eigen_Qjp[m][0] * Q[node[i].n_n[0]][0].ip + l_eigen_Qjp[m][1] * Q[node[i].n_n[0]][1].ip + l_eigen_Qjp[m][2] * Q[node[i].n_n[0]][2].ip + l_eigen_Qjp[m][3] * Q[node[i].n_n[0]][3].ip + l_eigen_Qjp[m][4] * Q[node[i].n_n[0]][4].ip;
			Qj_iplus[i][k] = l_eigen_Qjp[m][0] * Q[i][0].ip + l_eigen_Qjp[m][1] * Q[i][1].ip + l_eigen_Qjp[m][2] * Q[i][2].ip + l_eigen_Qjp[m][3] * Q[i][3].ip + l_eigen_Qjp[m][4] * Q[i][4].ip;
			Qj_iplus[node[i].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[i].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[i].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[i].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[i].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[i].n_n[2]][4].ip;
			Qj_iplus[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[node[i].n_n[2]].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[i].n_n[2]].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[i].n_n[2]].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[i].n_n[2]].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[i].n_n[2]].n_n[2]][4].ip;
			Qj_iplus[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_Qjp[m][0] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip + l_eigen_Qjp[m][1] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip + l_eigen_Qjp[m][2] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip + l_eigen_Qjp[m][3] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip + l_eigen_Qjp[m][4] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;
		}

		m++;
	}
}


void weno_solver_1()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_ip[i] - 0.5*(roe_u_ip[i] * roe_u_ip[i] + roe_v_ip[i] * roe_v_ip[i] + roe_w_ip[i] * roe_w_ip[i])));

	zeta_xi = metric[i].zeta_xip;
	zeta_yi = metric[i].zeta_yip;
	zeta_zi = metric[i].zeta_zip;

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_ip[i] + ky_bar * roe_v_ip[i] + kz_bar * roe_w_ip[i];
	phi_sq = 0.5*0.4*(roe_u_ip[i] * roe_u_ip[i] + roe_v_ip[i] * roe_v_ip[i] + roe_w_ip[i] * roe_w_ip[i]);
	alpha = roe_rho_ip[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_ip[i] * roe_a[i]);

	r_eigen_Qip[0][0] = kx_bar;
	r_eigen_Qip[0][1] = ky_bar;
	r_eigen_Qip[0][2] = kz_bar;
	r_eigen_Qip[0][3] = alpha;
	r_eigen_Qip[0][4] = alpha;

	r_eigen_Qip[1][0] = kx_bar * roe_u_ip[i];
	r_eigen_Qip[1][1] = ky_bar * roe_u_ip[i] - kz_bar * roe_rho_ip[i];
	r_eigen_Qip[1][2] = kz_bar * roe_u_ip[i] + ky_bar * roe_rho_ip[i];
	r_eigen_Qip[1][3] = alpha * (roe_u_ip[i] + kx_bar * roe_a[i]);
	r_eigen_Qip[1][4] = alpha * (roe_u_ip[i] - kx_bar * roe_a[i]);

	r_eigen_Qip[2][0] = kx_bar * roe_v_ip[i] + kz_bar * roe_rho_ip[i];
	r_eigen_Qip[2][1] = ky_bar * roe_v_ip[i];
	r_eigen_Qip[2][2] = kz_bar * roe_v_ip[i] - kx_bar * roe_rho_ip[i];
	r_eigen_Qip[2][3] = alpha * (roe_v_ip[i] + ky_bar * roe_a[i]);
	r_eigen_Qip[2][4] = alpha * (roe_v_ip[i] - ky_bar * roe_a[i]);

	r_eigen_Qip[3][0] = kx_bar * roe_w_ip[i] - ky_bar * roe_rho_ip[i];
	r_eigen_Qip[3][1] = ky_bar * roe_w_ip[i] + kx_bar * roe_rho_ip[i];
	r_eigen_Qip[3][2] = kz_bar * roe_w_ip[i];
	r_eigen_Qip[3][3] = alpha * (roe_w_ip[i] + kz_bar * roe_a[i]);
	r_eigen_Qip[3][4] = alpha * (roe_w_ip[i] - kz_bar * roe_a[i]);

	r_eigen_Qip[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_ip[i] * (kz_bar*roe_v_ip[i] - ky_bar * roe_w_ip[i]);
	r_eigen_Qip[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_ip[i] * (kx_bar*roe_w_ip[i] - kz_bar * roe_u_ip[i]);
	r_eigen_Qip[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_ip[i] * (ky_bar*roe_u_ip[i] - kx_bar * roe_v_ip[i]);
	r_eigen_Qip[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qip[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qip[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_ip[i] - ky_bar * roe_w_ip[i]) / roe_rho_ip[i]);
	l_eigen_Qip[0][1] = kx_bar * 0.4*roe_u_ip[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[0][2] = (kx_bar*0.4*roe_v_ip[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_ip[i]);
	l_eigen_Qip[0][3] = (kx_bar*0.4*roe_w_ip[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_ip[i]);
	l_eigen_Qip[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_ip[i] - kz_bar * roe_u_ip[i]) / roe_rho_ip[i]);
	l_eigen_Qip[1][1] = (ky_bar*0.4*roe_u_ip[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_ip[i];
	l_eigen_Qip[1][2] = ky_bar * 0.4*roe_v_ip[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[1][3] = (ky_bar*0.4*roe_w_ip[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_ip[i]);
	l_eigen_Qip[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_ip[i] - kx_bar * roe_v_ip[i]) / roe_rho_ip[i]);
	l_eigen_Qip[2][1] = kz_bar * 0.4*roe_u_ip[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_ip[i]);
	l_eigen_Qip[2][2] = kz_bar * 0.4*roe_v_ip[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_ip[i]);
	l_eigen_Qip[2][3] = kz_bar * 0.4*roe_w_ip[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qip[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qip[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qip[3][1] = -beta * (0.4*roe_u_ip[i] - kx_bar * roe_a[i]);
	l_eigen_Qip[3][2] = -beta * (0.4*roe_v_ip[i] - ky_bar * roe_a[i]);
	l_eigen_Qip[3][3] = -beta * (0.4*roe_w_ip[i] - kz_bar * roe_a[i]);
	l_eigen_Qip[3][4] = beta * 0.4;

	l_eigen_Qip[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qip[4][1] = -beta * (0.4*roe_u_ip[i] + kx_bar * roe_a[i]);
	l_eigen_Qip[4][2] = -beta * (0.4*roe_v_ip[i] + ky_bar * roe_a[i]);
	l_eigen_Qip[4][3] = -beta * (0.4*roe_w_ip[i] + kz_bar * roe_a[i]);
	l_eigen_Qip[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		Qi_iplus[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;
		Qi_iplus[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[node[i].n_n[1]].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[node[i].n_n[1]].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[node[i].n_n[1]].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[node[i].n_n[1]].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[node[i].n_n[1]].n_n[1]][4].ip;
		Qi_iplus[node[i].n_n[1]][k] = l_eigen_Qip[m][0] * Q[node[i].n_n[1]][0].ip + l_eigen_Qip[m][1] * Q[node[i].n_n[1]][1].ip + l_eigen_Qip[m][2] * Q[node[i].n_n[1]][2].ip + l_eigen_Qip[m][3] * Q[node[i].n_n[1]][3].ip + l_eigen_Qip[m][4] * Q[node[i].n_n[1]][4].ip;
		Qi_iplus[i][k] = l_eigen_Qip[m][0] * Q[i][0].ip + l_eigen_Qip[m][1] * Q[i][1].ip + l_eigen_Qip[m][2] * Q[i][2].ip + l_eigen_Qip[m][3] * Q[i][3].ip + l_eigen_Qip[m][4] * Q[i][4].ip;
		Qi_iplus[node[i].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[i].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[i].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[i].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[i].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[i].n_n[3]][4].ip;
		Qi_iplus[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[node[i].n_n[3]].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[node[i].n_n[3]].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[node[i].n_n[3]].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[node[i].n_n[3]].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[node[i].n_n[3]].n_n[3]][4].ip;
		Qi_iplus[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_Qip[m][0] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip + l_eigen_Qip[m][1] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip + l_eigen_Qip[m][2] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip + l_eigen_Qip[m][3] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip + l_eigen_Qip[m][4] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;


		m++;
	}
}

void weno_solver_2()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_jm[i] - 0.5*(roe_u_jm[i] * roe_u_jm[i] + roe_v_jm[i] * roe_v_jm[i] + roe_w_jm[i] * roe_w_jm[i])));

	zeta_xi = metric[i].eta_xjm;
	zeta_yi = metric[i].eta_yjm;
	zeta_zi = metric[i].eta_zjm;


	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_jm[i] + ky_bar * roe_v_jm[i] + kz_bar * roe_w_jm[i];
	phi_sq = 0.5*0.4*(roe_u_jm[i] * roe_u_jm[i] + roe_v_jm[i] * roe_v_jm[i] + roe_w_jm[i] * roe_w_jm[i]);
	alpha = roe_rho_jm[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_jm[i] * roe_a[i]);

	r_eigen_Qjm[0][0] = kx_bar;
	r_eigen_Qjm[0][1] = ky_bar;
	r_eigen_Qjm[0][2] = kz_bar;
	r_eigen_Qjm[0][3] = alpha;
	r_eigen_Qjm[0][4] = alpha;

	r_eigen_Qjm[1][0] = kx_bar * roe_u_jm[i];
	r_eigen_Qjm[1][1] = ky_bar * roe_u_jm[i] - kz_bar * roe_rho_jm[i];
	r_eigen_Qjm[1][2] = kz_bar * roe_u_jm[i] + ky_bar * roe_rho_jm[i];
	r_eigen_Qjm[1][3] = alpha * (roe_u_jm[i] + kx_bar * roe_a[i]);
	r_eigen_Qjm[1][4] = alpha * (roe_u_jm[i] - kx_bar * roe_a[i]);

	r_eigen_Qjm[2][0] = kx_bar * roe_v_jm[i] + kz_bar * roe_rho_jm[i];
	r_eigen_Qjm[2][1] = ky_bar * roe_v_jm[i];
	r_eigen_Qjm[2][2] = kz_bar * roe_v_jm[i] - kx_bar * roe_rho_jm[i];
	r_eigen_Qjm[2][3] = alpha * (roe_v_jm[i] + ky_bar * roe_a[i]);
	r_eigen_Qjm[2][4] = alpha * (roe_v_jm[i] - ky_bar * roe_a[i]);

	r_eigen_Qjm[3][0] = kx_bar * roe_w_jm[i] - ky_bar * roe_rho_jm[i];
	r_eigen_Qjm[3][1] = ky_bar * roe_w_jm[i] + kx_bar * roe_rho_jm[i];
	r_eigen_Qjm[3][2] = kz_bar * roe_w_jm[i];
	r_eigen_Qjm[3][3] = alpha * (roe_w_jm[i] + kz_bar * roe_a[i]);
	r_eigen_Qjm[3][4] = alpha * (roe_w_jm[i] - kz_bar * roe_a[i]);

	r_eigen_Qjm[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_jm[i] * (kz_bar*roe_v_jm[i] - ky_bar * roe_w_jm[i]);
	r_eigen_Qjm[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_jm[i] * (kx_bar*roe_w_jm[i] - kz_bar * roe_u_jm[i]);
	r_eigen_Qjm[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_jm[i] * (ky_bar*roe_u_jm[i] - kx_bar * roe_v_jm[i]);
	r_eigen_Qjm[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qjm[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qjm[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_jm[i] - ky_bar * roe_w_jm[i]) / roe_rho_jm[i]);
	l_eigen_Qjm[0][1] = kx_bar * 0.4*roe_u_jm[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[0][2] = (kx_bar*0.4*roe_v_jm[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_jm[i]);
	l_eigen_Qjm[0][3] = (kx_bar*0.4*roe_w_jm[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_jm[i]);
	l_eigen_Qjm[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_jm[i] - kz_bar * roe_u_jm[i]) / roe_rho_jm[i]);
	l_eigen_Qjm[1][1] = (ky_bar*0.4*roe_u_jm[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_jm[i];
	l_eigen_Qjm[1][2] = ky_bar * 0.4*roe_v_jm[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[1][3] = (ky_bar*0.4*roe_w_jm[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_jm[i]);
	l_eigen_Qjm[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_jm[i] - kx_bar * roe_v_jm[i]) / roe_rho_jm[i]);
	l_eigen_Qjm[2][1] = kz_bar * 0.4*roe_u_jm[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_jm[i]);
	l_eigen_Qjm[2][2] = kz_bar * 0.4*roe_v_jm[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_jm[i]);
	l_eigen_Qjm[2][3] = kz_bar * 0.4*roe_w_jm[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qjm[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qjm[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qjm[3][1] = -beta * (0.4*roe_u_jm[i] - kx_bar * roe_a[i]);
	l_eigen_Qjm[3][2] = -beta * (0.4*roe_v_jm[i] - ky_bar * roe_a[i]);
	l_eigen_Qjm[3][3] = -beta * (0.4*roe_w_jm[i] - kz_bar * roe_a[i]);
	l_eigen_Qjm[3][4] = beta * 0.4;

	l_eigen_Qjm[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qjm[4][1] = -beta * (0.4*roe_u_jm[i] + kx_bar * roe_a[i]);
	l_eigen_Qjm[4][2] = -beta * (0.4*roe_v_jm[i] + ky_bar * roe_a[i]);
	l_eigen_Qjm[4][3] = -beta * (0.4*roe_w_jm[i] + kz_bar * roe_a[i]);
	l_eigen_Qjm[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 2)
		{
			Qj_iminus_n[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][4].ip;
			Qj_iminus_n[node[node[i].n_n[0]].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[node[i].n_n[0]].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[i].n_n[0]].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[i].n_n[0]].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[i].n_n[0]].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[i].n_n[0]].n_n[0]][4].ip;
			Qj_iminus_n[node[i].n_n[0]][k] = l_eigen_Qjm[m][0] * Q[node[i].n_n[0]][0].ip + l_eigen_Qjm[m][1] * Q[node[i].n_n[0]][1].ip + l_eigen_Qjm[m][2] * Q[node[i].n_n[0]][2].ip + l_eigen_Qjm[m][3] * Q[node[i].n_n[0]][3].ip + l_eigen_Qjm[m][4] * Q[node[i].n_n[0]][4].ip;
			Qj_iminus_n[i][k] = l_eigen_Qjm[m][0] * Q[i][0].ip + l_eigen_Qjm[m][1] * Q[i][1].ip + l_eigen_Qjm[m][2] * Q[i][2].ip + l_eigen_Qjm[m][3] * Q[i][3].ip + l_eigen_Qjm[m][4] * Q[i][4].ip;
			Qj_iminus_n[node[i].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[i].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[i].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[i].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[i].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[i].n_n[2]][4].ip;
			Qj_iminus_n[node[node[i].n_n[2]].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[node[i].n_n[2]].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[i].n_n[2]].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[i].n_n[2]].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[i].n_n[2]].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[i].n_n[2]].n_n[2]][4].ip;
			Qj_iminus_n[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][k] = l_eigen_Qjm[m][0] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][0].ip + l_eigen_Qjm[m][1] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][1].ip + l_eigen_Qjm[m][2] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][2].ip + l_eigen_Qjm[m][3] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][3].ip + l_eigen_Qjm[m][4] * Q[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][4].ip;
		}

		m++;
	}
}

void weno_solver_3()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_im[i] - 0.5*(roe_u_im[i] * roe_u_im[i] + roe_v_im[i] * roe_v_im[i] + roe_w_im[i] * roe_w_im[i])));

	//	if(hk == 3)
	{
		zeta_xi = metric[i].zeta_xim;
		zeta_yi = metric[i].zeta_yim;
		zeta_zi = metric[i].zeta_zim;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_im[i] + ky_bar * roe_v_im[i] + kz_bar * roe_w_im[i];
	phi_sq = 0.5*0.4*(roe_u_im[i] * roe_u_im[i] + roe_v_im[i] * roe_v_im[i] + roe_w_im[i] * roe_w_im[i]);
	alpha = roe_rho_im[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_im[i] * roe_a[i]);

	r_eigen_Qim[0][0] = kx_bar;
	r_eigen_Qim[0][1] = ky_bar;
	r_eigen_Qim[0][2] = kz_bar;
	r_eigen_Qim[0][3] = alpha;
	r_eigen_Qim[0][4] = alpha;

	r_eigen_Qim[1][0] = kx_bar * roe_u_im[i];
	r_eigen_Qim[1][1] = ky_bar * roe_u_im[i] - kz_bar * roe_rho_im[i];
	r_eigen_Qim[1][2] = kz_bar * roe_u_im[i] + ky_bar * roe_rho_im[i];
	r_eigen_Qim[1][3] = alpha * (roe_u_im[i] + kx_bar * roe_a[i]);
	r_eigen_Qim[1][4] = alpha * (roe_u_im[i] - kx_bar * roe_a[i]);

	r_eigen_Qim[2][0] = kx_bar * roe_v_im[i] + kz_bar * roe_rho_im[i];
	r_eigen_Qim[2][1] = ky_bar * roe_v_im[i];
	r_eigen_Qim[2][2] = kz_bar * roe_v_im[i] - kx_bar * roe_rho_im[i];
	r_eigen_Qim[2][3] = alpha * (roe_v_im[i] + ky_bar * roe_a[i]);
	r_eigen_Qim[2][4] = alpha * (roe_v_im[i] - ky_bar * roe_a[i]);

	r_eigen_Qim[3][0] = kx_bar * roe_w_im[i] - ky_bar * roe_rho_im[i];
	r_eigen_Qim[3][1] = ky_bar * roe_w_im[i] + kx_bar * roe_rho_im[i];
	r_eigen_Qim[3][2] = kz_bar * roe_w_im[i];
	r_eigen_Qim[3][3] = alpha * (roe_w_im[i] + kz_bar * roe_a[i]);
	r_eigen_Qim[3][4] = alpha * (roe_w_im[i] - kz_bar * roe_a[i]);

	r_eigen_Qim[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_im[i] * (kz_bar*roe_v_im[i] - ky_bar * roe_w_im[i]);
	r_eigen_Qim[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_im[i] * (kx_bar*roe_w_im[i] - kz_bar * roe_u_im[i]);
	r_eigen_Qim[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_im[i] * (ky_bar*roe_u_im[i] - kx_bar * roe_v_im[i]);
	r_eigen_Qim[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qim[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qim[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_im[i] - ky_bar * roe_w_im[i]) / roe_rho_im[i]);
	l_eigen_Qim[0][1] = kx_bar * 0.4*roe_u_im[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[0][2] = (kx_bar*0.4*roe_v_im[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_im[i]);
	l_eigen_Qim[0][3] = (kx_bar*0.4*roe_w_im[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_im[i]);
	l_eigen_Qim[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_im[i] - kz_bar * roe_u_im[i]) / roe_rho_im[i]);
	l_eigen_Qim[1][1] = (ky_bar*0.4*roe_u_im[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_im[i];
	l_eigen_Qim[1][2] = ky_bar * 0.4*roe_v_im[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[1][3] = (ky_bar*0.4*roe_w_im[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_im[i]);
	l_eigen_Qim[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_im[i] - kx_bar * roe_v_im[i]) / roe_rho_im[i]);
	l_eigen_Qim[2][1] = kz_bar * 0.4*roe_u_im[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_im[i]);
	l_eigen_Qim[2][2] = kz_bar * 0.4*roe_v_im[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_im[i]);
	l_eigen_Qim[2][3] = kz_bar * 0.4*roe_w_im[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qim[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qim[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qim[3][1] = -beta * (0.4*roe_u_im[i] - kx_bar * roe_a[i]);
	l_eigen_Qim[3][2] = -beta * (0.4*roe_v_im[i] - ky_bar * roe_a[i]);
	l_eigen_Qim[3][3] = -beta * (0.4*roe_w_im[i] - kz_bar * roe_a[i]);
	l_eigen_Qim[3][4] = beta * 0.4;

	l_eigen_Qim[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qim[4][1] = -beta * (0.4*roe_u_im[i] + kx_bar * roe_a[i]);
	l_eigen_Qim[4][2] = -beta * (0.4*roe_v_im[i] + ky_bar * roe_a[i]);
	l_eigen_Qim[4][3] = -beta * (0.4*roe_w_im[i] + kz_bar * roe_a[i]);
	l_eigen_Qim[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if(hk == 3 )
		{
			Qi_iminus_n[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][4].ip;
			Qi_iminus_n[node[node[i].n_n[1]].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[node[i].n_n[1]].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[node[i].n_n[1]].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[node[i].n_n[1]].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[node[i].n_n[1]].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[node[i].n_n[1]].n_n[1]][4].ip;
			Qi_iminus_n[node[i].n_n[1]][k] = l_eigen_Qim[m][0] * Q[node[i].n_n[1]][0].ip + l_eigen_Qim[m][1] * Q[node[i].n_n[1]][1].ip + l_eigen_Qim[m][2] * Q[node[i].n_n[1]][2].ip + l_eigen_Qim[m][3] * Q[node[i].n_n[1]][3].ip + l_eigen_Qim[m][4] * Q[node[i].n_n[1]][4].ip;
			Qi_iminus_n[i][k] = l_eigen_Qim[m][0] * Q[i][0].ip + l_eigen_Qim[m][1] * Q[i][1].ip + l_eigen_Qim[m][2] * Q[i][2].ip + l_eigen_Qim[m][3] * Q[i][3].ip + l_eigen_Qim[m][4] * Q[i][4].ip;
			Qi_iminus_n[node[i].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[i].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[i].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[i].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[i].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[i].n_n[3]][4].ip;
			Qi_iminus_n[node[node[i].n_n[3]].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[node[i].n_n[3]].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[node[i].n_n[3]].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[node[i].n_n[3]].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[node[i].n_n[3]].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[node[i].n_n[3]].n_n[3]][4].ip;
			Qi_iminus_n[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][k] = l_eigen_Qim[m][0] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][0].ip + l_eigen_Qim[m][1] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][1].ip + l_eigen_Qim[m][2] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][2].ip + l_eigen_Qim[m][3] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][3].ip + l_eigen_Qim[m][4] * Q[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][4].ip;
		}
		m++;
	}
}

void weno_solver_4()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_kp[i] - 0.5*(roe_u_kp[i] * roe_u_kp[i] + roe_v_kp[i] * roe_v_kp[i] + roe_w_kp[i] * roe_w_kp[i])));

	//	if (hk == 4)
	{
		zeta_xi = metric[i].xi_xkp;
		zeta_yi = metric[i].xi_ykp;
		zeta_zi = metric[i].xi_zkp;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_kp[i] + ky_bar * roe_v_kp[i] + kz_bar * roe_w_kp[i];
	phi_sq = 0.5*0.4*(roe_u_kp[i] * roe_u_kp[i] + roe_v_kp[i] * roe_v_kp[i] + roe_w_kp[i] * roe_w_kp[i]);
	alpha = roe_rho_kp[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_kp[i] * roe_a[i]);

	r_eigen_Qkp[0][0] = kx_bar;
	r_eigen_Qkp[0][1] = ky_bar;
	r_eigen_Qkp[0][2] = kz_bar;
	r_eigen_Qkp[0][3] = alpha;
	r_eigen_Qkp[0][4] = alpha;

	r_eigen_Qkp[1][0] = kx_bar * roe_u_kp[i];
	r_eigen_Qkp[1][1] = ky_bar * roe_u_kp[i] - kz_bar * roe_rho_kp[i];
	r_eigen_Qkp[1][2] = kz_bar * roe_u_kp[i] + ky_bar * roe_rho_kp[i];
	r_eigen_Qkp[1][3] = alpha * (roe_u_kp[i] + kx_bar * roe_a[i]);
	r_eigen_Qkp[1][4] = alpha * (roe_u_kp[i] - kx_bar * roe_a[i]);

	r_eigen_Qkp[2][0] = kx_bar * roe_v_kp[i] + kz_bar * roe_rho_kp[i];
	r_eigen_Qkp[2][1] = ky_bar * roe_v_kp[i];
	r_eigen_Qkp[2][2] = kz_bar * roe_v_kp[i] - kx_bar * roe_rho_kp[i];
	r_eigen_Qkp[2][3] = alpha * (roe_v_kp[i] + ky_bar * roe_a[i]);
	r_eigen_Qkp[2][4] = alpha * (roe_v_kp[i] - ky_bar * roe_a[i]);

	r_eigen_Qkp[3][0] = kx_bar * roe_w_kp[i] - ky_bar * roe_rho_kp[i];
	r_eigen_Qkp[3][1] = ky_bar * roe_w_kp[i] + kx_bar * roe_rho_kp[i];
	r_eigen_Qkp[3][2] = kz_bar * roe_w_kp[i];
	r_eigen_Qkp[3][3] = alpha * (roe_w_kp[i] + kz_bar * roe_a[i]);
	r_eigen_Qkp[3][4] = alpha * (roe_w_kp[i] - kz_bar * roe_a[i]);

	r_eigen_Qkp[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_kp[i] * (kz_bar*roe_v_kp[i] - ky_bar * roe_w_kp[i]);
	r_eigen_Qkp[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_kp[i] * (kx_bar*roe_w_kp[i] - kz_bar * roe_u_kp[i]);
	r_eigen_Qkp[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_kp[i] * (ky_bar*roe_u_kp[i] - kx_bar * roe_v_kp[i]);
	r_eigen_Qkp[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qkp[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qkp[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_kp[i] - ky_bar * roe_w_kp[i]) / roe_rho_kp[i]);
	l_eigen_Qkp[0][1] = kx_bar * 0.4*roe_u_kp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[0][2] = (kx_bar*0.4*roe_v_kp[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_kp[i]);
	l_eigen_Qkp[0][3] = (kx_bar*0.4*roe_w_kp[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_kp[i]);
	l_eigen_Qkp[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_kp[i] - kz_bar * roe_u_kp[i]) / roe_rho_kp[i]);
	l_eigen_Qkp[1][1] = (ky_bar*0.4*roe_u_kp[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_kp[i];
	l_eigen_Qkp[1][2] = ky_bar * 0.4*roe_v_kp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[1][3] = (ky_bar*0.4*roe_w_kp[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_kp[i]);
	l_eigen_Qkp[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_kp[i] - kx_bar * roe_v_kp[i]) / roe_rho_kp[i]);
	l_eigen_Qkp[2][1] = kz_bar * 0.4*roe_u_kp[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_kp[i]);
	l_eigen_Qkp[2][2] = kz_bar * 0.4*roe_v_kp[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_kp[i]);
	l_eigen_Qkp[2][3] = kz_bar * 0.4*roe_w_kp[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkp[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkp[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qkp[3][1] = -beta * (0.4*roe_u_kp[i] - kx_bar * roe_a[i]);
	l_eigen_Qkp[3][2] = -beta * (0.4*roe_v_kp[i] - ky_bar * roe_a[i]);
	l_eigen_Qkp[3][3] = -beta * (0.4*roe_w_kp[i] - kz_bar * roe_a[i]);
	l_eigen_Qkp[3][4] = beta * 0.4;

	l_eigen_Qkp[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qkp[4][1] = -beta * (0.4*roe_u_kp[i] + kx_bar * roe_a[i]);
	l_eigen_Qkp[4][2] = -beta * (0.4*roe_v_kp[i] + ky_bar * roe_a[i]);
	l_eigen_Qkp[4][3] = -beta * (0.4*roe_w_kp[i] + kz_bar * roe_a[i]);
	l_eigen_Qkp[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 4)
		{
			Qk_iplus[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;
			Qk_iplus[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[node[i].n_n[4]].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[i].n_n[4]].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[i].n_n[4]].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[i].n_n[4]].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[i].n_n[4]].n_n[4]][4].ip;
			Qk_iplus[node[i].n_n[4]][k] = l_eigen_Qkp[m][0] * Q[node[i].n_n[4]][0].ip + l_eigen_Qkp[m][1] * Q[node[i].n_n[4]][1].ip + l_eigen_Qkp[m][2] * Q[node[i].n_n[4]][2].ip + l_eigen_Qkp[m][3] * Q[node[i].n_n[4]][3].ip + l_eigen_Qkp[m][4] * Q[node[i].n_n[4]][4].ip;
			Qk_iplus[i][k] = l_eigen_Qkp[m][0] * Q[i][0].ip + l_eigen_Qkp[m][1] * Q[i][1].ip + l_eigen_Qkp[m][2] * Q[i][2].ip + l_eigen_Qkp[m][3] * Q[i][3].ip + l_eigen_Qkp[m][4] * Q[i][4].ip;
			Qk_iplus[node[i].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[i].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[i].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[i].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[i].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[i].n_n[5]][4].ip;
			Qk_iplus[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[node[i].n_n[5]].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[i].n_n[5]].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[i].n_n[5]].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[i].n_n[5]].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[i].n_n[5]].n_n[5]][4].ip;
			Qk_iplus[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_Qkp[m][0] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip + l_eigen_Qkp[m][1] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip + l_eigen_Qkp[m][2] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip + l_eigen_Qkp[m][3] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip + l_eigen_Qkp[m][4] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;
		}

		m++;
	}
}


void weno_solver_5()
{
	int k, m;
	double alpha, beta, theta, kk;
	double kx_bar, ky_bar, kz_bar;
	double phi_sq, zeta_xi, zeta_yi, zeta_zi;

	roe_a[i] = sqrt(0.4*(roe_h_km[i] - 0.5*(roe_u_km[i] * roe_u_km[i] + roe_v_km[i] * roe_v_km[i] + roe_w_km[i] * roe_w_km[i])));

	//	if(hk == 5)
	{
		zeta_xi = metric[i].xi_xkm;
		zeta_yi = metric[i].xi_ykm;
		zeta_zi = metric[i].xi_zkm;
	}

	kk = sqrt(zeta_xi*zeta_xi + zeta_yi * zeta_yi + zeta_zi * zeta_zi);

	kx_bar = zeta_xi / kk;
	ky_bar = zeta_yi / kk;
	kz_bar = zeta_zi / kk;

	theta = kx_bar * roe_u_km[i] + ky_bar * roe_v_km[i] + kz_bar * roe_w_km[i];
	phi_sq = 0.5*0.4*(roe_u_km[i] * roe_u_km[i] + roe_v_km[i] * roe_v_km[i] + roe_w_km[i] * roe_w_km[i]);
	alpha = roe_rho_km[i] / (sqrt(2.0)*roe_a[i]);
	beta = 1.0 / (sqrt(2.0)*roe_rho_km[i] * roe_a[i]);

	r_eigen_Qkm[0][0] = kx_bar;
	r_eigen_Qkm[0][1] = ky_bar;
	r_eigen_Qkm[0][2] = kz_bar;
	r_eigen_Qkm[0][3] = alpha;
	r_eigen_Qkm[0][4] = alpha;

	r_eigen_Qkm[1][0] = kx_bar * roe_u_km[i];
	r_eigen_Qkm[1][1] = ky_bar * roe_u_km[i] - kz_bar * roe_rho_km[i];
	r_eigen_Qkm[1][2] = kz_bar * roe_u_km[i] + ky_bar * roe_rho_km[i];
	r_eigen_Qkm[1][3] = alpha * (roe_u_km[i] + kx_bar * roe_a[i]);
	r_eigen_Qkm[1][4] = alpha * (roe_u_km[i] - kx_bar * roe_a[i]);

	r_eigen_Qkm[2][0] = kx_bar * roe_v_km[i] + kz_bar * roe_rho_km[i];
	r_eigen_Qkm[2][1] = ky_bar * roe_v_km[i];
	r_eigen_Qkm[2][2] = kz_bar * roe_v_km[i] - kx_bar * roe_rho_km[i];
	r_eigen_Qkm[2][3] = alpha * (roe_v_km[i] + ky_bar * roe_a[i]);
	r_eigen_Qkm[2][4] = alpha * (roe_v_km[i] - ky_bar * roe_a[i]);

	r_eigen_Qkm[3][0] = kx_bar * roe_w_km[i] - ky_bar * roe_rho_km[i];
	r_eigen_Qkm[3][1] = ky_bar * roe_w_km[i] + kx_bar * roe_rho_km[i];
	r_eigen_Qkm[3][2] = kz_bar * roe_w_km[i];
	r_eigen_Qkm[3][3] = alpha * (roe_w_km[i] + kz_bar * roe_a[i]);
	r_eigen_Qkm[3][4] = alpha * (roe_w_km[i] - kz_bar * roe_a[i]);

	r_eigen_Qkm[4][0] = ((kx_bar*phi_sq) / 0.4) + roe_rho_km[i] * (kz_bar*roe_v_km[i] - ky_bar * roe_w_km[i]);
	r_eigen_Qkm[4][1] = ((ky_bar*phi_sq) / 0.4) + roe_rho_km[i] * (kx_bar*roe_w_km[i] - kz_bar * roe_u_km[i]);
	r_eigen_Qkm[4][2] = ((kz_bar*phi_sq) / 0.4) + roe_rho_km[i] * (ky_bar*roe_u_km[i] - kx_bar * roe_v_km[i]);
	r_eigen_Qkm[4][3] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) + theta * roe_a[i]);
	r_eigen_Qkm[4][4] = alpha * (((phi_sq + roe_a[i] * roe_a[i]) / (0.4)) - theta * roe_a[i]);

	/*********************************************/

	l_eigen_Qkm[0][0] = kx_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kz_bar*roe_v_km[i] - ky_bar * roe_w_km[i]) / roe_rho_km[i]);
	l_eigen_Qkm[0][1] = kx_bar * 0.4*roe_u_km[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[0][2] = (kx_bar*0.4*roe_v_km[i] / (roe_a[i] * roe_a[i])) + (kz_bar / roe_rho_km[i]);
	l_eigen_Qkm[0][3] = (kx_bar*0.4*roe_w_km[i] / (roe_a[i] * roe_a[i])) - (ky_bar / roe_rho_km[i]);
	l_eigen_Qkm[0][4] = -kx_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[1][0] = ky_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((kx_bar*roe_w_km[i] - kz_bar * roe_u_km[i]) / roe_rho_km[i]);
	l_eigen_Qkm[1][1] = (ky_bar*0.4*roe_u_km[i] / (roe_a[i] * roe_a[i])) - kz_bar / roe_rho_km[i];
	l_eigen_Qkm[1][2] = ky_bar * 0.4*roe_v_km[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[1][3] = (ky_bar*0.4*roe_w_km[i] / (roe_a[i] * roe_a[i])) + (kx_bar / roe_rho_km[i]);
	l_eigen_Qkm[1][4] = -ky_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[2][0] = kz_bar * (1.0 - (phi_sq / (roe_a[i] * roe_a[i]))) - ((ky_bar*roe_u_km[i] - kx_bar * roe_v_km[i]) / roe_rho_km[i]);
	l_eigen_Qkm[2][1] = kz_bar * 0.4*roe_u_km[i] / (roe_a[i] * roe_a[i]) + (ky_bar / roe_rho_km[i]);
	l_eigen_Qkm[2][2] = kz_bar * 0.4*roe_v_km[i] / (roe_a[i] * roe_a[i]) - (kx_bar / roe_rho_km[i]);
	l_eigen_Qkm[2][3] = kz_bar * 0.4*roe_w_km[i] / (roe_a[i] * roe_a[i]);
	l_eigen_Qkm[2][4] = -kz_bar * 0.4 / (roe_a[i] * roe_a[i]);

	l_eigen_Qkm[3][0] = beta * (phi_sq - theta * roe_a[i]);
	l_eigen_Qkm[3][1] = -beta * (0.4*roe_u_km[i] - kx_bar * roe_a[i]);
	l_eigen_Qkm[3][2] = -beta * (0.4*roe_v_km[i] - ky_bar * roe_a[i]);
	l_eigen_Qkm[3][3] = -beta * (0.4*roe_w_km[i] - kz_bar * roe_a[i]);
	l_eigen_Qkm[3][4] = beta * 0.4;

	l_eigen_Qkm[4][0] = beta * (phi_sq + theta * roe_a[i]);
	l_eigen_Qkm[4][1] = -beta * (0.4*roe_u_km[i] + kx_bar * roe_a[i]);
	l_eigen_Qkm[4][2] = -beta * (0.4*roe_v_km[i] + ky_bar * roe_a[i]);
	l_eigen_Qkm[4][3] = -beta * (0.4*roe_w_km[i] + kz_bar * roe_a[i]);
	l_eigen_Qkm[4][4] = beta * 0.4;

	/****************************Local Lax-Friedrichs(LLF) scheme*******************************/

	m = 0;
	for (k = 0; k<5; k++)
	{
		//	if (hk == 5)
		{
			Qk_iminus_n[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][4].ip;
			Qk_iminus_n[node[node[i].n_n[4]].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[node[i].n_n[4]].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[i].n_n[4]].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[i].n_n[4]].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[i].n_n[4]].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[i].n_n[4]].n_n[4]][4].ip;
			Qk_iminus_n[node[i].n_n[4]][k] = l_eigen_Qkm[m][0] * Q[node[i].n_n[4]][0].ip + l_eigen_Qkm[m][1] * Q[node[i].n_n[4]][1].ip + l_eigen_Qkm[m][2] * Q[node[i].n_n[4]][2].ip + l_eigen_Qkm[m][3] * Q[node[i].n_n[4]][3].ip + l_eigen_Qkm[m][4] * Q[node[i].n_n[4]][4].ip;
			Qk_iminus_n[i][k] = l_eigen_Qkm[m][0] * Q[i][0].ip + l_eigen_Qkm[m][1] * Q[i][1].ip + l_eigen_Qkm[m][2] * Q[i][2].ip + l_eigen_Qkm[m][3] * Q[i][3].ip + l_eigen_Qkm[m][4] * Q[i][4].ip;
			Qk_iminus_n[node[i].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[i].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[i].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[i].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[i].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[i].n_n[5]][4].ip;
			Qk_iminus_n[node[node[i].n_n[5]].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[node[i].n_n[5]].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[i].n_n[5]].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[i].n_n[5]].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[i].n_n[5]].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[i].n_n[5]].n_n[5]][4].ip;
			Qk_iminus_n[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][k] = l_eigen_Qkm[m][0] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][0].ip + l_eigen_Qkm[m][1] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][1].ip + l_eigen_Qkm[m][2] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][2].ip + l_eigen_Qkm[m][3] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][3].ip + l_eigen_Qkm[m][4] * Q[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][4].ip;
		}

		m++;
	}
}

