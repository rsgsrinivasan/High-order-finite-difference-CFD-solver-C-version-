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

using std::cout;
using std::cin;
using std::string;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::istringstream;

void viscousflux_variables()
{
	double t_zeta, t_eta, t_xi, u_zeta, u_eta, u_xi;
	double v_zeta, v_eta, v_xi, w_zeta, w_eta, w_xi;
	double U1_zeta, U2_zeta, U3_zeta, U4_zeta, U5_zeta, U1_eta, U2_eta, U3_eta, U4_eta, U5_eta, U1_xi, U2_xi, U3_xi, U4_xi, U5_xi;
	double d_temp, d_temp1;
	double stress_xy, stress_yz, stress_xz, stress_xx, stress_yy, stress_zz;
//	double tau_sgs_xx, tau_sgs_yy, tau_sgs_zz;
	double KS_x, KS_y, KS_z;
//	double TS_x, TS_y, TS_z;
	double EDDY_x, EDDY_y, EDDY_z, mod_stress;
	double DEL_ZETA, DEL_ETA, DEL_XI;
	double Q_x, Q_y, Q_z, DEL_X, DEL_Y, DEL_Z;
	double C_MTS, C_T, Rx, Ry, Rz;
	double stress_xy_star, stress_yz_star, stress_xz_star;
	double divergence_V, cross_V;

	t_zeta = 0.0;
	t_eta = 0.0;
	t_xi = 0.0;
	u_zeta = 0.0;
	u_eta = 0.0;
	u_xi = 0.0;
	v_zeta = 0.0;
	v_eta = 0.0;
	v_xi = 0.0;
	w_zeta = 0.0;
	w_eta = 0.0;
	w_xi = 0.0;
	U1_zeta = 0.0;
	U2_zeta = 0.0;
	U3_zeta = 0.0;
	U4_zeta = 0.0;
	U5_zeta = 0.0;
	U1_eta = 0.0;
	U2_eta = 0.0;
	U3_eta = 0.0;
	U4_eta = 0.0;
	U5_eta = 0.0;
	U1_xi = 0.0;
	U2_xi = 0.0;
	U3_xi = 0.0;
	U4_xi = 0.0;
	U5_xi = 0.0;
	d_temp = 0.0;
	d_temp1 = 0.0;

	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[1]].loc == 100))
	{
		U1_zeta = (1.0 / 6.0)*(11.0*U[i][j][0] - 18.0*U[node[i].n_n[3]][j][0] + 9.0*U[node[node[i].n_n[3]].n_n[3]][j][0] - 2.0*U[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][j][0])*(1.0 / del_eta);
		U2_zeta = (1.0 / 6.0)*(11.0*U[i][j][1] - 18.0*U[node[i].n_n[3]][j][1] + 9.0*U[node[node[i].n_n[3]].n_n[3]][j][1] - 2.0*U[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][j][1])*(1.0 / del_eta);
		U3_zeta = (1.0 / 6.0)*(11.0*U[i][j][2] - 18.0*U[node[i].n_n[3]][j][2] + 9.0*U[node[node[i].n_n[3]].n_n[3]][j][2] - 2.0*U[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][j][2])*(1.0 / del_eta);
		U4_zeta = (1.0 / 6.0)*(11.0*U[i][j][3] - 18.0*U[node[i].n_n[3]][j][3] + 9.0*U[node[node[i].n_n[3]].n_n[3]][j][3] - 2.0*U[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][j][3])*(1.0 / del_eta);
		U5_zeta = (1.0 / 6.0)*(11.0*U[i][j][4] - 18.0*U[node[i].n_n[3]][j][4] + 9.0*U[node[node[i].n_n[3]].n_n[3]][j][4] - 2.0*U[node[node[node[i].n_n[3]].n_n[3]].n_n[3]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[3]].loc == 100))
	{
		U1_zeta = (1.0 / 6.0)*(-11.0*U[i][j][0] + 18.0*U[node[i].n_n[1]][j][0] - 9.0*U[node[node[i].n_n[1]].n_n[1]][j][0] + 2.0*U[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][j][0])*(1.0 / del_eta);
		U2_zeta = (1.0 / 6.0)*(-11.0*U[i][j][1] + 18.0*U[node[i].n_n[1]][j][1] - 9.0*U[node[node[i].n_n[1]].n_n[1]][j][1] + 2.0*U[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][j][1])*(1.0 / del_eta);
		U3_zeta = (1.0 / 6.0)*(-11.0*U[i][j][2] + 18.0*U[node[i].n_n[1]][j][2] - 9.0*U[node[node[i].n_n[1]].n_n[1]][j][2] + 2.0*U[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][j][2])*(1.0 / del_eta);
		U4_zeta = (1.0 / 6.0)*(-11.0*U[i][j][3] + 18.0*U[node[i].n_n[1]][j][3] - 9.0*U[node[node[i].n_n[1]].n_n[1]][j][3] + 2.0*U[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][j][3])*(1.0 / del_eta);
		U5_zeta = (1.0 / 6.0)*(-11.0*U[i][j][4] + 18.0*U[node[i].n_n[1]][j][4] - 9.0*U[node[node[i].n_n[1]].n_n[1]][j][4] + 2.0*U[node[node[node[i].n_n[1]].n_n[1]].n_n[1]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[1]].loc > 0 && node[node[i].n_n[1]].loc <= 52 && (node[node[node[i].n_n[1]].n_n[1]].loc == 100))
	{
		U1_zeta = (1.0 / 6.0)*(2.0*U[node[i].n_n[1]][j][0] + 3.0*U[i][j][0] - 6.0*U[node[i].n_n[3]][j][0] + U[node[node[i].n_n[3]].n_n[3]][j][0])*(1.0 / del_eta);
		U2_zeta = (1.0 / 6.0)*(2.0*U[node[i].n_n[1]][j][1] + 3.0*U[i][j][1] - 6.0*U[node[i].n_n[3]][j][1] + U[node[node[i].n_n[3]].n_n[3]][j][1])*(1.0 / del_eta);
		U3_zeta = (1.0 / 6.0)*(2.0*U[node[i].n_n[1]][j][2] + 3.0*U[i][j][2] - 6.0*U[node[i].n_n[3]][j][2] + U[node[node[i].n_n[3]].n_n[3]][j][2])*(1.0 / del_eta);
		U4_zeta = (1.0 / 6.0)*(2.0*U[node[i].n_n[1]][j][3] + 3.0*U[i][j][3] - 6.0*U[node[i].n_n[3]][j][3] + U[node[node[i].n_n[3]].n_n[3]][j][3])*(1.0 / del_eta);
		U5_zeta = (1.0 / 6.0)*(2.0*U[node[i].n_n[1]][j][4] + 3.0*U[i][j][4] - 6.0*U[node[i].n_n[3]][j][4] + U[node[node[i].n_n[3]].n_n[3]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[3]].loc > 0 && node[node[i].n_n[3]].loc <= 52 && (node[node[node[i].n_n[3]].n_n[3]].loc == 100))
	{
		U1_zeta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[3]][j][0] - 3.0*U[i][j][0] + 6.0*U[node[i].n_n[1]][j][0] - U[node[node[i].n_n[1]].n_n[1]][j][0])*(1.0 / del_zeta);
		U2_zeta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[3]][j][1] - 3.0*U[i][j][1] + 6.0*U[node[i].n_n[1]][j][1] - U[node[node[i].n_n[1]].n_n[1]][j][1])*(1.0 / del_zeta);
		U3_zeta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[3]][j][2] - 3.0*U[i][j][2] + 6.0*U[node[i].n_n[1]][j][2] - U[node[node[i].n_n[1]].n_n[1]][j][2])*(1.0 / del_zeta);
		U4_zeta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[3]][j][3] - 3.0*U[i][j][3] + 6.0*U[node[i].n_n[1]][j][3] - U[node[node[i].n_n[1]].n_n[1]][j][3])*(1.0 / del_zeta);
		U5_zeta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[3]][j][4] - 3.0*U[i][j][4] + 6.0*U[node[i].n_n[1]][j][4] - U[node[node[i].n_n[1]].n_n[1]][j][4])*(1.0 / del_zeta);
	}
	if (node[i].loc == 0 && node[node[node[i].n_n[3]].n_n[3]].loc != 100 && node[node[i].n_n[3]].loc != 100 && node[node[node[i].n_n[1]].n_n[1]].loc != 100 && node[node[i].n_n[1]].loc != 100)
	{
		U1_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][0] - U[node[i].n_n[3]][j][0]) - (U[node[node[i].n_n[1]].n_n[1]][j][0] - U[node[node[i].n_n[3]].n_n[3]][j][0]))*(1.0 / del_eta);
		U2_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][1] - U[node[i].n_n[3]][j][1]) - (U[node[node[i].n_n[1]].n_n[1]][j][1] - U[node[node[i].n_n[3]].n_n[3]][j][1]))*(1.0 / del_eta);
		U3_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][2] - U[node[i].n_n[3]][j][2]) - (U[node[node[i].n_n[1]].n_n[1]][j][2] - U[node[node[i].n_n[3]].n_n[3]][j][2]))*(1.0 / del_eta);
		U4_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][3] - U[node[i].n_n[3]][j][3]) - (U[node[node[i].n_n[1]].n_n[1]][j][3] - U[node[node[i].n_n[3]].n_n[3]][j][3]))*(1.0 / del_eta);
		U5_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][4] - U[node[i].n_n[3]][j][4]) - (U[node[node[i].n_n[1]].n_n[1]][j][4] - U[node[node[i].n_n[3]].n_n[3]][j][4]))*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[node[node[i].n_n[3]].n_n[3]].loc != 100 && node[node[i].n_n[3]].loc != 100 && node[node[node[i].n_n[1]].n_n[1]].loc != 100 && node[node[i].n_n[1]].loc != 100)
	{
		U1_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][0] - U[node[i].n_n[3]][j][0]) - (U[node[node[i].n_n[1]].n_n[1]][j][0] - U[node[node[i].n_n[3]].n_n[3]][j][0]))*(1.0 / del_eta);
		U2_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][1] - U[node[i].n_n[3]][j][1]) - (U[node[node[i].n_n[1]].n_n[1]][j][1] - U[node[node[i].n_n[3]].n_n[3]][j][1]))*(1.0 / del_eta);
		U3_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][2] - U[node[i].n_n[3]][j][2]) - (U[node[node[i].n_n[1]].n_n[1]][j][2] - U[node[node[i].n_n[3]].n_n[3]][j][2]))*(1.0 / del_eta);
		U4_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][3] - U[node[i].n_n[3]][j][3]) - (U[node[node[i].n_n[1]].n_n[1]][j][3] - U[node[node[i].n_n[3]].n_n[3]][j][3]))*(1.0 / del_eta);
		U5_zeta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[1]][j][4] - U[node[i].n_n[3]][j][4]) - (U[node[node[i].n_n[1]].n_n[1]][j][4] - U[node[node[i].n_n[3]].n_n[3]][j][4]))*(1.0 / del_eta);
	}



	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[0]].loc == 100))
	{
		U1_eta = (1.0 / 6.0)*(11.0*U[i][j][0] - 18.0*U[node[i].n_n[2]][j][0] + 9.0*U[node[node[i].n_n[2]].n_n[2]][j][0] - 2.0*U[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][j][0])*(1.0 / del_eta);
		U2_eta = (1.0 / 6.0)*(11.0*U[i][j][1] - 18.0*U[node[i].n_n[2]][j][1] + 9.0*U[node[node[i].n_n[2]].n_n[2]][j][1] - 2.0*U[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][j][1])*(1.0 / del_eta);
		U3_eta = (1.0 / 6.0)*(11.0*U[i][j][2] - 18.0*U[node[i].n_n[2]][j][2] + 9.0*U[node[node[i].n_n[2]].n_n[2]][j][2] - 2.0*U[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][j][2])*(1.0 / del_eta);
		U4_eta = (1.0 / 6.0)*(11.0*U[i][j][3] - 18.0*U[node[i].n_n[2]][j][3] + 9.0*U[node[node[i].n_n[2]].n_n[2]][j][3] - 2.0*U[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][j][3])*(1.0 / del_eta);
		U5_eta = (1.0 / 6.0)*(11.0*U[i][j][4] - 18.0*U[node[i].n_n[2]][j][4] + 9.0*U[node[node[i].n_n[2]].n_n[2]][j][4] - 2.0*U[node[node[node[i].n_n[2]].n_n[2]].n_n[2]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[2]].loc == 100))
	{
		U1_eta = (1.0 / 6.0)*(-11.0*U[i][j][0] + 18.0*U[node[i].n_n[0]][j][0] - 9.0*U[node[node[i].n_n[0]].n_n[0]][j][0] + 2.0*U[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][j][0])*(1.0 / del_eta);
		U2_eta = (1.0 / 6.0)*(-11.0*U[i][j][1] + 18.0*U[node[i].n_n[0]][j][1] - 9.0*U[node[node[i].n_n[0]].n_n[0]][j][1] + 2.0*U[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][j][1])*(1.0 / del_eta);
		U3_eta = (1.0 / 6.0)*(-11.0*U[i][j][2] + 18.0*U[node[i].n_n[0]][j][2] - 9.0*U[node[node[i].n_n[0]].n_n[0]][j][2] + 2.0*U[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][j][2])*(1.0 / del_eta);
		U4_eta = (1.0 / 6.0)*(-11.0*U[i][j][3] + 18.0*U[node[i].n_n[0]][j][3] - 9.0*U[node[node[i].n_n[0]].n_n[0]][j][3] + 2.0*U[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][j][3])*(1.0 / del_eta);
		U5_eta = (1.0 / 6.0)*(-11.0*U[i][j][4] + 18.0*U[node[i].n_n[0]][j][4] - 9.0*U[node[node[i].n_n[0]].n_n[0]][j][4] + 2.0*U[node[node[node[i].n_n[0]].n_n[0]].n_n[0]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[0]].loc > 0 && node[node[i].n_n[0]].loc <= 52 && (node[node[node[i].n_n[0]].n_n[0]].loc == 100))
	{
		U1_eta = (1.0 / 6.0)*(2.0*U[node[i].n_n[0]][j][0] + 3.0*U[i][j][0] - 6.0*U[node[i].n_n[2]][j][0] + U[node[node[i].n_n[2]].n_n[2]][j][0])*(1.0 / del_eta);
		U2_eta = (1.0 / 6.0)*(2.0*U[node[i].n_n[0]][j][1] + 3.0*U[i][j][1] - 6.0*U[node[i].n_n[2]][j][1] + U[node[node[i].n_n[2]].n_n[2]][j][1])*(1.0 / del_eta);
		U3_eta = (1.0 / 6.0)*(2.0*U[node[i].n_n[0]][j][2] + 3.0*U[i][j][2] - 6.0*U[node[i].n_n[2]][j][2] + U[node[node[i].n_n[2]].n_n[2]][j][2])*(1.0 / del_eta);
		U4_eta = (1.0 / 6.0)*(2.0*U[node[i].n_n[0]][j][3] + 3.0*U[i][j][3] - 6.0*U[node[i].n_n[2]][j][3] + U[node[node[i].n_n[2]].n_n[2]][j][3])*(1.0 / del_eta);
		U5_eta = (1.0 / 6.0)*(2.0*U[node[i].n_n[0]][j][4] + 3.0*U[i][j][4] - 6.0*U[node[i].n_n[2]][j][4] + U[node[node[i].n_n[2]].n_n[2]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[2]].loc > 0 && node[node[i].n_n[2]].loc <= 52 && (node[node[node[i].n_n[2]].n_n[2]].loc == 100))
	{
		U1_eta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[2]][j][0] - 3.0*U[i][j][0] + 6.0*U[node[i].n_n[0]][j][0] - U[node[node[i].n_n[0]].n_n[0]][j][0])*(1.0 / del_zeta);
		U2_eta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[2]][j][1] - 3.0*U[i][j][1] + 6.0*U[node[i].n_n[0]][j][1] - U[node[node[i].n_n[0]].n_n[0]][j][1])*(1.0 / del_zeta);
		U3_eta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[2]][j][2] - 3.0*U[i][j][2] + 6.0*U[node[i].n_n[0]][j][2] - U[node[node[i].n_n[0]].n_n[0]][j][2])*(1.0 / del_zeta);
		U4_eta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[2]][j][3] - 3.0*U[i][j][3] + 6.0*U[node[i].n_n[0]][j][3] - U[node[node[i].n_n[0]].n_n[0]][j][3])*(1.0 / del_zeta);
		U5_eta = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[2]][j][4] - 3.0*U[i][j][4] + 6.0*U[node[i].n_n[0]][j][4] - U[node[node[i].n_n[0]].n_n[0]][j][4])*(1.0 / del_zeta);
	}
	if (node[i].loc == 0 && node[node[node[i].n_n[2]].n_n[2]].loc != 100 && node[node[i].n_n[2]].loc != 100 && node[node[node[i].n_n[0]].n_n[0]].loc != 100 && node[node[i].n_n[0]].loc != 100)
	{
		U1_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][0] - U[node[i].n_n[2]][j][0]) - (U[node[node[i].n_n[0]].n_n[0]][j][0] - U[node[node[i].n_n[2]].n_n[2]][j][0]))*(1.0 / del_eta);
		U2_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][1] - U[node[i].n_n[2]][j][1]) - (U[node[node[i].n_n[0]].n_n[0]][j][1] - U[node[node[i].n_n[2]].n_n[2]][j][1]))*(1.0 / del_eta);
		U3_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][2] - U[node[i].n_n[2]][j][2]) - (U[node[node[i].n_n[0]].n_n[0]][j][2] - U[node[node[i].n_n[2]].n_n[2]][j][2]))*(1.0 / del_eta);
		U4_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][3] - U[node[i].n_n[2]][j][3]) - (U[node[node[i].n_n[0]].n_n[0]][j][3] - U[node[node[i].n_n[2]].n_n[2]][j][3]))*(1.0 / del_eta);
		U5_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][4] - U[node[i].n_n[2]][j][4]) - (U[node[node[i].n_n[0]].n_n[0]][j][4] - U[node[node[i].n_n[2]].n_n[2]][j][4]))*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[node[node[i].n_n[2]].n_n[2]].loc != 100 && node[node[i].n_n[2]].loc != 100 && node[node[node[i].n_n[0]].n_n[0]].loc != 100 && node[node[i].n_n[0]].loc != 100)
	{
		U1_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][0] - U[node[i].n_n[2]][j][0]) - (U[node[node[i].n_n[0]].n_n[0]][j][0] - U[node[node[i].n_n[2]].n_n[2]][j][0]))*(1.0 / del_eta);
		U2_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][1] - U[node[i].n_n[2]][j][1]) - (U[node[node[i].n_n[0]].n_n[0]][j][1] - U[node[node[i].n_n[2]].n_n[2]][j][1]))*(1.0 / del_eta);
		U3_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][2] - U[node[i].n_n[2]][j][2]) - (U[node[node[i].n_n[0]].n_n[0]][j][2] - U[node[node[i].n_n[2]].n_n[2]][j][2]))*(1.0 / del_eta);
		U4_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][3] - U[node[i].n_n[2]][j][3]) - (U[node[node[i].n_n[0]].n_n[0]][j][3] - U[node[node[i].n_n[2]].n_n[2]][j][3]))*(1.0 / del_eta);
		U5_eta = (1.0 / 12.0)*(8.0*(U[node[i].n_n[0]][j][4] - U[node[i].n_n[2]][j][4]) - (U[node[node[i].n_n[0]].n_n[0]][j][4] - U[node[node[i].n_n[2]].n_n[2]][j][4]))*(1.0 / del_eta);
	}



	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[4]].loc == 100))
	{
		U1_xi = (1.0 / 6.0)*(11.0*U[i][j][0] - 18.0*U[node[i].n_n[5]][j][0] + 9.0*U[node[node[i].n_n[5]].n_n[5]][j][0] - 2.0*U[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][j][0])*(1.0 / del_eta);
		U2_xi = (1.0 / 6.0)*(11.0*U[i][j][1] - 18.0*U[node[i].n_n[5]][j][1] + 9.0*U[node[node[i].n_n[5]].n_n[5]][j][1] - 2.0*U[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][j][1])*(1.0 / del_eta);
		U3_xi = (1.0 / 6.0)*(11.0*U[i][j][2] - 18.0*U[node[i].n_n[5]][j][2] + 9.0*U[node[node[i].n_n[5]].n_n[5]][j][2] - 2.0*U[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][j][2])*(1.0 / del_eta);
		U4_xi = (1.0 / 6.0)*(11.0*U[i][j][3] - 18.0*U[node[i].n_n[5]][j][3] + 9.0*U[node[node[i].n_n[5]].n_n[5]][j][3] - 2.0*U[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][j][3])*(1.0 / del_eta);
		U5_xi = (1.0 / 6.0)*(11.0*U[i][j][4] - 18.0*U[node[i].n_n[5]][j][4] + 9.0*U[node[node[i].n_n[5]].n_n[5]][j][4] - 2.0*U[node[node[node[i].n_n[5]].n_n[5]].n_n[5]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[i].loc <= 52 && (node[node[i].n_n[5]].loc == 100))
	{
		U1_xi = (1.0 / 6.0)*(-11.0*U[i][j][0] + 18.0*U[node[i].n_n[4]][j][0] - 9.0*U[node[node[i].n_n[4]].n_n[4]][j][0] + 2.0*U[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][j][0])*(1.0 / del_eta);
		U2_xi = (1.0 / 6.0)*(-11.0*U[i][j][1] + 18.0*U[node[i].n_n[4]][j][1] - 9.0*U[node[node[i].n_n[4]].n_n[4]][j][1] + 2.0*U[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][j][1])*(1.0 / del_eta);
		U3_xi = (1.0 / 6.0)*(-11.0*U[i][j][2] + 18.0*U[node[i].n_n[4]][j][2] - 9.0*U[node[node[i].n_n[4]].n_n[4]][j][2] + 2.0*U[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][j][2])*(1.0 / del_eta);
		U4_xi = (1.0 / 6.0)*(-11.0*U[i][j][3] + 18.0*U[node[i].n_n[4]][j][3] - 9.0*U[node[node[i].n_n[4]].n_n[4]][j][3] + 2.0*U[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][j][3])*(1.0 / del_eta);
		U5_xi = (1.0 / 6.0)*(-11.0*U[i][j][4] + 18.0*U[node[i].n_n[4]][j][4] - 9.0*U[node[node[i].n_n[4]].n_n[4]][j][4] + 2.0*U[node[node[node[i].n_n[4]].n_n[4]].n_n[4]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[4]].loc > 0 && node[node[i].n_n[4]].loc <= 52 && (node[node[node[i].n_n[4]].n_n[4]].loc == 100))
	{
		U1_xi = (1.0 / 6.0)*(2.0*U[node[i].n_n[4]][j][0] + 3.0*U[i][j][0] - 6.0*U[node[i].n_n[5]][j][0] + U[node[node[i].n_n[5]].n_n[5]][j][0])*(1.0 / del_eta);
		U2_xi = (1.0 / 6.0)*(2.0*U[node[i].n_n[4]][j][1] + 3.0*U[i][j][1] - 6.0*U[node[i].n_n[5]][j][1] + U[node[node[i].n_n[5]].n_n[5]][j][1])*(1.0 / del_eta);
		U3_xi = (1.0 / 6.0)*(2.0*U[node[i].n_n[4]][j][2] + 3.0*U[i][j][2] - 6.0*U[node[i].n_n[5]][j][2] + U[node[node[i].n_n[5]].n_n[5]][j][2])*(1.0 / del_eta);
		U4_xi = (1.0 / 6.0)*(2.0*U[node[i].n_n[4]][j][3] + 3.0*U[i][j][3] - 6.0*U[node[i].n_n[5]][j][3] + U[node[node[i].n_n[5]].n_n[5]][j][3])*(1.0 / del_eta);
		U5_xi = (1.0 / 6.0)*(2.0*U[node[i].n_n[4]][j][4] + 3.0*U[i][j][4] - 6.0*U[node[i].n_n[5]][j][4] + U[node[node[i].n_n[5]].n_n[5]][j][4])*(1.0 / del_eta);
	}
	if (node[i].loc == 0 && node[node[i].n_n[5]].loc > 0 && node[node[i].n_n[5]].loc <= 52 && (node[node[node[i].n_n[5]].n_n[5]].loc == 100))
	{
		U1_xi = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[5]][j][0] - 3.0*U[i][j][0] + 6.0*U[node[i].n_n[4]][j][0] - U[node[node[i].n_n[4]].n_n[4]][j][0])*(1.0 / del_zeta);
		U2_xi = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[5]][j][1] - 3.0*U[i][j][1] + 6.0*U[node[i].n_n[4]][j][1] - U[node[node[i].n_n[4]].n_n[4]][j][1])*(1.0 / del_zeta);
		U3_xi = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[5]][j][2] - 3.0*U[i][j][2] + 6.0*U[node[i].n_n[4]][j][2] - U[node[node[i].n_n[4]].n_n[4]][j][2])*(1.0 / del_zeta);
		U4_xi = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[5]][j][3] - 3.0*U[i][j][3] + 6.0*U[node[i].n_n[4]][j][3] - U[node[node[i].n_n[4]].n_n[4]][j][3])*(1.0 / del_zeta);
		U5_xi = (1.0 / 6.0)*((-2.0)*U[node[i].n_n[5]][j][4] - 3.0*U[i][j][4] + 6.0*U[node[i].n_n[4]][j][4] - U[node[node[i].n_n[4]].n_n[4]][j][4])*(1.0 / del_zeta);
	}
	if (node[i].loc == 0 && node[node[node[i].n_n[5]].n_n[5]].loc != 100 && node[node[i].n_n[5]].loc != 100 && node[node[node[i].n_n[4]].n_n[4]].loc != 100 && node[node[i].n_n[4]].loc != 100)
	{
		U1_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][0] - U[node[i].n_n[5]][j][0]) - (U[node[node[i].n_n[4]].n_n[4]][j][0] - U[node[node[i].n_n[5]].n_n[5]][j][0]))*(1.0 / del_eta);
		U2_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][1] - U[node[i].n_n[5]][j][1]) - (U[node[node[i].n_n[4]].n_n[4]][j][1] - U[node[node[i].n_n[5]].n_n[5]][j][1]))*(1.0 / del_eta);
		U3_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][2] - U[node[i].n_n[5]][j][2]) - (U[node[node[i].n_n[4]].n_n[4]][j][2] - U[node[node[i].n_n[5]].n_n[5]][j][2]))*(1.0 / del_eta);
		U4_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][3] - U[node[i].n_n[5]][j][3]) - (U[node[node[i].n_n[4]].n_n[4]][j][3] - U[node[node[i].n_n[5]].n_n[5]][j][3]))*(1.0 / del_eta);
		U5_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][4] - U[node[i].n_n[5]][j][4]) - (U[node[node[i].n_n[4]].n_n[4]][j][4] - U[node[node[i].n_n[5]].n_n[5]][j][4]))*(1.0 / del_eta);
	}
	if (node[i].loc > 0 && node[node[node[i].n_n[5]].n_n[5]].loc != 100 && node[node[i].n_n[5]].loc != 100 && node[node[node[i].n_n[4]].n_n[4]].loc != 100 && node[node[i].n_n[4]].loc != 100)
	{
		U1_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][0] - U[node[i].n_n[5]][j][0]) - (U[node[node[i].n_n[4]].n_n[4]][j][0] - U[node[node[i].n_n[5]].n_n[5]][j][0]))*(1.0 / del_eta);
		U2_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][1] - U[node[i].n_n[5]][j][1]) - (U[node[node[i].n_n[4]].n_n[4]][j][1] - U[node[node[i].n_n[5]].n_n[5]][j][1]))*(1.0 / del_eta);
		U3_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][2] - U[node[i].n_n[5]][j][2]) - (U[node[node[i].n_n[4]].n_n[4]][j][2] - U[node[node[i].n_n[5]].n_n[5]][j][2]))*(1.0 / del_eta);
		U4_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][3] - U[node[i].n_n[5]][j][3]) - (U[node[node[i].n_n[4]].n_n[4]][j][3] - U[node[node[i].n_n[5]].n_n[5]][j][3]))*(1.0 / del_eta);
		U5_xi = (1.0 / 12.0)*(8.0*(U[node[i].n_n[4]][j][4] - U[node[i].n_n[5]][j][4]) - (U[node[node[i].n_n[4]].n_n[4]][j][4] - U[node[node[i].n_n[5]].n_n[5]][j][4]))*(1.0 / del_eta);
	}


	u_zeta = (1.0 / U[i][j][0])*(U2_zeta)-(U[i][j][1] / ((U[i][j][0] * U[i][j][0])))*(U1_zeta);
	u_eta = (1.0 / U[i][j][0])*(U2_eta)-(U[i][j][1] / ((U[i][j][0] * U[i][j][0])))*(U1_eta);
	u_xi = (1.0 / U[i][j][0])*(U2_xi)-(U[i][j][1] / ((U[i][j][0] * U[i][j][0])))*(U1_xi);

	v_zeta = (1.0 / U[i][j][0])*(U3_zeta)-(U[i][j][2] / ((U[i][j][0] * U[i][j][0])))*(U1_zeta);
	v_eta = (1.0 / U[i][j][0])*(U3_eta)-(U[i][j][2] / ((U[i][j][0] * U[i][j][0])))*(U1_eta);
	v_xi = (1.0 / U[i][j][0])*(U3_xi)-(U[i][j][2] / ((U[i][j][0] * U[i][j][0])))*(U1_xi);

	w_zeta = (1.0 / U[i][j][0])*(U4_zeta)-(U[i][j][3] / ((U[i][j][0] * U[i][j][0])))*(U1_zeta);
	w_eta = (1.0 / U[i][j][0])*(U4_eta)-(U[i][j][3] / ((U[i][j][0] * U[i][j][0])))*(U1_eta);
	w_xi = (1.0 / U[i][j][0])*(U4_xi)-(U[i][j][3] / ((U[i][j][0] * U[i][j][0])))*(U1_xi);

	t_zeta = (1.0 / U[i][j][0])*(U5_zeta)-(U[i][j][4] / ((U[i][j][0] * U[i][j][0])))*(U1_zeta)-(U[i][j][1] / U[i][j][0])*u_zeta - (U[i][j][2] / U[i][j][0])*v_zeta - (U[i][j][3] / U[i][j][0])*w_zeta;
	t_eta = (1.0 / U[i][j][0])*(U5_eta)-(U[i][j][4] / ((U[i][j][0] * U[i][j][0])))*(U1_eta)-(U[i][j][1] / U[i][j][0])*u_eta - (U[i][j][2] / U[i][j][0])*v_eta - (U[i][j][3] / U[i][j][0])*w_eta;
	t_xi = (1.0 / U[i][j][0])*(U5_xi)-(U[i][j][4] / ((U[i][j][0] * U[i][j][0])))*(U1_xi)-(U[i][j][1] / U[i][j][0])*u_xi - (U[i][j][2] / U[i][j][0])*v_xi - (U[i][j][3] / U[i][j][0])*w_xi;

	d_temp = (1.4*Mach*Mach*0.4*((U[i][j][4] / U[i][j][0]) - 0.5*(((U[i][j][1] * U[i][j][1]) / (U[i][j][0] * U[i][j][0])) + ((U[i][j][2] * U[i][j][2]) / (U[i][j][0] * U[i][j][0])) + ((U[i][j][3] * U[i][j][3]) / (U[i][j][0] * U[i][j][0])))));
	d_temp1 = sqrt(d_temp)*sqrt(d_temp)*sqrt(d_temp);
	/**************************************SUTHERLAND LAW***************************************************************************/
	mu[j][i] = (d_temp1)*((1.0 + (110.4 / Free_t)) / (d_temp + (110.4 / Free_t)));
	/********************************************************************************************************************************/
	/*****************************************HEAT FLUX******************************************************************************/
	qz[i] = (-1.0)*(mu[j][i] * 1.4 / (Pr*Reyl))*(t_zeta*metric[i].zeta_x + t_eta * metric[i].eta_x + t_xi * metric[i].xi_x);
	qe[i] = (-1.0)*(mu[j][i] * 1.4 / (Pr*Reyl))*(t_zeta*metric[i].zeta_y + t_eta * metric[i].eta_y + t_xi * metric[i].xi_y);
	qx[i] = (-1.0)*(mu[j][i] * 1.4 / (Pr*Reyl))*(t_zeta*metric[i].zeta_z + t_eta * metric[i].eta_z + t_xi * metric[i].xi_z);
	/*********************************************************************************************************************************/
	/***************************************************viscous shear stess terms******************************************************/
	tauzz[i] = (mu[j][i] / Reyl)*((4.0 / 3.0)*(u_zeta*metric[i].zeta_x + u_eta * metric[i].eta_x + u_xi * metric[i].xi_x) - (2.0 / 3.0)*(v_zeta*metric[i].zeta_y + v_eta * metric[i].eta_y + v_xi * metric[i].xi_y) - (2.0 / 3.0)*(metric[i].zeta_z*w_zeta + metric[i].eta_z*w_eta + metric[i].xi_z*w_xi));

	tauee[i] = (mu[j][i] / Reyl)*((4.0 / 3.0)*(v_zeta*metric[i].zeta_y + v_eta * metric[i].eta_y + v_xi * metric[i].xi_y) - (2.0 / 3.0)*(u_zeta*metric[i].zeta_x + u_eta * metric[i].eta_x + u_xi * metric[i].xi_x) - (2.0 / 3.0)*(metric[i].zeta_z*w_zeta + metric[i].eta_z*w_eta + metric[i].xi_z*w_xi));

	tauxx[i] = (mu[j][i] / Reyl)*((4.0 / 3.0)*(w_zeta*metric[i].zeta_z + w_eta * metric[i].eta_z + w_xi * metric[i].xi_z) - (2.0 / 3.0)*(u_zeta*metric[i].zeta_x + u_eta * metric[i].eta_x + u_xi * metric[i].xi_x) - (2.0 / 3.0)*(metric[i].zeta_y*v_zeta + metric[i].eta_y*v_eta + metric[i].xi_y*v_xi));

	tauze[i] = (mu[j][i] / Reyl)*(u_zeta*metric[i].zeta_y + u_eta * metric[i].eta_y + u_xi * metric[i].xi_y + v_zeta * metric[i].zeta_x + v_eta * metric[i].eta_x + v_xi * metric[i].xi_x);

	tauzx[i] = (mu[j][i] / Reyl)*(metric[i].zeta_z*u_zeta + metric[i].eta_z*u_eta + metric[i].xi_z*u_xi + metric[i].zeta_x*w_zeta + metric[i].eta_x*w_eta + metric[i].xi_x*w_xi);

	tauex[i] = (mu[j][i] / Reyl)*(metric[i].zeta_z*v_zeta + metric[i].eta_z*v_eta + metric[i].xi_z*v_xi + metric[i].zeta_y*w_zeta + metric[i].eta_y*w_eta + metric[i].xi_y*w_xi);

	/*******************************************DUCROS SENSOR*********************************************************************/
	divergence_V = u_zeta * metric[i].zeta_x + u_eta * metric[i].eta_x + u_xi * metric[i].xi_x + v_zeta * metric[i].zeta_y + v_eta * metric[i].eta_y + v_xi * metric[i].xi_y + w_zeta * metric[i].zeta_z + w_eta * metric[i].eta_z + w_xi * metric[i].xi_z;
	cross_V = pow(v_zeta*metric[i].zeta_z + v_eta * metric[i].eta_z + v_xi * metric[i].xi_z - (w_zeta*metric[i].zeta_y + w_eta * metric[i].eta_y + w_xi * metric[i].xi_y), 2.0) + \
		pow(u_zeta*metric[i].zeta_z + u_eta * metric[i].eta_z + u_xi * metric[i].xi_z - (w_zeta*metric[i].zeta_x + w_eta * metric[i].eta_x + w_xi * metric[i].xi_x), 2.0) + \
		pow(u_zeta*metric[i].zeta_y + u_eta * metric[i].eta_y + u_xi * metric[i].xi_y - (v_zeta*metric[i].zeta_x + v_eta * metric[i].eta_x + v_xi * metric[i].xi_x), 2.0);

	DUCROS[i] = (divergence_V*divergence_V) / ((divergence_V*divergence_V) + (cross_V)+10e-30);




	/****************************************************************************************************************************/
	/****************************************LARGE EDDY SIMULATION TERMS*********************************************************/
	/***********************************MIXED TIME SCALE EDDY VISCOSITY TERMS****************************************************/
	/*******************************************SGS MOMENTUM TERMS***************************************************************/
	C_MTS = 0.03;
	C_T = 10.0;

	KS_x = pow(((U[i][j][1] / U[i][j][0]) - ((0.25*U[node[i].n_n[3]][j][1] + 0.5*U[i][j][1] + 0.25*U[node[i].n_n[1]][j][1]) / (0.25*U[node[i].n_n[3]][j][0] + 0.5*U[i][j][0] + 0.25*U[node[i].n_n[1]][j][0]))), 2.0);
	KS_y = pow(((U[i][j][2] / U[i][j][0]) - ((0.25*U[node[i].n_n[2]][j][2] + 0.5*U[i][j][2] + 0.25*U[node[i].n_n[0]][j][2]) / (0.25*U[node[i].n_n[2]][j][0] + 0.5*U[i][j][0] + 0.25*U[node[i].n_n[0]][j][0]))), 2.0);
	KS_z = pow(((U[i][j][3] / U[i][j][0]) - ((0.25*U[node[i].n_n[5]][j][3] + 0.5*U[i][j][3] + 0.25*U[node[i].n_n[4]][j][3]) / (0.25*U[node[i].n_n[5]][j][0] + 0.5*U[i][j][0] + 0.25*U[node[i].n_n[4]][j][0]))), 2.0);

	stress_xy_star = (0.5*(u_zeta*metric[i].zeta_y + u_eta * metric[i].eta_y + u_xi * metric[i].xi_y + v_zeta * metric[i].zeta_x + v_eta * metric[i].eta_x + v_xi * metric[i].xi_x));
	stress_yz_star = (0.5*(v_zeta*metric[i].zeta_z + v_eta * metric[i].eta_z + v_xi * metric[i].xi_z + w_zeta * metric[i].zeta_y + w_eta * metric[i].eta_y + w_xi * metric[i].xi_y));
	stress_xz_star = (0.5*(w_zeta*metric[i].zeta_x + w_eta * metric[i].eta_x + w_xi * metric[i].xi_x + u_zeta * metric[i].zeta_z + u_eta * metric[i].eta_z + u_xi * metric[i].xi_z));

	stress_xy = (0.5*(u_zeta*metric[i].zeta_y + u_eta * metric[i].eta_y + u_xi * metric[i].xi_y + v_zeta * metric[i].zeta_x + v_eta * metric[i].eta_x + v_xi * metric[i].xi_x));
	stress_yz = (0.5*(v_zeta*metric[i].zeta_z + v_eta * metric[i].eta_z + v_xi * metric[i].xi_z + w_zeta * metric[i].zeta_y + w_eta * metric[i].eta_y + w_xi * metric[i].xi_y));
	stress_xz = (0.5*(w_zeta*metric[i].zeta_x + w_eta * metric[i].eta_x + w_xi * metric[i].xi_x + u_zeta * metric[i].zeta_z + u_eta * metric[i].eta_z + u_xi * metric[i].xi_z));

	stress_xx = (u_zeta*metric[i].zeta_x + u_eta * metric[i].eta_x + u_xi * metric[i].xi_x);
	stress_yy = (v_zeta*metric[i].zeta_y + v_eta * metric[i].eta_y + v_xi * metric[i].xi_y);
	stress_zz = (w_zeta*metric[i].zeta_z + w_eta * metric[i].eta_z + w_xi * metric[i].xi_z);

	mod_stress = sqrt(2.0*((stress_xy*stress_xy + stress_yz * stress_yz + stress_xz * stress_xz) + (2.0 / 3.0)*(stress_xx*stress_xx + stress_yy * stress_yy + stress_zz * stress_zz)));

	DEL_ZETA = (1.0 / cbrt(det[i]));
	DEL_ETA = (1.0 / cbrt(det[i]));
	DEL_XI = (1.0 / cbrt(det[i]));

	/**********Be carefull with Rx, Ry, Rz******************************/
	/************THE RX Ry and Ry are different from paper to avoid zero divide*************/
	/************STUDY PAPER BEFORE TOUCHING Eddy viscosity ******************************/

	if (mod_stress != 0.0)
	{
		Rx = (DEL_ZETA*(mod_stress)) / sqrt(KS_x);
		Ry = (DEL_ZETA*(mod_stress)) / sqrt(KS_y);
		Rz = (DEL_ZETA*(mod_stress)) / sqrt(KS_z);

		EDDY_x = (C_MTS / (1.0 + (1.0*Rx / (C_T))))*DEL_ZETA*sqrt(KS_x);
		EDDY_y = (C_MTS / (1.0 + (1.0*Ry / (C_T))))*DEL_ZETA*sqrt(KS_y);
		EDDY_z = (C_MTS / (1.0 + (1.0*Rz / (C_T))))*DEL_ZETA*sqrt(KS_z);
	}
	if (mod_stress == 0.0)
	{
		EDDY_x = (C_MTS)*DEL_ZETA*sqrt(KS_x);
		EDDY_y = (C_MTS)*DEL_ZETA*sqrt(KS_y);
		EDDY_z = (C_MTS)*DEL_ZETA*sqrt(KS_z);
	}
	TAU_SGS_XY[i] = -2.0*U[i][j][0] * EDDY_x*(stress_xy_star);
	TAU_SGS_YZ[i] = -2.0*U[i][j][0] * EDDY_y*(stress_yz_star);
	TAU_SGS_XZ[i] = -2.0*U[i][j][0] * EDDY_z*(stress_xz_star);

	TAU_SGS_XX[i] = -2.0*U[i][j][0] * EDDY_x*stress_xx;
	TAU_SGS_YY[i] = -2.0*U[i][j][0] * EDDY_y*stress_yy;
	TAU_SGS_ZZ[i] = -2.0*U[i][j][0] * EDDY_z*stress_zz;
	/**********************************SGS ENERGY FLUX TERMS*******************************************************************/

	Q_x = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_x)*(t_zeta*metric[i].zeta_x + t_eta * metric[i].eta_x + t_xi * metric[i].xi_x);
	Q_y = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_y)*(t_zeta*metric[i].zeta_y + t_eta * metric[i].eta_y + t_xi * metric[i].xi_y);
	Q_z = (1.0 / (1.4*0.4*Mach*Mach))*(-1.0)*(EDDY_z)*(t_zeta*metric[i].zeta_z + t_eta * metric[i].eta_z + t_xi * metric[i].xi_z);

	DEL_X = TAU_SGS_XX[i] * u[j][i] + TAU_SGS_XY[i] * u[j][i] + TAU_SGS_XZ[i] * u[j][i];
	DEL_Y = TAU_SGS_XY[i] * v[j][i] + TAU_SGS_YY[i] * v[j][i] + TAU_SGS_YZ[i] * v[j][i];
	DEL_Z = TAU_SGS_XZ[i] * w[j][i] + TAU_SGS_YZ[i] * w[j][i] + TAU_SGS_ZZ[i] * w[j][i];

	H_SGS_X[i] = Q_x + DEL_X;
	H_SGS_Y[i] = Q_y + DEL_Y;
	H_SGS_Z[i] = Q_z + DEL_Z;

	D_SGS_X[i] = 0.4*Q_x;
	D_SGS_Y[i] = 0.4*Q_y;
	D_SGS_Z[i] = 0.4*Q_z;

}


