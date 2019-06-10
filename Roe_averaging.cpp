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

void roe_average()
{

	double sqr_b1, sqr_i, sqr_b3, sqr_c0, sqr_c2, sqr_d4, sqr_d5;
	sqr_i = sqrt(rho[j][i]);
	sqr_b1 = sqrt(rho[j][node[i].n_n[1]]);
	sqr_b3 = sqrt(rho[j][node[i].n_n[3]]);
	sqr_c0 = sqrt(rho[j][node[i].n_n[0]]);
	sqr_c2 = sqrt(rho[j][node[i].n_n[2]]);
	sqr_d4 = sqrt(rho[j][node[i].n_n[4]]);
	sqr_d5 = sqrt(rho[j][node[i].n_n[5]]);


	roe_rho_ip[i] = sqr_b1 * sqr_i;
	roe_u_ip[i] = (sqr_b1*u[j][node[i].n_n[1]] + sqr_i * u[j][i]) / (sqr_b1 + sqr_i);
	roe_v_ip[i] = (sqr_b1*v[j][node[i].n_n[1]] + sqr_i * v[j][i]) / (sqr_b1 + sqr_i);
	roe_w_ip[i] = (sqr_b1*w[j][node[i].n_n[1]] + sqr_i * w[j][i]) / (sqr_b1 + sqr_i);
	roe_h_ip[i] = (sqr_b1*(((p[j][node[i].n_n[1]] / 0.4) + 0.5*rho[j][node[i].n_n[1]] * (u[j][node[i].n_n[1]] * u[j][node[i].n_n[1]] + v[j][node[i].n_n[1]] * v[j][node[i].n_n[1]] + w[j][node[i].n_n[1]] * w[j][node[i].n_n[1]]) + p[j][node[i].n_n[1]]) / rho[j][node[i].n_n[1]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_b1 + sqr_i);

	roe_rho_im[i] = sqr_b3 * sqr_i;
	roe_u_im[i] = (sqr_b3*u[j][node[i].n_n[3]] + sqr_i * u[j][i]) / (sqr_b3 + sqr_i);
	roe_v_im[i] = (sqr_b3*v[j][node[i].n_n[3]] + sqr_i * v[j][i]) / (sqr_b3 + sqr_i);
	roe_w_im[i] = (sqr_b3*w[j][node[i].n_n[3]] + sqr_i * w[j][i]) / (sqr_b3 + sqr_i);
	roe_h_im[i] = (sqr_b3*(((p[j][node[i].n_n[3]] / 0.4) + 0.5*rho[j][node[i].n_n[3]] * (u[j][node[i].n_n[3]] * u[j][node[i].n_n[3]] + v[j][node[i].n_n[3]] * v[j][node[i].n_n[3]] + w[j][node[i].n_n[3]] * w[j][node[i].n_n[3]]) + p[j][node[i].n_n[3]]) / rho[j][node[i].n_n[3]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_b3 + sqr_i);

	roe_rho_jp[i] = sqr_c0 * sqr_i;
	roe_u_jp[i] = (sqr_c0*u[j][node[i].n_n[0]] + sqr_i * u[j][i]) / (sqr_c0 + sqr_i);
	roe_v_jp[i] = (sqr_c0*v[j][node[i].n_n[0]] + sqr_i * v[j][i]) / (sqr_c0 + sqr_i);
	roe_w_jp[i] = (sqr_c0*w[j][node[i].n_n[0]] + sqr_i * w[j][i]) / (sqr_c0 + sqr_i);
	roe_h_jp[i] = (sqr_c0*(((p[j][node[i].n_n[0]] / 0.4) + 0.5*rho[j][node[i].n_n[0]] * (u[j][node[i].n_n[0]] * u[j][node[i].n_n[0]] + v[j][node[i].n_n[0]] * v[j][node[i].n_n[0]] + w[j][node[i].n_n[0]] * w[j][node[i].n_n[0]]) + p[j][node[i].n_n[0]]) / rho[j][node[i].n_n[0]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_c0 + sqr_i);

	roe_rho_jm[i] = sqr_c2 * sqr_i;
	roe_u_jm[i] = (sqr_c2*u[j][node[i].n_n[2]] + sqr_i * u[j][i]) / (sqr_c2 + sqr_i);
	roe_v_jm[i] = (sqr_c2*v[j][node[i].n_n[2]] + sqr_i * v[j][i]) / (sqr_c2 + sqr_i);
	roe_w_jm[i] = (sqr_c2*w[j][node[i].n_n[2]] + sqr_i * w[j][i]) / (sqr_c2 + sqr_i);
	roe_h_jm[i] = (sqr_c2*(((p[j][node[i].n_n[2]] / 0.4) + 0.5*rho[j][node[i].n_n[2]] * (u[j][node[i].n_n[2]] * u[j][node[i].n_n[2]] + v[j][node[i].n_n[2]] * v[j][node[i].n_n[2]] + w[j][node[i].n_n[2]] * w[j][node[i].n_n[2]]) + p[j][node[i].n_n[2]]) / rho[j][node[i].n_n[2]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_c2 + sqr_i);

	roe_rho_kp[i] = sqr_d4 * sqr_i;
	roe_u_kp[i] = (sqr_d4*u[j][node[i].n_n[4]] + sqr_i * u[j][i]) / (sqr_d4 + sqr_i);
	roe_v_kp[i] = (sqr_d4*v[j][node[i].n_n[4]] + sqr_i * v[j][i]) / (sqr_d4 + sqr_i);
	roe_w_kp[i] = (sqr_d4*w[j][node[i].n_n[4]] + sqr_i * w[j][i]) / (sqr_d4 + sqr_i);
	roe_h_kp[i] = (sqr_d4*(((p[j][node[i].n_n[4]] / 0.4) + 0.5*rho[j][node[i].n_n[4]] * (u[j][node[i].n_n[4]] * u[j][node[i].n_n[4]] + v[j][node[i].n_n[4]] * v[j][node[i].n_n[4]] + w[j][node[i].n_n[4]] * w[j][node[i].n_n[4]]) + p[j][node[i].n_n[4]]) / rho[j][node[i].n_n[4]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_d4 + sqr_i);

	roe_rho_km[i] = sqr_d5 * sqr_i;
	roe_u_km[i] = (sqr_d5*u[j][node[i].n_n[5]] + sqr_i * u[j][i]) / (sqr_d5 + sqr_i);
	roe_v_km[i] = (sqr_d5*v[j][node[i].n_n[5]] + sqr_i * v[j][i]) / (sqr_d5 + sqr_i);
	roe_w_km[i] = (sqr_d5*w[j][node[i].n_n[5]] + sqr_i * w[j][i]) / (sqr_d5 + sqr_i);
	roe_h_km[i] = (sqr_d5*(((p[j][node[i].n_n[5]] / 0.4) + 0.5*rho[j][node[i].n_n[5]] * (u[j][node[i].n_n[5]] * u[j][node[i].n_n[5]] + v[j][node[i].n_n[5]] * v[j][node[i].n_n[5]] + w[j][node[i].n_n[5]] * w[j][node[i].n_n[5]]) + p[j][node[i].n_n[5]]) / rho[j][node[i].n_n[5]]) + \
		sqr_i*(((p[j][i] / 0.4) + 0.5*rho[j][i] * (u[j][i] * u[j][i] + v[j][i] * v[j][i] + w[j][i] * w[j][i]) + p[j][i]) / rho[j][i])) / (sqr_d5 + sqr_i);


	if (i <= sd_node)
	{
		diver[i][0].u = u[0][i];
		diver[i][0].v = v[0][i];
		diver[i][0].e = e[0][i];
	}
}


