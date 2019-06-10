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

void metric_term()
{
	int step, M, ci, repeat, k, gar;
	double cx[2000], cy[2000], cz[2000], d[2000], a[2000], b[2000], d_dash[2000], cx_dash[2000], cy_dash[2000], cz_dash[2000];
	double c1[2000], c2[2000], c3[2000], c1_dash[2000], c2_dash[2000], c3_dash[2000];
	double d1[2000], d2[2000], d3[2000], d1_dash[2000], d2_dash[2000], d3_dash[2000];

	for (h = 0; h<all_bou_node; h++)
	{
		i = all_boundary_nodes[h];

		repeat = 0;
		if ((node[i].n_n[1] == 0) && node[i].ID != 10)
		{
			if (node[i].n_n[1] == 0)
			{
				ci = 1;
				k = 3;
			}
			step = 0;
			gar = 0;
			while (i != 0)
			{
				if (step > 1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[3]].n_n[3]].x - (7.0 / 9.0)*node[node[i].n_n[3]].x + (7.0 / 9.0)*node[node[i].n_n[1]].x + (1.0 / 36.0)*node[node[node[i].n_n[1]].n_n[1]].x);
					cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[3]].n_n[3]].y - (7.0 / 9.0)*node[node[i].n_n[3]].y + (7.0 / 9.0)*node[node[i].n_n[1]].y + (1.0 / 36.0)*node[node[node[i].n_n[1]].n_n[1]].y);
					cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[3]].n_n[3]].z - (7.0 / 9.0)*node[node[i].n_n[3]].z + (7.0 / 9.0)*node[node[i].n_n[1]].z + (1.0 / 36.0)*node[node[node[i].n_n[1]].n_n[1]].z);
				}

				if (step == 0 && node[i].n_n[1] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					cx[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[3]].x - (5.0)*node[node[node[i].n_n[3]].n_n[3]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].x);
					cy[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[3]].y - (5.0)*node[node[node[i].n_n[3]].n_n[3]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].y);
					cz[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[3]].z - (5.0)*node[node[node[i].n_n[3]].n_n[3]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].z);
				}
				if (step == 1)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[1]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[3]].x + (7.0 / 33.0)*node[node[node[i].n_n[3]].n_n[3]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].x);
					cy[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[1]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[3]].y + (7.0 / 33.0)*node[node[node[i].n_n[3]].n_n[3]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].y);
					cz[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[1]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[3]].z + (7.0 / 33.0)*node[node[node[i].n_n[3]].n_n[3]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].z);
				}

				if (step > 1 && node[node[i].n_n[3]].n_n[3] == 0 && node[i].n_n[3] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[3]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[1]].x + (7.0 / 33.0)*node[node[node[i].n_n[1]].n_n[1]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[3]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[1]].y + (7.0 / 33.0)*node[node[node[i].n_n[1]].n_n[1]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[3]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[1]].z + (7.0 / 33.0)*node[node[node[i].n_n[1]].n_n[1]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].z);
				}
				if (step != 0 && node[i].n_n[3] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[1]].x - (5.0)*node[node[node[i].n_n[1]].n_n[1]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[1]].y - (5.0)*node[node[node[i].n_n[1]].n_n[1]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[1]].z - (5.0)*node[node[node[i].n_n[1]].n_n[1]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].z);
				}

				if (step >1 && (node[i].ID == 10 || node[node[node[i].n_n[1]].n_n[1]].ID == 10 || node[node[i].n_n[1]].ID == 10))
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					if (node[i].e[2] != 0 && (CD[node[i].e[2]].connect[1] == gar || CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1] == gar || node[i].ID == 0))
					{
						cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].x - (7.0 / 9.0)*node[CD[node[i].e[1]].connect[0]].x + (7.0 / 9.0)*node[CD[node[i].e[2]].connect[1]].x + (1.0 / 36.0)*node[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].x);
						cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].y - (7.0 / 9.0)*node[CD[node[i].e[1]].connect[0]].y + (7.0 / 9.0)*node[CD[node[i].e[2]].connect[1]].y + (1.0 / 36.0)*node[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].y);
						cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].z - (7.0 / 9.0)*node[CD[node[i].e[1]].connect[0]].z + (7.0 / 9.0)*node[CD[node[i].e[2]].connect[1]].z + (1.0 / 36.0)*node[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].z);
					}
					if (node[i].e[3] != 0 && (CD[node[i].e[3]].connect[5] == gar || CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5] == gar || node[i].ID == 0))
					{
						cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].x - (7.0 / 9.0)*node[CD[node[i].e[0]].connect[4]].x + (7.0 / 9.0)*node[CD[node[i].e[3]].connect[5]].x + (1.0 / 36.0)*node[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].x);
						cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].y - (7.0 / 9.0)*node[CD[node[i].e[0]].connect[4]].y + (7.0 / 9.0)*node[CD[node[i].e[3]].connect[5]].y + (1.0 / 36.0)*node[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].y);
						cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].z - (7.0 / 9.0)*node[CD[node[i].e[0]].connect[4]].z + (7.0 / 9.0)*node[CD[node[i].e[3]].connect[5]].z + (1.0 / 36.0)*node[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].z);
					}


					if (node[i].e[7] != 0 && (CD[node[i].e[7]].connect[6] == gar || CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6] == gar || node[i].ID == 0))
					{
						cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].x - (7.0 / 9.0)*node[CD[node[i].e[4]].connect[7]].x + (7.0 / 9.0)*node[CD[node[i].e[7]].connect[6]].x + (1.0 / 36.0)*node[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].x);
						cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].y - (7.0 / 9.0)*node[CD[node[i].e[4]].connect[7]].y + (7.0 / 9.0)*node[CD[node[i].e[7]].connect[6]].y + (1.0 / 36.0)*node[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].y);
						cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].z - (7.0 / 9.0)*node[CD[node[i].e[4]].connect[7]].z + (7.0 / 9.0)*node[CD[node[i].e[7]].connect[6]].z + (1.0 / 36.0)*node[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].z);
					}
					if (node[i].e[6] != 0 && (CD[node[i].e[6]].connect[2] == gar || CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2] == gar || node[i].ID == 0))
					{
						cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].x - (7.0 / 9.0)*node[CD[node[i].e[5]].connect[3]].x + (7.0 / 9.0)*node[CD[node[i].e[6]].connect[2]].x + (1.0 / 36.0)*node[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].x);
						cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].y - (7.0 / 9.0)*node[CD[node[i].e[5]].connect[3]].y + (7.0 / 9.0)*node[CD[node[i].e[6]].connect[2]].y + (1.0 / 36.0)*node[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].y);
						cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].z - (7.0 / 9.0)*node[CD[node[i].e[5]].connect[3]].z + (7.0 / 9.0)*node[CD[node[i].e[6]].connect[2]].z + (1.0 / 36.0)*node[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].z);
					}
				}

				if (node[node[i].n_n[k]].ID == 10)
				{
					gar = i;
				}
				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;
			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);

			cx_dash[M - 1] = cx[M - 1] - ((cx[M] * b[M - 1]) / d[M]);
			cx_dash[M - 2] = cx[M - 2] - ((cx_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cy_dash[M - 1] = cy[M - 1] - ((cy[M] * b[M - 1]) / d[M]);
			cy_dash[M - 2] = cy[M - 2] - ((cy_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cz_dash[M - 1] = cz[M - 1] - ((cz[M] * b[M - 1]) / d[M]);
			cz_dash[M - 2] = cz[M - 2] - ((cz_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				cx_dash[i] = cx[i] - ((cx_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cy_dash[i] = cy[i] - ((cy_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cz_dash[i] = cz[i] - ((cz_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jacobian[i].x_zeta = (cx_dash[step] - a[step] * jacobian[node[i].n_n[ci]].x_zeta) / d_dash[step];
					jacobian[i].y_zeta = (cy_dash[step] - a[step] * jacobian[node[i].n_n[ci]].y_zeta) / d_dash[step];
					jacobian[i].z_zeta = (cz_dash[step] - a[step] * jacobian[node[i].n_n[ci]].z_zeta) / d_dash[step];

					if (node[i].ID == 10)
					{
						if (node[i].e[2] != 0 && CD[node[i].e[2]].connect[1] == gar)
						{
							jacobian[i].x_zeta = (cx_dash[step] - a[step] * jacobian[CD[node[i].e[2]].connect[1]].x_zeta) / d_dash[step];
							jacobian[i].y_zeta = (cy_dash[step] - a[step] * jacobian[CD[node[i].e[2]].connect[1]].y_zeta) / d_dash[step];
							jacobian[i].z_zeta = (cz_dash[step] - a[step] * jacobian[CD[node[i].e[2]].connect[1]].z_zeta) / d_dash[step];
						}
						if (node[i].e[3] != 0 && CD[node[i].e[3]].connect[5] == gar)
						{
							jacobian[i].x_zeta = (cx_dash[step] - a[step] * jacobian[CD[node[i].e[3]].connect[5]].x_zeta) / d_dash[step];
							jacobian[i].y_zeta = (cy_dash[step] - a[step] * jacobian[CD[node[i].e[3]].connect[5]].y_zeta) / d_dash[step];
							jacobian[i].z_zeta = (cz_dash[step] - a[step] * jacobian[CD[node[i].e[3]].connect[5]].z_zeta) / d_dash[step];
						}
						if (node[i].e[7] != 0 && CD[node[i].e[7]].connect[6] == gar)
						{
							jacobian[i].x_zeta = (cx_dash[step] - a[step] * jacobian[CD[node[i].e[7]].connect[6]].x_zeta) / d_dash[step];
							jacobian[i].y_zeta = (cy_dash[step] - a[step] * jacobian[CD[node[i].e[7]].connect[6]].y_zeta) / d_dash[step];
							jacobian[i].z_zeta = (cz_dash[step] - a[step] * jacobian[CD[node[i].e[7]].connect[6]].z_zeta) / d_dash[step];
						}
						if (node[i].e[6] != 0 && CD[node[i].e[6]].connect[2] == gar)
						{
							jacobian[i].x_zeta = (cx_dash[step] - a[step] * jacobian[CD[node[i].e[6]].connect[2]].x_zeta) / d_dash[step];
							jacobian[i].y_zeta = (cy_dash[step] - a[step] * jacobian[CD[node[i].e[6]].connect[2]].y_zeta) / d_dash[step];
							jacobian[i].z_zeta = (cz_dash[step] - a[step] * jacobian[CD[node[i].e[6]].connect[2]].z_zeta) / d_dash[step];
						}
					}
				}
				if (step == 0)
				{
					jacobian[i].x_zeta = cx_dash[step] / d_dash[0];
					jacobian[i].y_zeta = cy_dash[step] / d_dash[0];
					jacobian[i].z_zeta = cz_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jacobian[i].x_zeta = (cx[step] - a[M] * jacobian[node[i].n_n[ci]].x_zeta) / d[M];
					jacobian[i].y_zeta = (cy[step] - a[M] * jacobian[node[i].n_n[ci]].y_zeta) / d[M];
					jacobian[i].z_zeta = (cz[step] - a[M] * jacobian[node[i].n_n[ci]].z_zeta) / d[M];
				}

				i = node[i].n_n[k];

				step++;
			}
		}
		/******************************************************************y-derivative*******************************************************************/
		i = all_boundary_nodes[h];
		if ((node[i].n_n[0] == 0) && node[i].ID != 10)
		{
			if (node[i].n_n[0] == 0)
			{
				ci = 0;
				k = 2;
			}
			step = 0;
			while (i != 0)
			{
				if (step>1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[2]].n_n[2]].x - (7.0 / 9.0)*node[node[i].n_n[2]].x + (7.0 / 9.0)*node[node[i].n_n[0]].x + (1.0 / 36.0)*node[node[node[i].n_n[0]].n_n[0]].x);
					cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[2]].n_n[2]].y - (7.0 / 9.0)*node[node[i].n_n[2]].y + (7.0 / 9.0)*node[node[i].n_n[0]].y + (1.0 / 36.0)*node[node[node[i].n_n[0]].n_n[0]].y);
					cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[2]].n_n[2]].z - (7.0 / 9.0)*node[node[i].n_n[2]].z + (7.0 / 9.0)*node[node[i].n_n[0]].z + (1.0 / 36.0)*node[node[node[i].n_n[0]].n_n[0]].z);
				}

				if (step == 0 && node[i].n_n[0] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					cx[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[2]].x - (5.0)*node[node[node[i].n_n[2]].n_n[2]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].x);
					cy[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[2]].y - (5.0)*node[node[node[i].n_n[2]].n_n[2]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].y);
					cz[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[2]].z - (5.0)*node[node[node[i].n_n[2]].n_n[2]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].z);
				}

				if (step == 1)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[0]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[2]].x + (7.0 / 33.0)*node[node[node[i].n_n[2]].n_n[2]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].x);
					cy[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[0]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[2]].y + (7.0 / 33.0)*node[node[node[i].n_n[2]].n_n[2]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].y);
					cz[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[0]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[2]].z + (7.0 / 33.0)*node[node[node[i].n_n[2]].n_n[2]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].z);
				}

				if (step > 1 && node[node[i].n_n[2]].n_n[2] == 0 && node[i].n_n[2] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[2]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[0]].x + (7.0 / 33.0)*node[node[node[i].n_n[0]].n_n[0]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[2]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[0]].y + (7.0 / 33.0)*node[node[node[i].n_n[0]].n_n[0]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[2]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[0]].z + (7.0 / 33.0)*node[node[node[i].n_n[0]].n_n[0]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].z);
				}

				if (step != 0 && node[i].n_n[2] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[0]].x - (5.0)*node[node[node[i].n_n[0]].n_n[0]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[0]].y - (5.0)*node[node[node[i].n_n[0]].n_n[0]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[0]].z - (5.0)*node[node[node[i].n_n[0]].n_n[0]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].z);
				}

				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;

			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);

			cx_dash[M - 1] = cx[M - 1] - ((cx[M] * b[M - 1]) / d[M]);
			cx_dash[M - 2] = cx[M - 2] - ((cx_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cy_dash[M - 1] = cy[M - 1] - ((cy[M] * b[M - 1]) / d[M]);
			cy_dash[M - 2] = cy[M - 2] - ((cy_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cz_dash[M - 1] = cz[M - 1] - ((cz[M] * b[M - 1]) / d[M]);
			cz_dash[M - 2] = cz[M - 2] - ((cz_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				cx_dash[i] = cx[i] - ((cx_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cy_dash[i] = cy[i] - ((cy_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cz_dash[i] = cz[i] - ((cz_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jacobian[i].x_eta = (cx_dash[step] - a[step] * jacobian[node[i].n_n[ci]].x_eta) / d_dash[step];
					jacobian[i].y_eta = (cy_dash[step] - a[step] * jacobian[node[i].n_n[ci]].y_eta) / d_dash[step];
					jacobian[i].z_eta = (cz_dash[step] - a[step] * jacobian[node[i].n_n[ci]].z_eta) / d_dash[step];
				}
				if (step == 0)
				{
					jacobian[i].x_eta = cx_dash[step] / d_dash[0];
					jacobian[i].y_eta = cy_dash[step] / d_dash[0];
					jacobian[i].z_eta = cz_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jacobian[i].x_eta = (cx[step] - a[M] * jacobian[node[i].n_n[ci]].x_eta) / d[M];
					jacobian[i].y_eta = (cy[step] - a[M] * jacobian[node[i].n_n[ci]].y_eta) / d[M];
					jacobian[i].z_eta = (cz[step] - a[M] * jacobian[node[i].n_n[ci]].z_eta) / d[M];
				}
				i = node[i].n_n[k];
				step++;
			}
		}
		/******************************************************************z-derivative*******************************************************************/
		i = all_boundary_nodes[h];
		if ((node[i].n_n[4] == 0) && (node[i].loc == 30 || node[i].loc == 31 || node[i].loc <= 27 || node[i].loc == 28))
		{
			if (node[i].n_n[4] == 0)
			{
				ci = 4;
				k = 5;
			}
			step = 0;
			while (i != 0)
			{
				if (step>1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					cx[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[5]].n_n[5]].x - (7.0 / 9.0)*node[node[i].n_n[5]].x + (7.0 / 9.0)*node[node[i].n_n[4]].x + (1.0 / 36.0)*node[node[node[i].n_n[4]].n_n[4]].x);
					cy[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[5]].n_n[5]].y - (7.0 / 9.0)*node[node[i].n_n[5]].y + (7.0 / 9.0)*node[node[i].n_n[4]].y + (1.0 / 36.0)*node[node[node[i].n_n[4]].n_n[4]].y);
					cz[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*node[node[node[i].n_n[5]].n_n[5]].z - (7.0 / 9.0)*node[node[i].n_n[5]].z + (7.0 / 9.0)*node[node[i].n_n[4]].z + (1.0 / 36.0)*node[node[node[i].n_n[4]].n_n[4]].z);
				}

				if (step == 0 && node[i].n_n[4] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					cx[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[5]].x - (5.0)*node[node[node[i].n_n[5]].n_n[5]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].x);
					cy[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[5]].y - (5.0)*node[node[node[i].n_n[5]].n_n[5]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].y);
					cz[step] = (1.0 / del_zeta)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[5]].z - (5.0)*node[node[node[i].n_n[5]].n_n[5]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].z);
				}

				if (step == 1)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[4]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[5]].x + (7.0 / 33.0)*node[node[node[i].n_n[5]].n_n[5]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].x);
					cy[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[4]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[5]].y + (7.0 / 33.0)*node[node[node[i].n_n[5]].n_n[5]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].y);
					cz[step] = (1.0 / del_zeta)*((20.0 / 33.0)*node[node[i].n_n[4]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[5]].z + (7.0 / 33.0)*node[node[node[i].n_n[5]].n_n[5]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].z);
				}

				if (step > 1 && node[node[i].n_n[5]].n_n[5] == 0 && node[i].n_n[5] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[5]].x + (35.0 / 132.0)*node[i].x - (34.0 / 33.0)*node[node[i].n_n[4]].x + (7.0 / 33.0)*node[node[node[i].n_n[4]].n_n[4]].x - (2.0 / 33.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].x + (1.0 / 132.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[5]].y + (35.0 / 132.0)*node[i].y - (34.0 / 33.0)*node[node[i].n_n[4]].y + (7.0 / 33.0)*node[node[node[i].n_n[4]].n_n[4]].y - (2.0 / 33.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].y + (1.0 / 132.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*node[node[i].n_n[5]].z + (35.0 / 132.0)*node[i].z - (34.0 / 33.0)*node[node[i].n_n[4]].z + (7.0 / 33.0)*node[node[node[i].n_n[4]].n_n[4]].z - (2.0 / 33.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].z + (1.0 / 132.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].z);
				}

				if (step != 0 && node[i].n_n[5] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					cx[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].x + (5.0 / 12.0)*node[node[i].n_n[4]].x - (5.0)*node[node[node[i].n_n[4]].n_n[4]].x + (5.0 / 3.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].x - (5.0 / 12.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].x + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].x);
					cy[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].y + (5.0 / 12.0)*node[node[i].n_n[4]].y - (5.0)*node[node[node[i].n_n[4]].n_n[4]].y + (5.0 / 3.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].y - (5.0 / 12.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].y + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].y);
					cz[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*node[i].z + (5.0 / 12.0)*node[node[i].n_n[4]].z - (5.0)*node[node[node[i].n_n[4]].n_n[4]].z + (5.0 / 3.0)*node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].z - (5.0 / 12.0)*node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].z + (1.0 / 20.0)*node[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].z);
				}

				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;

			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);

			cx_dash[M - 1] = cx[M - 1] - ((cx[M] * b[M - 1]) / d[M]);
			cx_dash[M - 2] = cx[M - 2] - ((cx_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cy_dash[M - 1] = cy[M - 1] - ((cy[M] * b[M - 1]) / d[M]);
			cy_dash[M - 2] = cy[M - 2] - ((cy_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			cz_dash[M - 1] = cz[M - 1] - ((cz[M] * b[M - 1]) / d[M]);
			cz_dash[M - 2] = cz[M - 2] - ((cz_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				cx_dash[i] = cx[i] - ((cx_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cy_dash[i] = cy[i] - ((cy_dash[i + 1] * b[i]) / d_dash[i + 1]);
				cz_dash[i] = cz[i] - ((cz_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jacobian[i].x_xi = (cx_dash[step] - a[step] * jacobian[node[i].n_n[ci]].x_xi) / d_dash[step];
					jacobian[i].y_xi = (cy_dash[step] - a[step] * jacobian[node[i].n_n[ci]].y_xi) / d_dash[step];
					jacobian[i].z_xi = (cz_dash[step] - a[step] * jacobian[node[i].n_n[ci]].z_xi) / d_dash[step];
				}
				if (step == 0)
				{
					jacobian[i].x_xi = cx_dash[step] / d_dash[0];
					jacobian[i].y_xi = cy_dash[step] / d_dash[0];
					jacobian[i].z_xi = cz_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jacobian[i].x_xi = (cx[step] - a[M] * jacobian[node[i].n_n[ci]].x_xi) / d[M];
					jacobian[i].y_xi = (cy[step] - a[M] * jacobian[node[i].n_n[ci]].y_xi) / d[M];
					jacobian[i].z_xi = (cz[step] - a[M] * jacobian[node[i].n_n[ci]].z_xi) / d[M];
				}

				i = node[i].n_n[k];
				step++;
			}
		}
	}

	for (i = 1; i <= NUMNP; i++)
	{
		det[i] = (((jacobian[i].x_zeta)*((jacobian[i].y_eta*jacobian[i].z_xi) - (jacobian[i].y_xi*jacobian[i].z_eta))) + ((jacobian[i].x_eta)*((jacobian[i].y_xi*jacobian[i].z_zeta) - (jacobian[i].y_zeta*jacobian[i].z_xi))) + ((jacobian[i].x_xi)*((jacobian[i].y_zeta*jacobian[i].z_eta) - (jacobian[i].y_eta*jacobian[i].z_zeta))));

		jac[i].zeta_1 = jacobian[i].x_zeta*node[i].y;
		jac[i].zeta_2 = jacobian[i].y_zeta*node[i].z;
		jac[i].zeta_3 = jacobian[i].z_zeta*node[i].x;

		jac[i].eta_1 = jacobian[i].x_eta*node[i].y;
		jac[i].eta_2 = jacobian[i].y_eta*node[i].z;
		jac[i].eta_3 = jacobian[i].z_eta*node[i].x;

		jac[i].xi_1 = jacobian[i].x_xi*node[i].y;
		jac[i].xi_2 = jacobian[i].y_xi*node[i].z;
		jac[i].xi_3 = jacobian[i].z_xi*node[i].x;
	}


	/*******************************************************************************************************************************/
	/***************************************THOMAS AND LOMBARD(metric derivatives in conserved form)*******************************/
	for (h = 0; h<all_bou_node; h++)
	{
		i = all_boundary_nodes[h];
		repeat = 0;
		if ((node[i].n_n[1] == 0) && node[i].ID != 10)
		{
			if (node[i].n_n[1] == 0)
			{
				ci = 1;
				k = 3;
			}
			step = 0;
			gar = 0;
			while (i != 0)
			{
				if (step>1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_1 - (7.0 / 9.0)*jac[node[i].n_n[3]].xi_1 + (7.0 / 9.0)*jac[node[i].n_n[1]].xi_1 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_1);
					c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_2 - (7.0 / 9.0)*jac[node[i].n_n[3]].xi_2 + (7.0 / 9.0)*jac[node[i].n_n[1]].xi_2 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_2);
					c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_3 - (7.0 / 9.0)*jac[node[i].n_n[3]].xi_3 + (7.0 / 9.0)*jac[node[i].n_n[1]].xi_3 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_3);
					d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_1 - (7.0 / 9.0)*jac[node[i].n_n[3]].eta_1 + (7.0 / 9.0)*jac[node[i].n_n[1]].eta_1 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_1);
					d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_2 - (7.0 / 9.0)*jac[node[i].n_n[3]].eta_2 + (7.0 / 9.0)*jac[node[i].n_n[1]].eta_2 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_2);
					d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_3 - (7.0 / 9.0)*jac[node[i].n_n[3]].eta_3 + (7.0 / 9.0)*jac[node[i].n_n[1]].eta_3 + (1.0 / 36.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_3);
				}

				if (step == 0 && node[i].n_n[1] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					c1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_1 + (5.0 / 12.0)*jac[node[i].n_n[3]].xi_1 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_1);
					c2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_2 + (5.0 / 12.0)*jac[node[i].n_n[3]].xi_2 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_2);
					c3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_3 + (5.0 / 12.0)*jac[node[i].n_n[3]].xi_3 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_3);
					d1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_1 + (5.0 / 12.0)*jac[node[i].n_n[3]].eta_1 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_1);
					d2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_2 + (5.0 / 12.0)*jac[node[i].n_n[3]].eta_2 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_2);
					d3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_3 + (5.0 / 12.0)*jac[node[i].n_n[3]].eta_3 - (5.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_3);

				}

				if (step == 1)
				{

					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].xi_1 + (35.0 / 132.0)*jac[i].xi_1 - (34.0 / 33.0)*jac[node[i].n_n[3]].xi_1 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_1);
					c2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].xi_2 + (35.0 / 132.0)*jac[i].xi_2 - (34.0 / 33.0)*jac[node[i].n_n[3]].xi_2 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_2);
					c3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].xi_3 + (35.0 / 132.0)*jac[i].xi_3 - (34.0 / 33.0)*jac[node[i].n_n[3]].xi_3 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].xi_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].xi_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].xi_3);
					d1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].eta_1 + (35.0 / 132.0)*jac[i].eta_1 - (34.0 / 33.0)*jac[node[i].n_n[3]].eta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_1);
					d2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].eta_2 + (35.0 / 132.0)*jac[i].eta_2 - (34.0 / 33.0)*jac[node[i].n_n[3]].eta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_2);
					d3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[1]].eta_3 + (35.0 / 132.0)*jac[i].eta_3 - (34.0 / 33.0)*jac[node[i].n_n[3]].eta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[3]].n_n[3]].eta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].eta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[3]].n_n[3]].n_n[3]].n_n[3]].eta_3);
				}

				if (step > 1 && node[node[i].n_n[3]].n_n[3] == 0 && node[i].n_n[3] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].xi_1 + (35.0 / 132.0)*jac[i].xi_1 - (34.0 / 33.0)*jac[node[i].n_n[1]].xi_1 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].xi_2 + (35.0 / 132.0)*jac[i].xi_2 - (34.0 / 33.0)*jac[node[i].n_n[1]].xi_2 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].xi_3 + (35.0 / 132.0)*jac[i].xi_3 - (34.0 / 33.0)*jac[node[i].n_n[1]].xi_3 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].eta_1 + (35.0 / 132.0)*jac[i].eta_1 - (34.0 / 33.0)*jac[node[i].n_n[1]].eta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].eta_2 + (35.0 / 132.0)*jac[i].eta_2 - (34.0 / 33.0)*jac[node[i].n_n[1]].eta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[3]].eta_3 + (35.0 / 132.0)*jac[i].eta_3 - (34.0 / 33.0)*jac[node[i].n_n[1]].eta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_3);
				}
				if (step != 0 && node[i].n_n[3] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_1 + (5.0 / 12.0)*jac[node[i].n_n[1]].xi_1 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_2 + (5.0 / 12.0)*jac[node[i].n_n[1]].xi_2 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_3 + (5.0 / 12.0)*jac[node[i].n_n[1]].xi_3 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].xi_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].xi_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].xi_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_1 + (5.0 / 12.0)*jac[node[i].n_n[1]].eta_1 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_2 + (5.0 / 12.0)*jac[node[i].n_n[1]].eta_2 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_3 + (5.0 / 12.0)*jac[node[i].n_n[1]].eta_3 - (5.0)*jac[node[node[i].n_n[1]].n_n[1]].eta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].eta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[1]].n_n[1]].n_n[1]].n_n[1]].n_n[1]].eta_3);

				}

				if (step >1 && (node[i].ID == 10 || node[node[node[i].n_n[3]].n_n[3]].ID == 10 || node[node[node[i].n_n[1]].n_n[1]].ID == 10 || node[node[i].n_n[3]].ID == 10 || node[node[i].n_n[1]].ID == 10))
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					if (node[i].e[2] != 0 && (CD[node[i].e[2]].connect[1] == gar || CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1] == gar || node[i].ID == 0))
					{
						c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].xi_1 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].xi_1 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].xi_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].xi_1);
						c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].xi_2 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].xi_2 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].xi_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].xi_2);
						c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].xi_3 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].xi_3 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].xi_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].xi_3);
						d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].eta_1 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].eta_1 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].eta_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].eta_1);
						d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].eta_2 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].eta_2 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].eta_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].eta_2);
						d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[1]].connect[0]].e[1]].connect[0]].eta_3 - (7.0 / 9.0)*jac[CD[node[i].e[1]].connect[0]].eta_3 + (7.0 / 9.0)*jac[CD[node[i].e[2]].connect[1]].eta_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[2]].connect[1]].e[2]].connect[1]].eta_3);
					}
					if (node[i].e[3] != 0 && (CD[node[i].e[3]].connect[5] == gar || CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5] == gar || node[i].ID == 0))
					{
						c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].xi_1 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].xi_1 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].xi_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].xi_1);
						c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].xi_2 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].xi_2 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].xi_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].xi_2);
						c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].xi_3 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].xi_3 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].xi_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].xi_3);
						d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].eta_1 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].eta_1 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].eta_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].eta_1);
						d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].eta_2 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].eta_2 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].eta_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].eta_2);
						d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[0]].connect[4]].e[0]].connect[4]].eta_3 - (7.0 / 9.0)*jac[CD[node[i].e[0]].connect[4]].eta_3 + (7.0 / 9.0)*jac[CD[node[i].e[3]].connect[5]].eta_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[3]].connect[5]].e[3]].connect[5]].eta_3);
					}
					if (node[i].e[7] != 0 && (CD[node[i].e[7]].connect[6] == gar || CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6] == gar || node[i].ID == 0))
					{
						c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].xi_1 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].xi_1 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].xi_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].xi_1);
						c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].xi_2 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].xi_2 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].xi_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].xi_2);
						c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].xi_3 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].xi_3 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].xi_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].xi_3);
						d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].eta_1 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].eta_1 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].eta_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].eta_1);
						d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].eta_2 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].eta_2 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].eta_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].eta_2);
						d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[4]].connect[7]].e[4]].connect[7]].eta_3 - (7.0 / 9.0)*jac[CD[node[i].e[4]].connect[7]].eta_3 + (7.0 / 9.0)*jac[CD[node[i].e[7]].connect[6]].eta_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[7]].connect[6]].e[7]].connect[6]].eta_3);
					}
					if (node[i].e[6] != 0 && (CD[node[i].e[6]].connect[2] == gar || CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2] == gar || node[i].ID == 0))
					{
						c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].xi_1 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].xi_1 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].xi_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].xi_1);
						c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].xi_2 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].xi_2 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].xi_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].xi_2);
						c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].xi_3 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].xi_3 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].xi_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].xi_3);
						d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].eta_1 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].eta_1 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].eta_1 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].eta_1);
						d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].eta_2 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].eta_2 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].eta_2 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].eta_2);
						d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[CD[node[CD[node[i].e[5]].connect[3]].e[5]].connect[3]].eta_3 - (7.0 / 9.0)*jac[CD[node[i].e[5]].connect[3]].eta_3 + (7.0 / 9.0)*jac[CD[node[i].e[6]].connect[2]].eta_3 + (1.0 / 36.0)*jac[CD[node[CD[node[i].e[6]].connect[2]].e[6]].connect[2]].eta_3);
					}
				}

				if (node[i].ID == 10)
				{
					gar = i;
				}
				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;
			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);

			c1_dash[M - 1] = c1[M - 1] - ((c1[M] * b[M - 1]) / d[M]);
			c1_dash[M - 2] = c1[M - 2] - ((c1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c2_dash[M - 1] = c2[M - 1] - ((c2[M] * b[M - 1]) / d[M]);
			c2_dash[M - 2] = c2[M - 2] - ((c2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c3_dash[M - 1] = c3[M - 1] - ((c3[M] * b[M - 1]) / d[M]);
			c3_dash[M - 2] = c3[M - 2] - ((c3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d1_dash[M - 1] = d1[M - 1] - ((d1[M] * b[M - 1]) / d[M]);
			d1_dash[M - 2] = d1[M - 2] - ((d1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d2_dash[M - 1] = d2[M - 1] - ((d2[M] * b[M - 1]) / d[M]);
			d2_dash[M - 2] = d2[M - 2] - ((d2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d3_dash[M - 1] = d3[M - 1] - ((d3[M] * b[M - 1]) / d[M]);
			d3_dash[M - 2] = d3[M - 2] - ((d3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				c1_dash[i] = c1[i] - ((c1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c2_dash[i] = c2[i] - ((c2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c3_dash[i] = c3[i] - ((c3_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d1_dash[i] = d1[i] - ((d1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d2_dash[i] = d2[i] - ((d2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d3_dash[i] = d3[i] - ((d3_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jaco[i].zeta_1 = (c1_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_1) / d_dash[step];
					jaco[i].zeta_2 = (c2_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_2) / d_dash[step];
					jaco[i].zeta_3 = (c3_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_3) / d_dash[step];
					jaco[i].zeta_4 = (d1_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_4) / d_dash[step];
					jaco[i].zeta_5 = (d2_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_5) / d_dash[step];
					jaco[i].zeta_6 = (d3_dash[step] - a[step] * jaco[node[i].n_n[ci]].zeta_6) / d_dash[step];

					if (node[i].ID == 10)
					{
						if (node[i].e[2] != 0 && CD[node[i].e[2]].connect[1] == gar)
						{
							jaco[i].zeta_1 = (c1_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_1) / d_dash[step];
							jaco[i].zeta_2 = (c2_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_2) / d_dash[step];
							jaco[i].zeta_3 = (c3_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_3) / d_dash[step];
							jaco[i].zeta_4 = (d1_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_4) / d_dash[step];
							jaco[i].zeta_5 = (d2_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_5) / d_dash[step];
							jaco[i].zeta_6 = (d3_dash[step] - a[step] * jaco[CD[node[i].e[2]].connect[1]].zeta_6) / d_dash[step];
						}
						if (node[i].e[3] != 0 && CD[node[i].e[3]].connect[5] == gar)
						{
							jaco[i].zeta_1 = (c1_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_1) / d_dash[step];
							jaco[i].zeta_2 = (c2_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_2) / d_dash[step];
							jaco[i].zeta_3 = (c3_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_3) / d_dash[step];
							jaco[i].zeta_4 = (d1_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_4) / d_dash[step];
							jaco[i].zeta_5 = (d2_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_5) / d_dash[step];
							jaco[i].zeta_6 = (d3_dash[step] - a[step] * jaco[CD[node[i].e[3]].connect[5]].zeta_6) / d_dash[step];
						}
						if (node[i].e[7] != 0 && CD[node[i].e[7]].connect[6] == gar)
						{
							jaco[i].zeta_1 = (c1_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_1) / d_dash[step];
							jaco[i].zeta_2 = (c2_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_2) / d_dash[step];
							jaco[i].zeta_3 = (c3_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_3) / d_dash[step];
							jaco[i].zeta_4 = (d1_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_4) / d_dash[step];
							jaco[i].zeta_5 = (d2_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_5) / d_dash[step];
							jaco[i].zeta_6 = (d3_dash[step] - a[step] * jaco[CD[node[i].e[7]].connect[6]].zeta_6) / d_dash[step];
						}
						if (node[i].e[6] != 0 && CD[node[i].e[6]].connect[2] == gar)
						{
							jaco[i].zeta_1 = (c1_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_1) / d_dash[step];
							jaco[i].zeta_2 = (c2_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_2) / d_dash[step];
							jaco[i].zeta_3 = (c3_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_3) / d_dash[step];
							jaco[i].zeta_4 = (d1_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_4) / d_dash[step];
							jaco[i].zeta_5 = (d2_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_5) / d_dash[step];
							jaco[i].zeta_6 = (d3_dash[step] - a[step] * jaco[CD[node[i].e[6]].connect[2]].zeta_6) / d_dash[step];
						}
					}
				}
				if (step == 0)
				{
					jaco[i].zeta_1 = c1_dash[step] / d_dash[0];
					jaco[i].zeta_2 = c2_dash[step] / d_dash[0];
					jaco[i].zeta_3 = c3_dash[step] / d_dash[0];
					jaco[i].zeta_4 = d1_dash[step] / d_dash[0];
					jaco[i].zeta_5 = d2_dash[step] / d_dash[0];
					jaco[i].zeta_6 = d3_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jaco[i].zeta_1 = (c1[step] - a[M] * jaco[node[i].n_n[ci]].zeta_1) / d[M];
					jaco[i].zeta_2 = (c2[step] - a[M] * jaco[node[i].n_n[ci]].zeta_2) / d[M];
					jaco[i].zeta_3 = (c3[step] - a[M] * jaco[node[i].n_n[ci]].zeta_3) / d[M];
					jaco[i].zeta_4 = (d1[step] - a[M] * jaco[node[i].n_n[ci]].zeta_4) / d[M];
					jaco[i].zeta_5 = (d2[step] - a[M] * jaco[node[i].n_n[ci]].zeta_5) / d[M];
					jaco[i].zeta_6 = (d3[step] - a[M] * jaco[node[i].n_n[ci]].zeta_6) / d[M];
				}

				i = node[i].n_n[k];

				step++;
			}
		}
		/******************************************************************y-derivative*******************************************************************/
		i = all_boundary_nodes[h];
		if ((node[i].n_n[0] == 0) && node[i].ID != 10)
		{
			if (node[i].n_n[0] == 0)
			{
				ci = 0;
				k = 2;
			}
			step = 0;
			while (i != 0)
			{
				if (step>1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_1 - (7.0 / 9.0)*jac[node[i].n_n[2]].xi_1 + (7.0 / 9.0)*jac[node[i].n_n[0]].xi_1 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_1);
					c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_2 - (7.0 / 9.0)*jac[node[i].n_n[2]].xi_2 + (7.0 / 9.0)*jac[node[i].n_n[0]].xi_2 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_2);
					c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_3 - (7.0 / 9.0)*jac[node[i].n_n[2]].xi_3 + (7.0 / 9.0)*jac[node[i].n_n[0]].xi_3 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_3);
					d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_1 - (7.0 / 9.0)*jac[node[i].n_n[2]].zeta_1 + (7.0 / 9.0)*jac[node[i].n_n[0]].zeta_1 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_2 - (7.0 / 9.0)*jac[node[i].n_n[2]].zeta_2 + (7.0 / 9.0)*jac[node[i].n_n[0]].zeta_2 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_3 - (7.0 / 9.0)*jac[node[i].n_n[2]].zeta_3 + (7.0 / 9.0)*jac[node[i].n_n[0]].zeta_3 + (1.0 / 36.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_3);
				}

				if (step == 0 && node[i].n_n[0] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					c1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_1 + (5.0 / 12.0)*jac[node[i].n_n[2]].xi_1 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_1);
					c2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_2 + (5.0 / 12.0)*jac[node[i].n_n[2]].xi_2 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_2);
					c3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].xi_3 + (5.0 / 12.0)*jac[node[i].n_n[2]].xi_3 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_3);
					d1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_1 + (5.0 / 12.0)*jac[node[i].n_n[2]].zeta_1 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_2 + (5.0 / 12.0)*jac[node[i].n_n[2]].zeta_2 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_3 + (5.0 / 12.0)*jac[node[i].n_n[2]].zeta_3 - (5.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_3);

				}
				if (step == 1)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].xi_1 + (35.0 / 132.0)*jac[i].xi_1 - (34.0 / 33.0)*jac[node[i].n_n[2]].xi_1 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_1);
					c2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].xi_2 + (35.0 / 132.0)*jac[i].xi_2 - (34.0 / 33.0)*jac[node[i].n_n[2]].xi_2 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_2);
					c3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].xi_3 + (35.0 / 132.0)*jac[i].xi_3 - (34.0 / 33.0)*jac[node[i].n_n[2]].xi_3 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].xi_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].xi_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].xi_3);
					d1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].zeta_1 + (35.0 / 132.0)*jac[i].zeta_1 - (34.0 / 33.0)*jac[node[i].n_n[2]].zeta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].zeta_2 + (35.0 / 132.0)*jac[i].zeta_2 - (34.0 / 33.0)*jac[node[i].n_n[2]].zeta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[0]].zeta_3 + (35.0 / 132.0)*jac[i].zeta_3 - (34.0 / 33.0)*jac[node[i].n_n[2]].zeta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[2]].n_n[2]].zeta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].zeta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[2]].n_n[2]].n_n[2]].n_n[2]].zeta_3);

				}

				if (step > 1 && node[node[i].n_n[2]].n_n[2] == 0 && node[i].n_n[2] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].xi_1 + (35.0 / 132.0)*jac[i].xi_1 - (34.0 / 33.0)*jac[node[i].n_n[0]].xi_1 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].xi_2 + (35.0 / 132.0)*jac[i].xi_2 - (34.0 / 33.0)*jac[node[i].n_n[0]].xi_2 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].xi_3 + (35.0 / 132.0)*jac[i].xi_3 - (34.0 / 33.0)*jac[node[i].n_n[0]].xi_3 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].zeta_1 + (35.0 / 132.0)*jac[i].zeta_1 - (34.0 / 33.0)*jac[node[i].n_n[0]].zeta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].zeta_2 + (35.0 / 132.0)*jac[i].zeta_2 - (34.0 / 33.0)*jac[node[i].n_n[0]].zeta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[2]].zeta_3 + (35.0 / 132.0)*jac[i].zeta_3 - (34.0 / 33.0)*jac[node[i].n_n[0]].zeta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_3);

				}

				if (step != 0 && node[i].n_n[2] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_1 + (5.0 / 12.0)*jac[node[i].n_n[0]].xi_1 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_2 + (5.0 / 12.0)*jac[node[i].n_n[0]].xi_2 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].xi_3 + (5.0 / 12.0)*jac[node[i].n_n[0]].xi_3 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].xi_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].xi_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].xi_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_1 + (5.0 / 12.0)*jac[node[i].n_n[0]].zeta_1 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_2 + (5.0 / 12.0)*jac[node[i].n_n[0]].zeta_2 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_3 + (5.0 / 12.0)*jac[node[i].n_n[0]].zeta_3 - (5.0)*jac[node[node[i].n_n[0]].n_n[0]].zeta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].zeta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[0]].n_n[0]].n_n[0]].n_n[0]].n_n[0]].zeta_3);
				}

				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;

			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);
			c1_dash[M - 1] = c1[M - 1] - ((c1[M] * b[M - 1]) / d[M]);
			c1_dash[M - 2] = c1[M - 2] - ((c1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c2_dash[M - 1] = c2[M - 1] - ((c2[M] * b[M - 1]) / d[M]);
			c2_dash[M - 2] = c2[M - 2] - ((c2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c3_dash[M - 1] = c3[M - 1] - ((c3[M] * b[M - 1]) / d[M]);
			c3_dash[M - 2] = c3[M - 2] - ((c3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d1_dash[M - 1] = d1[M - 1] - ((d1[M] * b[M - 1]) / d[M]);
			d1_dash[M - 2] = d1[M - 2] - ((d1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d2_dash[M - 1] = d2[M - 1] - ((d2[M] * b[M - 1]) / d[M]);
			d2_dash[M - 2] = d2[M - 2] - ((d2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d3_dash[M - 1] = d3[M - 1] - ((d3[M] * b[M - 1]) / d[M]);
			d3_dash[M - 2] = d3[M - 2] - ((d3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				c1_dash[i] = c1[i] - ((c1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c2_dash[i] = c2[i] - ((c2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c3_dash[i] = c3[i] - ((c3_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d1_dash[i] = d1[i] - ((d1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d2_dash[i] = d2[i] - ((d2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d3_dash[i] = d3[i] - ((d3_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jaco[i].eta_1 = (c1_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_1) / d_dash[step];
					jaco[i].eta_2 = (c2_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_2) / d_dash[step];
					jaco[i].eta_3 = (c3_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_3) / d_dash[step];
					jaco[i].eta_4 = (d1_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_4) / d_dash[step];
					jaco[i].eta_5 = (d2_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_5) / d_dash[step];
					jaco[i].eta_6 = (d3_dash[step] - a[step] * jaco[node[i].n_n[ci]].eta_6) / d_dash[step];
				}
				if (step == 0)
				{
					jaco[i].eta_1 = c1_dash[step] / d_dash[0];
					jaco[i].eta_2 = c2_dash[step] / d_dash[0];
					jaco[i].eta_3 = c3_dash[step] / d_dash[0];
					jaco[i].eta_4 = d1_dash[step] / d_dash[0];
					jaco[i].eta_5 = d2_dash[step] / d_dash[0];
					jaco[i].eta_6 = d3_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jaco[i].eta_1 = (c1[step] - a[M] * jaco[node[i].n_n[ci]].eta_1) / d[M];
					jaco[i].eta_2 = (c2[step] - a[M] * jaco[node[i].n_n[ci]].eta_2) / d[M];
					jaco[i].eta_3 = (c3[step] - a[M] * jaco[node[i].n_n[ci]].eta_3) / d[M];
					jaco[i].eta_4 = (d1[step] - a[M] * jaco[node[i].n_n[ci]].eta_4) / d[M];
					jaco[i].eta_5 = (d2[step] - a[M] * jaco[node[i].n_n[ci]].eta_5) / d[M];
					jaco[i].eta_6 = (d3[step] - a[M] * jaco[node[i].n_n[ci]].eta_6) / d[M];
				}

				i = node[i].n_n[k];
				step++;
			}
		}
		/******************************************************************z-derivative*******************************************************************/
		i = all_boundary_nodes[h];
		if ((node[i].n_n[4] == 0) && (node[i].loc == 30 || node[i].loc == 31 || node[i].loc <= 27 || node[i].loc == 28))
		{
			if (node[i].n_n[4] == 0)
			{
				ci = 4;
				k = 5;
			}
			step = 0;
			while (i != 0)
			{
				if (step>1)
				{
					b[step] = 1.0 / 3.0;
					d[step] = 1.0;
					a[step] = 1.0 / 3.0;
					c1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_1 - (7.0 / 9.0)*jac[node[i].n_n[5]].eta_1 + (7.0 / 9.0)*jac[node[i].n_n[4]].eta_1 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_1);
					c2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_2 - (7.0 / 9.0)*jac[node[i].n_n[5]].eta_2 + (7.0 / 9.0)*jac[node[i].n_n[4]].eta_2 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_2);
					c3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_3 - (7.0 / 9.0)*jac[node[i].n_n[5]].eta_3 + (7.0 / 9.0)*jac[node[i].n_n[4]].eta_3 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_3);
					d1[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_1 - (7.0 / 9.0)*jac[node[i].n_n[5]].zeta_1 + (7.0 / 9.0)*jac[node[i].n_n[4]].zeta_1 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_2 - (7.0 / 9.0)*jac[node[i].n_n[5]].zeta_2 + (7.0 / 9.0)*jac[node[i].n_n[4]].zeta_2 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((-1.0 / 36.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_3 - (7.0 / 9.0)*jac[node[i].n_n[5]].zeta_3 + (7.0 / 9.0)*jac[node[i].n_n[4]].zeta_3 + (1.0 / 36.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_3);
				}

				if (step == 0 && node[i].n_n[4] == 0)
				{
					b[step] = 5.0;
					d[step] = 1.0;
					a[step] = 0.0;
					c1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_1 + (5.0 / 12.0)*jac[node[i].n_n[5]].eta_1 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_1);
					c2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_2 + (5.0 / 12.0)*jac[node[i].n_n[5]].eta_2 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_2);
					c3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].eta_3 + (5.0 / 12.0)*jac[node[i].n_n[5]].eta_3 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_3);
					d1[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_1 + (5.0 / 12.0)*jac[node[i].n_n[5]].zeta_1 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_2 + (5.0 / 12.0)*jac[node[i].n_n[5]].zeta_2 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((197.0 / 60.0)*jac[i].zeta_3 + (5.0 / 12.0)*jac[node[i].n_n[5]].zeta_3 - (5.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_3);

				}

				if (step == 1)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].eta_1 + (35.0 / 132.0)*jac[i].eta_1 - (34.0 / 33.0)*jac[node[i].n_n[5]].eta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_1);
					c2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].eta_2 + (35.0 / 132.0)*jac[i].eta_2 - (34.0 / 33.0)*jac[node[i].n_n[5]].eta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_2);
					c3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].eta_3 + (35.0 / 132.0)*jac[i].eta_3 - (34.0 / 33.0)*jac[node[i].n_n[5]].eta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].eta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].eta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].eta_3);
					d1[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].zeta_1 + (35.0 / 132.0)*jac[i].zeta_1 - (34.0 / 33.0)*jac[node[i].n_n[5]].zeta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_1);
					d2[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].zeta_2 + (35.0 / 132.0)*jac[i].zeta_2 - (34.0 / 33.0)*jac[node[i].n_n[5]].zeta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_2);
					d3[step] = (1.0 / del_zeta)*((20.0 / 33.0)*jac[node[i].n_n[4]].zeta_3 + (35.0 / 132.0)*jac[i].zeta_3 - (34.0 / 33.0)*jac[node[i].n_n[5]].zeta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[5]].n_n[5]].zeta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].zeta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[5]].n_n[5]].n_n[5]].n_n[5]].zeta_3);
				}

				if (step > 1 && node[node[i].n_n[5]].n_n[5] == 0 && node[i].n_n[5] != 0)
				{
					b[step] = 2.0 / 11.0;
					d[step] = 1.0;
					a[step] = 2.0 / 11.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].eta_1 + (35.0 / 132.0)*jac[i].eta_1 - (34.0 / 33.0)*jac[node[i].n_n[4]].eta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].eta_2 + (35.0 / 132.0)*jac[i].eta_2 - (34.0 / 33.0)*jac[node[i].n_n[4]].eta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].eta_3 + (35.0 / 132.0)*jac[i].eta_3 - (34.0 / 33.0)*jac[node[i].n_n[4]].eta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].zeta_1 + (35.0 / 132.0)*jac[i].zeta_1 - (34.0 / 33.0)*jac[node[i].n_n[4]].zeta_1 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_1 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_1 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].zeta_2 + (35.0 / 132.0)*jac[i].zeta_2 - (34.0 / 33.0)*jac[node[i].n_n[4]].zeta_2 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_2 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_2 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((20.0 / 33.0)*jac[node[i].n_n[5]].zeta_3 + (35.0 / 132.0)*jac[i].zeta_3 - (34.0 / 33.0)*jac[node[i].n_n[4]].zeta_3 + (7.0 / 33.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_3 - (2.0 / 33.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_3 + (1.0 / 132.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_3);
				}

				if (step != 0 && node[i].n_n[5] == 0)
				{
					b[step] = 0.0;
					d[step] = 1.0;
					a[step] = 5.0;
					c1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_1 + (5.0 / 12.0)*jac[node[i].n_n[4]].eta_1 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_1);
					c2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_2 + (5.0 / 12.0)*jac[node[i].n_n[4]].eta_2 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_2);
					c3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].eta_3 + (5.0 / 12.0)*jac[node[i].n_n[4]].eta_3 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].eta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].eta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].eta_3);
					d1[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_1 + (5.0 / 12.0)*jac[node[i].n_n[4]].zeta_1 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_1 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_1 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_1 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_1);
					d2[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_2 + (5.0 / 12.0)*jac[node[i].n_n[4]].zeta_2 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_2 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_2 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_2 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_2);
					d3[step] = (1.0 / del_zeta)*(-1.0)*((197.0 / 60.0)*jac[i].zeta_3 + (5.0 / 12.0)*jac[node[i].n_n[4]].zeta_3 - (5.0)*jac[node[node[i].n_n[4]].n_n[4]].zeta_3 + (5.0 / 3.0)*jac[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].zeta_3 - (5.0 / 12.0)*jac[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_3 + (1.0 / 20.0)*jac[node[node[node[node[node[i].n_n[4]].n_n[4]].n_n[4]].n_n[4]].n_n[4]].zeta_3);
				}

				i = node[i].n_n[k];
				step++;
			}

			M = step - 1;

			d_dash[M - 1] = d[M - 1] - ((b[M - 1] * a[M]) / d[M]);
			d_dash[M - 2] = d[M - 2] - ((b[M - 2] * a[M - 1]) / d_dash[M - 1]);
			c1_dash[M - 1] = c1[M - 1] - ((c1[M] * b[M - 1]) / d[M]);
			c1_dash[M - 2] = c1[M - 2] - ((c1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c2_dash[M - 1] = c2[M - 1] - ((c2[M] * b[M - 1]) / d[M]);
			c2_dash[M - 2] = c2[M - 2] - ((c2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			c3_dash[M - 1] = c3[M - 1] - ((c3[M] * b[M - 1]) / d[M]);
			c3_dash[M - 2] = c3[M - 2] - ((c3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d1_dash[M - 1] = d1[M - 1] - ((d1[M] * b[M - 1]) / d[M]);
			d1_dash[M - 2] = d1[M - 2] - ((d1_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d2_dash[M - 1] = d2[M - 1] - ((d2[M] * b[M - 1]) / d[M]);
			d2_dash[M - 2] = d2[M - 2] - ((d2_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			d3_dash[M - 1] = d3[M - 1] - ((d3[M] * b[M - 1]) / d[M]);
			d3_dash[M - 2] = d3[M - 2] - ((d3_dash[M - 1] * b[M - 2]) / d_dash[M - 1]);

			for (i = M - 3; i >= 0; i--)
			{
				d_dash[i] = d[i] - ((b[i] * a[i + 1]) / d_dash[i + 1]);
				c1_dash[i] = c1[i] - ((c1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c2_dash[i] = c2[i] - ((c2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				c3_dash[i] = c3[i] - ((c3_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d1_dash[i] = d1[i] - ((d1_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d2_dash[i] = d2[i] - ((d2_dash[i + 1] * b[i]) / d_dash[i + 1]);
				d3_dash[i] = d3[i] - ((d3_dash[i + 1] * b[i]) / d_dash[i + 1]);
			}

			i = all_boundary_nodes[h];
			step = 0;
			while (i != 0)
			{
				if (step > 0 && step < M)
				{
					jaco[i].xi_1 = (c1_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_1) / d_dash[step];
					jaco[i].xi_2 = (c2_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_2) / d_dash[step];
					jaco[i].xi_3 = (c3_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_3) / d_dash[step];
					jaco[i].xi_4 = (d1_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_4) / d_dash[step];
					jaco[i].xi_5 = (d2_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_5) / d_dash[step];
					jaco[i].xi_6 = (d3_dash[step] - a[step] * jaco[node[i].n_n[ci]].xi_6) / d_dash[step];
				}
				if (step == 0)
				{
					jaco[i].xi_1 = c1_dash[step] / d_dash[0];
					jaco[i].xi_2 = c2_dash[step] / d_dash[0];
					jaco[i].xi_3 = c3_dash[step] / d_dash[0];
					jaco[i].xi_4 = d1_dash[step] / d_dash[0];
					jaco[i].xi_5 = d2_dash[step] / d_dash[0];
					jaco[i].xi_6 = d3_dash[step] / d_dash[0];
				}

				if (step == M)
				{
					jaco[i].xi_1 = (c1[step] - a[M] * jaco[node[i].n_n[ci]].xi_1) / d[M];
					jaco[i].xi_2 = (c2[step] - a[M] * jaco[node[i].n_n[ci]].xi_2) / d[M];
					jaco[i].xi_3 = (c3[step] - a[M] * jaco[node[i].n_n[ci]].xi_3) / d[M];
					jaco[i].xi_4 = (d1[step] - a[M] * jaco[node[i].n_n[ci]].xi_4) / d[M];
					jaco[i].xi_5 = (d2[step] - a[M] * jaco[node[i].n_n[ci]].xi_5) / d[M];
					jaco[i].xi_6 = (d3[step] - a[M] * jaco[node[i].n_n[ci]].xi_6) / d[M];
				}

				i = node[i].n_n[k];
				step++;
			}
		}
	}

	for (i = 1; i <= NUMNP; i++)
	{
		//det[i] = (((jacobian[i].x_zeta)*(jacobian[i].y_eta*jacobian[i].z_xi-jacobian[i].y_xi*jacobian[i].z_eta))+((jacobian[i].x_eta)*(jacobian[i].y_xi*jacobian[i].z_zeta-jacobian[i].y_zeta*jacobian[i].z_xi))+((jacobian[i].x_xi)*(jacobian[i].y_zeta*jacobian[i].z_eta-jacobian[i].y_eta*jacobian[i].z_zeta)));
		/*		metric[i].zeta_x = (1.0/det[i])*(jaco[i].xi_2-jaco[i].eta_2);
		metric[i].zeta_y = (1.0/det[i])*(jaco[i].xi_3-jaco[i].eta_3);
		metric[i].zeta_z = (1.0/det[i])*(jaco[i].xi_1-jaco[i].eta_1);
		metric[i].eta_x = (1.0/det[i])*(jaco[i].zeta_2-jaco[i].xi_5);
		metric[i].eta_y = (1.0/det[i])*(jaco[i].zeta_3-jaco[i].xi_6);
		metric[i].eta_z = (1.0/det[i])*(jaco[i].zeta_1-jaco[i].xi_4);
		metric[i].xi_x = (1.0/det[i])*(jaco[i].eta_5-jaco[i].zeta_5);
		metric[i].xi_y = (1.0/det[i])*(jaco[i].eta_6-jaco[i].zeta_6);
		metric[i].xi_z = (1.0/det[i])*(jaco[i].eta_4-jaco[i].zeta_4);
		*/
		metric[i].zeta_x = (1.0 / det[i])*(jacobian[i].y_eta*jacobian[i].z_xi - jacobian[i].y_xi*jacobian[i].z_eta);
		metric[i].zeta_y = (1.0 / det[i])*(jacobian[i].z_eta*jacobian[i].x_xi - jacobian[i].z_xi*jacobian[i].x_eta);
		metric[i].zeta_z = (1.0 / det[i])*(jacobian[i].x_eta*jacobian[i].y_xi - jacobian[i].x_xi*jacobian[i].y_eta);
		metric[i].eta_x = (1.0 / det[i])*(jacobian[i].y_xi*jacobian[i].z_zeta - jacobian[i].y_zeta*jacobian[i].z_xi);
		metric[i].eta_y = (1.0 / det[i])*(jacobian[i].z_xi*jacobian[i].x_zeta - jacobian[i].z_zeta*jacobian[i].x_xi);
		metric[i].eta_z = (1.0 / det[i])*(jacobian[i].x_xi*jacobian[i].y_zeta - jacobian[i].x_zeta*jacobian[i].y_xi);
		metric[i].xi_x = (1.0 / det[i])*(jacobian[i].y_zeta*jacobian[i].z_eta - jacobian[i].y_eta*jacobian[i].z_zeta);
		metric[i].xi_y = (1.0 / det[i])*(jacobian[i].z_zeta*jacobian[i].x_eta - jacobian[i].z_eta*jacobian[i].x_zeta);
		metric[i].xi_z = (1.0 / det[i])*(jacobian[i].x_zeta*jacobian[i].y_eta - jacobian[i].x_eta*jacobian[i].y_zeta);


	}

}



