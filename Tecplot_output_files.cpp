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

void nodefile()
{
	/***********************************************************************************************************/
	sprintf(filename2, "nodefile.plt");
	outFile.open(filename2);

	outFile << "TITLE = Node file\n";
	outFile << "VARIABLES = X, Y, Z, Proc, Loc,\n";
	outFile << "ZONE NODES= " << NUMNP << " ELEMENTS= " << NELEM << " DATAPACKING=POINT, ZONETYPE=FEBRICK \n";
	for (i = 0; i<NUMNP; i++)
	{
		outFile << node[i + 1].x << "\t" << node[i + 1].y << "\t" << node[i + 1].z << "\t" << node[i + 1].proc << "\t" << node[i + 1].loc << "\n";
	}
	outFile << "\n\n\n";

	for (i = 1; i <= NELEM; i++)
	{
		outFile << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[4] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\t" << CD[i].connect[7] << "\n";
	}
	outFile.close();
	
}



void metric_file()
{
	/***********************************************************************************************************/
	sprintf(filename2, "Grid_Transformation.plt");
	outFile.open(filename2);

	outFile << "TITLE = Node file\n";
	outFile << "VARIABLES = X, Y, Z, Proc, Loc, det \n";
	outFile << "ZONE NODES= " << NUMNP << " ELEMENTS= " << NELEM << " DATAPACKING=POINT, ZONETYPE=FEBRICK \n";
	for (i = 0; i<NUMNP; i++)
	{
		outFile << node[i + 1].x << "\t" << node[i + 1].y << "\t" << node[i + 1].z << "\t" << node[i + 1].proc << "\t" << node[i + 1].loc << "\t" << 1.0 / det[i + 1] << "\n";
	}
	outFile << "\n\n\n";

	for (i = 1; i <= NELEM; i++)
	{
		outFile << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[4] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\t" << CD[i].connect[7] << "\n";
	}
	outFile.close();

}


void proc_request()
{
	ofstream fp;
	sprintf(filename2, "proc_request.plt");
	fp.open(filename2);
	fp << "TITLE = Node file\n";
	fp << "VARIABLES = X, Y, Z, Req, Proc, loc,\n";
	fp << "ZONE NODES= " << NUMNP << " ELEMENTS= " << NELEM << " DATAPACKING=POINT, ZONETYPE=FEBRICK\n";
	for (i = 0; i<NUMNP; i++)
	{
		fp << node[i + 1].x << "\t" << node[i + 1].y << "\t" << node[i + 1].z << "\t" << node[i + 1].req << "\t" << node[i + 1].proc << "\t" << node[i + 1].loc << "\n";
	}
	fp << "\n\n\n";

	for (i = 1; i <= NELEM; i++)
	{
		fp << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[4] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\t" << CD[i].connect[7] << "\n";
	}
	fp.close();
}

void writepltfile()
{
	char writeplt[50] = {};
	ofstream fp;
	int tem = sprintf(writeplt, "nodefile_%d.dat", iter);
	fp.open(writeplt);
	fp << "TITLE = Node file\n";
	fp << "VARIABLES = X, Y, Z, U, V, W, RHO, P, T, e, mu, loc, Proc, DUCROS,\n";
	fp << "ZONE NODES= " << NUMNP << " ELEMENTS= " << NELEM << " DATAPACKING=POINT, ZONETYPE=FEBRICK, SOLUTIONTIME= " << iteration << "\n";
	for (i = 0; i<NUMNP; i++)
	{
		fp << node[i + 1].x << "\t" << node[i + 1].y << "\t" << node[i + 1].z << "\t" << u[0][i + 1] << "\t" << v[0][i + 1] << "\t" << w[0][i + 1] << "\t" << rho[0][i + 1] << "\t" << p[0][i + 1] << "\t" << t[0][i + 1] << "\t" << e[0][i + 1] << "\t" << mu[0][i + 1] << "\t" << node[i + 1].loc << "\t" << node[i + 1].proc << "\t" << DUCROS[i + 1] << "\n";
	}
	fp << "\n\n\n";

	for (i = 1; i <= NELEM; i++)
	{
		fp << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[4] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\t" << CD[i].connect[7] << "\n";
	}

	fp.close();
}

void write_restart_file()
{
	char writerestart[50] = {};
	int tem = sprintf(writerestart, "restart_file_%d.neu", restart_num);
	ofstream fp;
	fp.open(writerestart);

	for (i = 0; i<NUMNP; i++)
	{
		fp << u[0][i + 1] << "\t" << v[0][i + 1] << "\t" << w[0][i + 1] << "\t" << rho[0][i + 1] << "\t" << p[0][i + 1] << "\t" << t[0][i + 1] << "\t" << e[0][i + 1] << "\t" << mu[0][i + 1] << "\t" << iter << "\n";
	}
	fp.close();
}

void writepltfile_jacob()
{
	ofstream fp;
	fp.open("nodefilejacob.plt");
	
	fp << "TITLE = Node file\n";
	fp << "VARIABLES = X, Y, Z, x_zeta, x_eta, y_zeta, y_eta, zeta_x, eta_x, zeta_y, eta_y, det, \n";
	fp << "ZONE NODES= " << NUMNP << " ELEMENTS= " << NELEM << " DATAPACKING=POINT, ZONETYPE=FEBRICK\n";
	for (i = 0; i<NUMNP; i++)
	{
		fp << node[i + 1].x << node[i + 1].y << "\t" << node[i + 1].z << "\t" << jacobian[i + 1].x_zeta << "\t" << jacobian[i + 1].x_eta << "\t" << jacobian[i + 1].y_zeta << "\t" << jacobian[i + 1].y_eta << "\t" << metric[i + 1].zeta_x << "\t" << metric[i + 1].eta_x << "\t" << metric[i + 1].zeta_y << "\t" << metric[i + 1].eta_y << "\t" << det[i + 1] << "\n";
	}
	fp << "\n\n\n";

	for (i = 1; i <= NELEM; i++)
	{
		fp << CD[i].connect[0] << "\t" << CD[i].connect[1] << "\t" << CD[i].connect[5] << "\t" << CD[i].connect[4] << "\t" << CD[i].connect[3] << "\t" << CD[i].connect[2] << "\t" << CD[i].connect[6] << "\t" << CD[i].connect[7] << "\n";
	}

	fp.close();
}


