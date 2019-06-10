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
using std::cerr;

void calculate_mesh()
{
	i = 0;
	int k = 0;
	int x;
	ifstream infile;
	vector<int> variables;

	infile.open(filename);
	if (!infile.is_open())
	{
		cerr << "* _____ File not found!* \n";
		exit(1);
	}
	else
	{
		while (infile.good())
		{
			getline(infile, line2, '\n');
			if (i == 6)
			{
				istringstream streamA(line2);
				while (streamA >> x)
				{
					variables.push_back(x);
				}
				break;
			}
			i = i + 1;
		}
	}
	infile.close();

	NUMNP = variables[0];
	NELEM = variables[1];
	NGRPS = variables[2];
	NBSETS = variables[3];
}

/*
void AIRFOIL_calculate_mesh()
{
	i = 0;
	NUMNP=0;
	NELEM=0;
	fp = fopen(filename,"rt");
	while(fgets(line, 150, fp) != NULL)
	{
		i++;
		if(i == 2)
		{
			sscanf(line, "%d %d %d", &IMAX, &KMAX, &JMAX);
			IMAX = IMAX-1;
			NUMNP = (IMAX)*JMAX*KMAX;
			NELEM = (IMAX)*(JMAX-1)*(KMAX-1);
			printf("Number of nodal points = %d\nNumber of Elements = %d\n", NUMNP, NELEM);
			break;
		}		
	}
	NBSETS = 4;
}
*/

void restart_file()
{
	ifstream fp;
	int count1, gar,k;
	double x;
	itera1 = -1;
	itera2 = -1;
	string fileline;
	vector<int> variables;
	variables.resize(11);
	count1 = 1;
	fp.open("restart_file_1.neu");
	if (fp.is_open())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		for (k=1; k<=9; k++)
		{
			streamA >> x;
			variables[k] = x;
		}
		fp.close();
		itera = variables[9];
		itera1 = itera;
	}
	
	count1 = 1;
	fp.open("restart_file_2.neu");
	if (fp.is_open())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		for (k = 1; k<=9; k++)
		{
			streamA >> x;
			variables[k] = x;
		}
		fp.close();
		itera = variables[9];
		itera2 = itera;
	}
	

	if (itera1 > itera2 && itera1 >= 0 && itera2 >= 0)
	{
		restart_num = 1;
		gar1 = itera1;
	}
	else if (itera2 > itera1 && itera1 >= 0 && itera2 >= 0)
	{
		restart_num = 2;
		gar1 = itera2;
	}
	gar = restart_num;

	sprintf(filename2, "restart_file_%d.neu", restart_num);
	fp.open(filename2);
	if (fp.is_open())
	{
		i = 1;
		count1 = 1;
		while (fp.good())
		{
			getline(fp, fileline, '\n');
			istringstream streamA(fileline);
			for (k=1; k<=9; k++)
			{
				streamA >> x;
				variables[k] = x;
			}
			itera = variables[9];
			i++;
			if (gar1 != itera)
			{
				i--;
				if (myrank == 0)
				{
					printf("breaking\n");
				}
				if(i > NUMNP  )
				{
					break;
				} 
			}
		}
		itera = -1;				
		if (i <= NUMNP && gar == 1)
		{
			restart_num++;
		}
		if (i <= NUMNP && gar == 2)
		{
			restart_num--;
		}
	}
}


void boundary_nodes_count()
{
	ifstream fp;
	int k;
	no_of_nodes = 0;
	string fileline;
	fp.open(filename);

	no_of_nodes = 0;
	inl = 0;
	out = 0;
	wal = 0;
	bou = 0;
	out_mem = 1;
	out_nl = 0;
	out_old = 0;
	inl_mem = 1;
	inl_nl = 0;
	inl_old = 0;
	wal_mem = 1;
	wal_nl = 0;
	wal_old = 0;
	bou_mem = 1;
	bou_nl = 0;
	bou_old = 0;
	all_bou = 0;
	bound = 0;
	cor = 0;

	while (fp.good())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		streamA >> garbage1;
		if (strcmp(garbage1, "BOUNDARY") == 0)
		{
			break;
		}
	}
//	getline(fp, line, '\n');

	for (k = 0; k<NBSETS; k++)
	{
		while (fp.good())
		{
			getline(fp, fileline, '\n');
			istringstream streamA(fileline);
			streamA >> garbage1;
			streamA >> gar1;
			streamA >> NE;

			if (strncmp(garbage1, "inlet", 5) == 0)
			{
				no_of_nodes = no_of_nodes + NE;
			}

			else if (strncmp(garbage1, "outlet", 6) == 0)
			{
				no_of_nodes = no_of_nodes + NE;
			}

			else if (strncmp(garbage1, "wall", 4) == 0)
			{
				no_of_nodes = no_of_nodes + NE;
			}

			else if (strncmp(garbage1, "boundary", 8) == 0)
			{
				no_of_nodes = no_of_nodes + NE;

			}


			if (strncmp(garbage1, "inlet", 5) == 0)
			{
				inl_old = inl_old + NE;
				IBCODE = 1;
				if (inl_mem == inl_nl)
				{
					inlet_node.resize((NE + 10) + (inl_old));
					inl_mem++;
				}
				else
				{
					inlet_node.resize(NE + 10);
				}
				inl_nl++;
			}

			else if (strncmp(garbage1, "outlet", 6) == 0)
			{
				out_old = out_old + NE;
				IBCODE = 2;
				if (out_mem == out_nl)
				{
					outlet_node.resize((NE + 10) + (out_old));
					out_mem++;
				}
				else
				{
					outlet_node.resize(NE + 10);
				}
				out_nl++;
			}

			else if (strncmp(garbage1, "wall", 4) == 0)
			{
				//wal_mem++;
				wal_old = wal_old + NE;
				IBCODE = 3;
				if (wal_mem == wal_nl)
				{
					wall_node.resize((NE + 10) + (wal_old));
					wal_mem++;
				}
				else
				{
					wall_node.resize(NE + 10);
				}
				wal_nl++;
			}

			else if (strncmp(garbage1, "boundary", 8) == 0)
			{

				bou_old = bou_old + NE;
				IBCODE = 4;
				if (bou_mem == bou_nl)
				{
					boundary_node.resize((NE + 10) + (bou_old));
					bou_mem++;
				}
				else
				{
					boundary_node.resize(NE + 10);
				}
				bou_nl++;
			}
		}		
	}
	fp.close();
}


void read_mesh_file()
{
	ifstream fp;
	i = 1;
	j = 0;
	int k;
	fp.open(filename);
	string fileline;

	while (fp.good())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		streamA >> garbage1;
		if (strcmp(garbage1, "NODAL") == 0)
		{
			break;
		}
	}

	/*********************Reading coordinates****************************************************/
	i =1;
	while (fp.good())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		streamA >> garbage1;
		if (strcmp(garbage1, "ENDOFSECTION") == 0)
		{
			break;
		}
		else
		{
			for (h = 1; h <= 3; h++)
			{
				if (h == 1)
					streamA >> node[i].x;
				else if (h == 2)
					streamA >> node[i].y;
				else if (h == 3)
					streamA >> node[i].z;
			}
			i++;
		}
	}
	getline(fp, fileline, '\n');
	
	/*****************************storing node numbers of each element****************************/
	i = 1;
	j = 1;
	while (fp.good())
	{
		for (j=1; j <=2; j++ )
		{
			getline(fp, fileline, '\n');
			istringstream streamA(fileline);
			streamA >> garbage1;
			stringstream geek(garbage1); 
			if (strcmp(garbage1, "ENDOFSECTION") == 0)
			{
				break;
			}
			else
			{
				for (h = 1; h < 10; h++)
				{
					if (j == 1 && h == 1)
					{
						streamA >> garbage1;
						streamA >> garbage1;
					}
					if (icemcfd == 1)
					{ 
						if (h == 1 && j == 2)
							geek >> CD[i].connect[2];
						else if (h == 3 && j == 1)
							streamA >> CD[i].connect[4];
						else if (h == 4 && j == 1)
							streamA >> CD[i].connect[5];
						else if (h == 5 && j == 1)
							streamA >> CD[i].connect[7];
						else if (h == 6 && j == 1)
							streamA >> CD[i].connect[6];
						else if (h == 7 && j == 1)
							streamA >> CD[i].connect[0];
						else if (h == 8 && j == 1)
							streamA >> CD[i].connect[1];
						else if (h == 9 && j == 1)
							streamA >> CD[i].connect[3];
					}
					else if (gambit == 1)
					{
						if (h == 1 && j == 2)
							geek >> CD[i].connect[6];
						else if (h == 3 && j == 1)
							streamA >> CD[i].connect[0];
						else if (h == 4 && j == 1)
							streamA >> CD[i].connect[1];
						else if (h == 5 && j == 1)
							streamA >> CD[i].connect[3];
						else if (h == 6 && j == 1)
							streamA >> CD[i].connect[2];
						else if (h == 7 && j == 1)
							streamA >> CD[i].connect[4];
						else if (h == 8 && j == 1)
							streamA >> CD[i].connect[5];
						else if (h == 9 && j == 1)
							streamA >> CD[i].connect[7];
					}
				}							
			}
		}
		temp_node[CD[i].connect[0]].e[node[CD[i].connect[0]].val++] = i;
		temp_node[CD[i].connect[1]].e[node[CD[i].connect[1]].val++] = i;
		temp_node[CD[i].connect[2]].e[node[CD[i].connect[2]].val++] = i;
		temp_node[CD[i].connect[3]].e[node[CD[i].connect[3]].val++] = i;
		temp_node[CD[i].connect[4]].e[node[CD[i].connect[4]].val++] = i;
		temp_node[CD[i].connect[5]].e[node[CD[i].connect[5]].val++] = i;
		temp_node[CD[i].connect[6]].e[node[CD[i].connect[6]].val++] = i;
		temp_node[CD[i].connect[7]].e[node[CD[i].connect[7]].val++] = i;

		i++;
		if (i == NELEM+1)
		{
			break;
		}
	}
	/**********************************************************************************************/
	/******************************************Reading boundary nodes from file*************************************/
	inl = 0;
	out = 0;
	wal = 0;
	bou = 0;
	out_mem = 1;
	out_nl = 0;
	out_old = 0;
	inl_mem = 1;
	inl_nl = 0;
	inl_old = 0;
	wal_mem = 1;
	wal_nl = 0;
	wal_old = 0;
	bou_mem = 1;
	bou_nl = 0;
	bou_old = 0;
	all_bou = 0;
	bound = 0;
	cor = 0;

	while (fp.good())
	{
		getline(fp, fileline, '\n');
		istringstream streamA(fileline);
		streamA >> garbage1;
		if (strcmp(garbage1, "BOUNDARY") == 0)
		{
			break;
		}
	}
	//getline(fp, fileline, '\n');

	for (k = 0; k<NBSETS; k++)
	{
		
		getline(fp, fileline, '\n');
		if (k>0)
		{
			getline(fp, fileline, '\n');
		}
		istringstream streamA(fileline);
		streamA >> garbage1;
		streamA >> gar1;
		streamA >> NE;

		if (strcmp(garbage1, "inlet") == 0)
		{
			inl_old = inl_old + NE;
			IBCODE = 1;
		}

		else if (strcmp(garbage1, "outlet") == 0)
		{
			out_old = out_old + NE;
			IBCODE = 2;
		}

		else if (strcmp(garbage1, "wall") == 0)
		{
			wal_old = wal_old + NE;
			IBCODE = 3;
		}

		else if (strcmp(garbage1, "boundary") == 0)
		{
			bou_old = bou_old + NE;
			IBCODE = 4;
		}


		switch (IBCODE)
		{
		case 1:
			while (fp.good())
			{
				getline(fp, fileline, '\n');
				istringstream streamA(fileline);
				streamA >> garbage1;

				if (strcmp(garbage1, "ENDOFSECTION") == 0)
				{
					break;
				}
				else
				{
					stringstream geek(garbage1);
					geek >> inlet_node[inl];
				}
				if (node[inlet_node[inl]].ID == 0)
				{
					node[inlet_node[inl]].ID = 1;
				}
				else
				{
					if (node[inlet_node[inl]].ID != 3)
					{
						node[inlet_node[inl]].ID = 1;
					}
				}
				inl++;
			}
			break;

		case 2:
			while (fp.good())
			{
				getline(fp, fileline, '\n');
				istringstream streamA(fileline);
				streamA >> garbage1;

				if (strcmp(garbage1, "ENDOFSECTION") == 0)
				{
					break;
				}
				else
				{
					stringstream geek(garbage1);
					geek >> outlet_node[out];
				}
				if (node[outlet_node[out]].ID == 0)
				{
					node[outlet_node[out]].ID = 2;
				}
				else
				{
					cor++;
				}
				out++;
			}
			break;

		case 3:
			while (fp.good())
			{
				getline(fp, fileline, '\n');
				istringstream streamA(fileline);
				streamA >> garbage1;

				if (strcmp(garbage1, "ENDOFSECTION") == 0)
				{
					break;
				}
				else
				{
					stringstream geek(garbage1);
					geek >> wall_node[wal];
				}
				if (node[wall_node[wal]].ID == 0)
				{
					node[wall_node[wal]].ID = 3;
				}
				else
				{
					node[wall_node[wal]].ID=3;
					cor++;
				}
				wal++;
			}
			break;

		case 4:
			while (fp.good())
			{
				getline(fp, fileline, '\n');
				istringstream streamA(fileline);
				streamA >> garbage1;

				if (strcmp(garbage1, "ENDOFSECTION") == 0)
				{
					break;
				}
				else
				{
					stringstream geek(garbage1);
					geek >> boundary_node[bou];
				}
				if (node[boundary_node[bou]].ID == 0)
				{
					node[boundary_node[bou]].ID = 4;
				}
				else
				{
					cor++;
				}
				bou++;
			}
			IBCODE = 10;
			break;
		}			
	}
	fp.close();
}





