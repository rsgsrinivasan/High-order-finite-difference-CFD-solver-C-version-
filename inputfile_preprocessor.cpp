#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<strings.h>
#include<string.h>
#include<malloc.h>
#include<limits.h>

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

typedef struct
{
	int n_n[6], inl, wal, bou, out;
	double x, y, z;
} MNODES;

typedef struct
{
	int connect[10];
} ELEMS;

int inputfile_preprocessor()
{
	int T_i, T_j, T_k, NUMNP, NELEM, GAR1, GAR2, GAR3, GAR4, GAR5;
	int ELEM1[10], POS1[10], IB;
	//int inl, wal,bou, out;
	int tot_inl, tot_out, tot_bou, tot_wal;
	char line[200], GARBAGE[20];
	int *inlet, *outlet, *boundary, *wall;
	FILE *fp, *fp2;
	MNODES *node;
	ELEMS *CD; 
	T_i = 0;
	
	char ofname[128];
//	printf("Enter .cfx5 file name\n");
//	scanf("%123s",inputfilename);

	strcpy(ofname,inputfilename);
	strcat(inputfilename,".cfx5");
	strcat(ofname,".neu");
	
	fp = fopen(inputfilename,"rt");
	while(fgets(line, 150, fp) != NULL)
	{
		T_i++;
		if(T_i == 3)
		{
			sscanf(line, "%d %d %d %d %d %d %d", &NUMNP, &GAR1, &GAR2, &NELEM, &GAR3, &GAR4, &GAR5);
			printf("reading CFX5 file\n");
			//printf("Number of nodal points = %d\nNumber of Elements = %d\n", NUMNP, NELEM);
			break;
		}		
	}
	
	node = (MNODES*)malloc(((NUMNP+10))*sizeof(MNODES));
	CD = (ELEMS*)malloc(((NELEM+10))*sizeof(ELEMS));
	
	for (T_i=1; T_i<=NUMNP; T_i++)
	{
		node[T_i].x = 0.0;
		node[T_i].y = 0.0;
		node[T_i].z = 0.0;
		node[T_i].inl = 0;
		node[T_i].out = 0;
		node[T_i].wal = 0;
		node[T_i].bou = 0; 
	}
	
	T_i = 1;
	T_j = 0;
	while(fgets(line, 150, fp) != NULL)
	{		
		sscanf(line, "%lf %lf %lf", &node[T_i].x, &node[T_i].y, &node[T_i].z);
		T_i++;	
		if (T_i>NUMNP)
		{
			break;
		}
	}
	
	T_i = 1;
	T_j = 0;
	while(fgets(line, 150, fp) != NULL)
	{
		sscanf(line, "%d %d %d %d %d %d %d %d", &CD[T_i].connect[1], &CD[T_i].connect[2], &CD[T_i].connect[3], &CD[T_i].connect[4], &CD[T_i].connect[5], &CD[T_i].connect[6], &CD[T_i].connect[7], &CD[T_i].connect[8]);
		T_i++;
		if (T_i > NELEM)
		{
			break;
		}
	}
	
	/***********************************************************************************************************************************************/
	
	tot_inl = 0;
	tot_out = 0;
	tot_bou = 0;
	tot_wal = 0;
	
	GAR5 = 0;
	while(fgets(line, 150, fp) != NULL)
	{		
		GAR1 = 0;
		sscanf(line,"%d %s", &GAR1, &GARBAGE);
		if (strncmp(GARBAGE,"INLET",5)==0)
		{
			GAR5++;
			inl = GAR1;
			inlet = (int *)malloc(((inl*3))*sizeof(int));
			IB = 1;
		}
		
		if (strncmp(GARBAGE,"OUTLET",6)==0)
		{
			GAR5++;
			out = GAR1;
			outlet = (int *)malloc(((out*3))*sizeof(int));
			IB = 2;
		}
		
		if (strncmp(GARBAGE,"WALL",4)==0)
		{
			GAR5++;
			wal = GAR1;
			wall = (int *)malloc(((wal*3))*sizeof(int));
			IB = 3;
		}
		
		if (strncmp(GARBAGE,"BOUNDARY",8)==0)
		{
			GAR5++;
			bou = GAR1;
			boundary = (int *)malloc(((bou*3))*sizeof(int));
			IB = 4;
		}
		
		switch (IB)
		{
			case 1:
				GAR1 = 1;
				T_i = 0;
				GAR2 = 0;
				while(fgets(line, 150, fp) != NULL)
				{
					for (T_k=1; T_k<=5; T_k++)
					{
						ELEM1[T_k] = NELEM+10;
						POS1[T_k] = 10;
					}
					sscanf(line,"%d %d %d %d %d %d %d %d %d %d", &ELEM1[1], &POS1[1], &ELEM1[2], &POS1[2], &ELEM1[3], &POS1[3], &ELEM1[4], &POS1[4], &ELEM1[5], &POS1[5]);
					for (T_k = 1; T_k<=5; T_k++)
					{
						T_i = ELEM1[T_k];
						if (POS1[T_k] == 1 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].inl = 1;
							}
							if (node[CD[T_i].connect[3]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].inl = 1;
							}
							if (node[CD[T_i].connect[7]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].inl = 1;
							}
							if (node[CD[T_i].connect[5]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].inl = 1;
							}						
						}
						
						if (POS1[T_k] == 2 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].inl = 1;
							}
							if (node[CD[T_i].connect[8]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].inl = 1;
							}
							if (node[CD[T_i].connect[6]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].inl = 1;
							}
							if (node[CD[T_i].connect[2]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].inl = 1;
							}						
						}
						
						if (POS1[T_k] == 3 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].inl = 1;
							}
							if (node[CD[T_i].connect[5]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].inl = 1;
							}
							if (node[CD[T_i].connect[6]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].inl = 1;
							}
							if (node[CD[T_i].connect[2]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].inl = 1;
							}						
						}
						
						if (POS1[T_k] == 6 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[7]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].inl = 1;
							}
							if (node[CD[T_i].connect[8]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].inl = 1;
							}
							if (node[CD[T_i].connect[6]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].inl = 1;
							}
							if (node[CD[T_i].connect[5]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].inl = 1;
							}						
						}
						
						if (POS1[T_k] == 5 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].inl = 1;
							}
							if (node[CD[T_i].connect[3]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].inl = 1;
							}
							if (node[CD[T_i].connect[1]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].inl = 1;
							}
							if (node[CD[T_i].connect[2]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].inl = 1;
							}						
						}
						
						if (POS1[T_k] == 4 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].inl = 1;
							}
							if (node[CD[T_i].connect[8]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].inl = 1;
							}
							if (node[CD[T_i].connect[7]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].inl = 1;
							}
							if (node[CD[T_i].connect[3]].inl == 0)
							{
								inlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].inl = 1;
							}						
						}
					}
					
					if (GAR2 == inl)
					{
						tot_inl = GAR1-1;
						break;						
					}
				}	
			break;
			
			case 2:
				GAR1 = 1;
				T_i = 0;
				GAR2 = 0;
				while(fgets(line, 150, fp) != NULL)
				{
					for (T_k=1; T_k<=5; T_k++)
					{
						ELEM1[T_k] = NELEM+10;
						POS1[T_k] = 10;
					}
					sscanf(line,"%d %d %d %d %d %d %d %d %d %d", &ELEM1[1], &POS1[1], &ELEM1[2], &POS1[2], &ELEM1[3], &POS1[3], &ELEM1[4], &POS1[4], &ELEM1[5], &POS1[5]);
					for (T_k = 1; T_k<=5; T_k++)
					{
						T_i = ELEM1[T_k];
						if (POS1[T_k] == 1 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].out = 1;
							}
							if (node[CD[T_i].connect[3]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].out = 1;
							}
							if (node[CD[T_i].connect[7]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].out = 1;
							}
							if (node[CD[T_i].connect[5]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].out = 1;
							}						
						}
						
						if (POS1[T_k] == 2 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].out = 1;
							}
							if (node[CD[T_i].connect[8]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].out = 1;
							}
							if (node[CD[T_i].connect[6]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].out = 1;
							}
							if (node[CD[T_i].connect[2]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].out = 1;
							}						
						}
						
						if (POS1[T_k] == 3 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].out = 1;
							}
							if (node[CD[T_i].connect[5]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].out = 1;
							}
							if (node[CD[T_i].connect[6]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].out = 1;
							}
							if (node[CD[T_i].connect[2]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].out = 1;
							}						
						}
						
						if (POS1[T_k] == 6 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[7]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].out = 1;
							}
							if (node[CD[T_i].connect[8]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].out = 1;
							}
							if (node[CD[T_i].connect[6]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].out = 1;
							}
							if (node[CD[T_i].connect[5]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].out = 1;
							}						
						}
						
						if (POS1[T_k] == 5 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].out = 1;
							}
							if (node[CD[T_i].connect[3]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].out = 1;
							}
							if (node[CD[T_i].connect[1]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].out = 1;
							}
							if (node[CD[T_i].connect[2]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].out = 1;
							}						
						}
						
						if (POS1[T_k] == 4 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].out = 1;
							}
							if (node[CD[T_i].connect[8]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].out = 1;
							}
							if (node[CD[T_i].connect[7]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].out = 1;
							}
							if (node[CD[T_i].connect[3]].out == 0)
							{
								outlet[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].out = 1;
							}						
						}
					}
					
					if (GAR2 == out)
					{
						tot_out = GAR1-1;
						break;						
					}
				}	
			break;
			
			case 3:
				GAR1 = 1;
				T_i = 0;
				GAR2 = 0;
				
				while(fgets(line, 150, fp) != NULL)
				{
					for (T_k=1; T_k<=5; T_k++)
					{
						ELEM1[T_k] = NELEM+10;
						POS1[T_k] = 10;
					}
					sscanf(line,"%d %d %d %d %d %d %d %d %d %d", &ELEM1[1], &POS1[1], &ELEM1[2], &POS1[2], &ELEM1[3], &POS1[3], &ELEM1[4], &POS1[4], &ELEM1[5], &POS1[5]);
					for (T_k = 1; T_k<=5; T_k++)
					{
						T_i = ELEM1[T_k];
						if (POS1[T_k] == 1 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].wal = 1;
							}
							if (node[CD[T_i].connect[3]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].wal = 1;
							}
							if (node[CD[T_i].connect[7]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].wal = 1;
							}
							if (node[CD[T_i].connect[5]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].wal = 1;
							}						
						}
						
						if (POS1[T_k] == 2 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].wal = 1;
							}
							if (node[CD[T_i].connect[8]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].wal = 1;
							}
							if (node[CD[T_i].connect[6]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].wal = 1;
							}
							if (node[CD[T_i].connect[2]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].wal = 1;
							}						
						}
						
						if (POS1[T_k] == 3 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[1]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].wal = 1;
							}
							if (node[CD[T_i].connect[5]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].wal = 1;
							}
							if (node[CD[T_i].connect[6]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].wal = 1;
							}
							if (node[CD[T_i].connect[2]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].wal = 1;
							}						
						}
						
						if (POS1[T_k] == 6 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[7]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].wal = 1;
							}
							if (node[CD[T_i].connect[8]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].wal = 1;
							}
							if (node[CD[T_i].connect[6]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].wal = 1;
							}
							if (node[CD[T_i].connect[5]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].wal = 1;
							}						
						}
						
						if (POS1[T_k] == 5 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].wal = 1;
							}
							if (node[CD[T_i].connect[3]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].wal = 1;
							}
							if (node[CD[T_i].connect[1]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].wal = 1;
							}
							if (node[CD[T_i].connect[2]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].wal = 1;
							}						
						}
						
						if (POS1[T_k] == 4 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2 ++;
							if (node[CD[T_i].connect[4]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].wal = 1;
							}
							if (node[CD[T_i].connect[8]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].wal = 1;
							}
							if (node[CD[T_i].connect[7]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].wal = 1;
							}
							if (node[CD[T_i].connect[3]].wal == 0)
							{
								wall[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].wal = 1;
							}						
						}
					}
					
					if (GAR2 == wal)
					{
						tot_wal = GAR1-1;
						break;						
					}
				}	
			break;
			
			case 4:
				GAR1 = 1;
				T_i = 0;
				GAR2 = 0;
				
				while(fgets(line, 150, fp) != NULL)
				{
					for (T_k=1; T_k<=5; T_k++)
					{
						ELEM1[T_k] = NELEM+10;
						POS1[T_k] = 10;
					}
					sscanf(line,"%d %d %d %d %d %d %d %d %d %d", &ELEM1[1], &POS1[1], &ELEM1[2], &POS1[2], &ELEM1[3], &POS1[3], &ELEM1[4], &POS1[4], &ELEM1[5], &POS1[5]);
					for (T_k = 1; T_k<=5; T_k++)
					{
						T_i = ELEM1[T_k];
						if (POS1[T_k] == 1 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2++;
							if (node[CD[T_i].connect[1]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].bou = 1;
							}
							if (node[CD[T_i].connect[3]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].bou = 1;
							}
							if (node[CD[T_i].connect[7]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].bou = 1;
							}
							if (node[CD[T_i].connect[5]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].bou = 1;
							}						
						}
						
						if (POS1[T_k] == 2 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2++;
							if (node[CD[T_i].connect[4]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].bou = 1;
							}
							if (node[CD[T_i].connect[8]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].bou = 1;
							}
							if (node[CD[T_i].connect[6]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].bou = 1;
							}
							if (node[CD[T_i].connect[2]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].bou = 1;
							}						
						}
						
						if (POS1[T_k] == 3 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2++;
							if (node[CD[T_i].connect[1]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].bou = 1;
							}
							if (node[CD[T_i].connect[5]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].bou = 1;
							}
							if (node[CD[T_i].connect[6]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].bou = 1;
							}
							if (node[CD[T_i].connect[2]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].bou = 1;
							}						
						}
						
						if (POS1[T_k] == 6 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2++;
							if (node[CD[T_i].connect[7]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].bou = 1;
							}
							if (node[CD[T_i].connect[8]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].bou = 1;
							}
							if (node[CD[T_i].connect[6]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[6];
								node[CD[T_i].connect[6]].bou = 1;
							}
							if (node[CD[T_i].connect[5]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[5];
								node[CD[T_i].connect[5]].bou = 1;
							}						
						}
						
						if (POS1[T_k] == 5 && POS1[T_k] < 7 && ELEM1[T_k]<=NELEM)
						{
							GAR2++;
							if (node[CD[T_i].connect[4]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].bou = 1;
							}
							if (node[CD[T_i].connect[3]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].bou = 1;
							}
							if (node[CD[T_i].connect[1]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[1];
								node[CD[T_i].connect[1]].bou = 1;
							}
							if (node[CD[T_i].connect[2]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[2];
								node[CD[T_i].connect[2]].bou = 1;
							}						
						}
						
						if (POS1[T_k] == 4 && POS1[T_k] < 7)
						{
							GAR2++;
							if (node[CD[T_i].connect[4]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[4];
								node[CD[T_i].connect[4]].bou = 1;
							}
							if (node[CD[T_i].connect[8]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[8];
								node[CD[T_i].connect[8]].bou = 1;
							}
							if (node[CD[T_i].connect[7]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[7];
								node[CD[T_i].connect[7]].bou = 1;
							}
							if (node[CD[T_i].connect[3]].bou == 0)
							{
								boundary[GAR1++] = CD[T_i].connect[3];
								node[CD[T_i].connect[3]].bou = 1;
							}						
						}
					}
					
					if (GAR2 == bou)
					{
						tot_bou = GAR1-1;
						break;						
					}
				}	
			break;			
		}
		
		if (GAR5 == 4)
		{
			break;
		}		
	}
	
	fclose(fp);
	
	outFile.open(ofname);
	
	outFile << "\t\t\tCONTROL INFO 2.4.6\n";
	outFile << "**\tGAMBIT NEUTRAL FILE\n";
	outFile << "file\n";
	outFile << "PROGRAM:\tGAMBIT VERSION:\t2.4.6\n\n";
	outFile << "\t\tNUMNP\t\tNELEM\t\tNGRPS\t\tNBSETS\t\tNDFCD\t\tNDFVL\n";
	outFile << "\t\t" << NUMNP << "\t\t" << NELEM << "\t\t100\t\t100\t\t100\t\t100\n";
	outFile << "ENDOFSECTION\n";
	outFile << "\tNODAL COORDINATES 2.4.6\n";
	
	for(T_i=1;T_i<=NUMNP;T_i++)
	{
		outFile << "\t\t\t" << T_i << "\t" << node[T_i].x << "\t" << node[T_i].y << "\t" << node[T_i].z << "\n";
	}		
	outFile << "ENDOFSECTION\n";
	outFile << "\t\tELEMENTS/CELLS 2.4.6\n";
	for(T_i=1;T_i<=NELEM;T_i++)
	{
		outFile << "\t\t" << T_i << "\t4\t8\t\t\t" << CD[T_i].connect[1] << "\t" << CD[T_i].connect[2] << "\t" << CD[T_i].connect[3] << "\t" << CD[T_i].connect[4] << "\t" << CD[T_i].connect[5] << "\t" << CD[T_i].connect[6] << "\t" << CD[T_i].connect[7] << "\n";
		outFile << "\t\t\t\t\t\t\t" << CD[T_i].connect[8] << "\n";
	}
	outFile << "ENDOFSECTION\n";
	outFile << "\t\tELEMENT GROUP 2.4.6\n";
	outFile << "GROUP:\t\t\t\t1 ELEMENTS:\t\t\t22\tMATERIAL:\t\t2\tNFLAGS:\t\t\t1\n";
	outFile << "\t\t\t\t\t\t\tfluid\n";
	outFile << "\t\t0\n";
	outFile << "ENDOFSECTION\n";
	
	outFile << " BOUNDARY CONDITIONS 2.4.6\n";
	outFile << "\t\t\t\t\t\t\t\t\tinlet\t\t\t0\t\t\t" << tot_inl << "\t\t\t0\t\t\t24\n";
	for(T_i=1;T_i<=tot_inl;T_i++)
	{
		outFile << "\t\t\t" << inlet[T_i] << "\n";
	}
	outFile << "ENDOFSECTION\n";
	
	outFile << " BOUNDARY CONDITIONS 2.4.6\n";
	outFile << "\t\t\t\t\t\t\t\t\toutlet\t\t\t0\t\t\t" << tot_out << "\t\t\t0\t\t\t24\n";
	for(T_i=1;T_i<=tot_out;T_i++)
	{
		outFile << "\t\t\t" << outlet[T_i] << "\n";
	}
	outFile << "ENDOFSECTION\n";
	
	outFile << " BOUNDARY CONDITIONS 2.4.6\n";
	outFile << "\t\t\t\t\t\t\t\t\twall\t\t\t0\t\t\t" << tot_wal << "\t\t\t0\t\t\t24\n";
	for(T_i=1;T_i<=tot_wal;T_i++)
	{
		outFile << "\t\t\t" << wall[T_i] << "\n";
	}
	outFile << "ENDOFSECTION\n";
	
	outFile << " BOUNDARY CONDITIONS 2.4.6\n";
	outFile << "\t\t\t\t\t\t\t\t\tboundary\t\t\t0\t\t\t" << tot_bou << "\t\t\t0\t\t\t24\n";
	for(T_i=1;T_i<=tot_bou;T_i++)
	{
		outFile << "\t\t\t" << boundary[T_i] << "\n";
	}
	outFile << "ENDOFSECTION\n";
	
	outFile.close();
	
	printf("File preprocessing***************************************done\n");
	
}
