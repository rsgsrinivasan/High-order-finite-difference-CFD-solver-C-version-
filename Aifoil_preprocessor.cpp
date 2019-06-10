/****************************************************AIFOIL PREPROCESSOR**************************************************************************/
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

void Airfoil_preprocessor()
{
	int k;
	inlet_node = (int *)malloc((NUMNP)*sizeof(int));
	outlet_node = (int *)malloc((NUMNP)*sizeof(int));
	wall_node = (int *)malloc((NUMNP)*sizeof(int));
	boundary_node = (int *)malloc((NUMNP)*sizeof(int));
	
	for (i=1; i<=NUMNP; i++)
	{
		node[i].x = 0.0;
		node[i].y = 0.0;
		node[i].z = 0.0;
	}
	
	i = 1;
	while(fgets(line, 150, fp) != NULL)
	{	
		
		sscanf(line, "%lf", &TEMP[i]);
		i++;	
		if (i>((IMAX+1)*JMAX*KMAX)*3)
		{
			break;
		}
	}
	fclose(fp);
	
	h = 1;
	hh = 1;
	for (dim = 1; dim <= 3; dim++ )
	{
		for (j = 1; j<=JMAX; j++)
		{
			for (k=1; k<=KMAX; k++)
			{
				for (i = 1; i<=IMAX+1; i++)
				{
					if (dim == 1)
					{
						ordered_node[i][j][k].x = TEMP[hh];
						hh++;
					}
					if (dim == 2)
					{
						ordered_node[i][j][k].z = TEMP[hh];
						hh++;
					}
					if (dim == 3)
					{
						ordered_node[i][j][k].y = TEMP[hh];
						hh++;
					}
				}
			}
		}
	}
	
	hh = 1;
	for (k=1; k<=KMAX; k++)
	{
		for (j=1; j<=JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				node[hh].x = ordered_node[i][j][k].x;
				node[hh].y = ordered_node[i][j][k].y;
				node[hh].z = ordered_node[i][j][k].z;
				ordered_node[i][j][k].global = hh;	
				hh++;
			}
		}
	}
	
	hh = 1;
	inl = 0;
	out = 0;
	wal = 0;
	bou = 0;
	for (k=1; k<=KMAX; k++)
	{
		for (j=1; j<=JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				node[ordered_node[i][j][k].global].x = ordered_node[i][j][k].x;
				node[ordered_node[i][j][k].global].y = ordered_node[i][j][k].y;
				node[ordered_node[i][j][k].global].z = ordered_node[i][j][k].z;
				node[ordered_node[i][j][k].global].loc = 0;
				
				if (j != 1 && j < JMAX)
				{
					node[ordered_node[i][j][k].global].n_n[0] = ordered_node[i][j+1][k].global;
					node[ordered_node[i][j][k].global].n_n[2] = ordered_node[i][j-1][k].global;			
				}
				
				if (i != 1 && i < IMAX)
				{
					node[ordered_node[i][j][k].global].n_n[1] = ordered_node[i+1][j][k].global;
					node[ordered_node[i][j][k].global].n_n[3] = ordered_node[i-1][j][k].global;
				}
				if (k != 1 && k < KMAX)
				{
					node[ordered_node[i][j][k].global].n_n[4] = ordered_node[i][j][k+1].global;
					node[ordered_node[i][j][k].global].n_n[5] = ordered_node[i][j][k-1].global;
				}
				
				if (j == 1 )
				{
					node[ordered_node[i][j][k].global].n_n[2] = 0;
					node[ordered_node[i][j][k].global].n_n[0] = ordered_node[i][j+1][k].global;
					wall_node[wal] = ordered_node[i][j][k].global;
					node[wall_node[wal]].ID = 3;
					wal++;
				}					
				if (i == 1 )
				{
					node[ordered_node[i][j][k].global].n_n[3] = ordered_node[IMAX][j][k].global;
					node[ordered_node[i][j][k].global].n_n[1] = ordered_node[i+1][j][k].global;
				}
				if (k == 1)
				{
					node[ordered_node[i][j][k].global].n_n[5] = 0;
					node[ordered_node[i][j][k].global].n_n[4] = ordered_node[i][j][k+1].global;
					boundary_node[bou] = ordered_node[i][j][k].global;
					node[boundary_node[bou]].ID = 4;
					bou++;
				}
				
				if (j == JMAX )
				{
					node[ordered_node[i][j][k].global].n_n[0] = 0;
					node[ordered_node[i][j][k].global].n_n[2] = ordered_node[i][j-1][k].global;
					if(node[ordered_node[i][j][k].global].x <=0)
					{
						inlet_node[inl] = ordered_node[i][j][k].global;
						node[inlet_node[inl]].ID = 1;
						inl++;
					}
					if(node[ordered_node[i][j][k].global].x >=0)
					{
						outlet_node[out] = ordered_node[i][j][k].global;
						node[outlet_node[out]].ID = 2;
						out++;
					}
				}					
				if (i == IMAX )
				{
					node[ordered_node[i][j][k].global].n_n[1] = ordered_node[1][j][k].global;
					node[ordered_node[i][j][k].global].n_n[3] = ordered_node[i-1][j][k].global;
				}
				if (k == KMAX)
				{
					node[ordered_node[i][j][k].global].n_n[4] = 0;
					node[ordered_node[i][j][k].global].n_n[5] = ordered_node[i][j][k-1].global;						
					boundary_node[bou] = ordered_node[i][j][k].global;
					node[boundary_node[bou]].ID = 4;
					bou++;
				}				
			}
		}
	}
	
	for (i=1; i<=NUMNP; i++)
	{
		node[i].corner_ID = 0;
		if(node[i].n_n[0] == 0)
		{
			node[i].loc = 1;
		}
		if(node[i].n_n[1] == 0)
		{
			node[i].loc = 2;
		}
		if(node[i].n_n[2] == 0)
		{
			node[i].loc = 3;
		}
		if(node[i].n_n[3] == 0)
		{
			node[i].loc = 4;
		}
		if(node[i].n_n[4] == 0)
		{
			node[i].loc = 6;
		}
		if(node[i].n_n[5] == 0)
		{
			node[i].loc = 5;
		}
		
		if (node[i].n_n[2] == 0 && node[i].n_n[4] == 0)
		{
			node[i].loc = 15;
		}
		if (node[i].n_n[2] == 0 && node[i].n_n[3] == 0)
		{
			node[i].loc = 16;
		}
		if (node[i].n_n[2] == 0 && node[i].n_n[5] == 0)
		{
			node[i].loc = 17;
		}
		if (node[i].n_n[2] == 0 && node[i].n_n[1] == 0)
		{
			node[i].loc = 18;
		}
		if (node[i].n_n[0] == 0 && node[i].n_n[4] == 0)
		{
			node[i].loc = 19;
		}
		if (node[i].n_n[0] == 0 && node[i].n_n[3] == 0)
		{
			node[i].loc = 20;
		}
		if (node[i].n_n[0] == 0 && node[i].n_n[5] == 0)
		{
			node[i].loc = 21;
		}
		if (node[i].n_n[0] == 0 && node[i].n_n[1] == 0)
		{
			node[i].loc = 22;
		}			
	}
	
	fp = fopen(ofname,"w");
	fprintf(fp,"TITLE = \"Node file\"\n");
	fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\",\n");
	fprintf(fp, "ZONE I=%d, J=%d, , K=%d, DATAPACKING=POINT\n",IMAX,JMAX,KMAX);
	for (k=1; k<=(IMAX)*JMAX*KMAX; k++)
	{
			fprintf(fp, "%lf\t%lf\t%lf\n",node[k].x,node[k].y,node[k].z);
	}
	fclose(fp);
	
	fp = fopen("tecplot_file.dat","w");
	fprintf(fp,"TITLE = \"Node file\"\n");
	fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\",\n");
	fprintf(fp, "ZONE I=%d, J=%d, , K=%d, DATAPACKING=POINT\n",IMAX,JMAX,KMAX);
	for (k=1; k<=KMAX; k++)
	{
		for (j=1; j<=JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				fprintf(fp, "%lf\t%lf\t%lf\n",ordered_node[i][j][k].x,ordered_node[i][j][k].y,ordered_node[i][j][k].z);
			}
		}
	}
	fclose(fp);
	
	/*****************************************************Mesh Neutral write************************************************************/
	fp = fopen(Mesh_file,"w");
	fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",NUMNP,NELEM ,IMAX,JMAX,KMAX,inl,out,bou,wal);
	for (i=1; i<=NUMNP; i++)
	{
		fprintf(fp, "%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", node[i].x, node[i].y, node[i].z, node[i].n_n[0], node[i].n_n[1], node[i].n_n[2], node[i].n_n[3], node[i].n_n[4], node[i].n_n[5], node[i].loc);
	}
	for (k=1; k<KMAX; k++)
	{
		for (j=1; j<JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				if (i < IMAX)
				{
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ordered_node[i+1][j][k].global, ordered_node[i+1][j][k+1].global, ordered_node[i][j][k+1].global, ordered_node[i][j][k].global, ordered_node[i+1][j+1][k].global, ordered_node[i+1][j+1][k+1].global, ordered_node[i][j+1][k+1].global, ordered_node[i][j+1][k].global);
				}
				if (i == IMAX)
				{
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ordered_node[1][j][k].global, ordered_node[1][j][k+1].global, ordered_node[i][j][k+1].global, ordered_node[i][j][k].global, ordered_node[1][j+1][k].global, ordered_node[1][j+1][k+1].global, ordered_node[i][j+1][k+1].global, ordered_node[i][j+1][k].global);
				}
			}
		}
	}
	for (i=0; i<inl; i++)
	{
		fprintf(fp, "%d\n", inlet_node[i]);
	}
	for (i=0; i<out; i++)
	{
		fprintf(fp, "%d\n", outlet_node[i]);
	}
	for (i=0; i<bou; i++)
	{
		fprintf(fp, "%d\n", boundary_node[i]);
	}
	for (i=0; i<wal; i++)
	{
		fprintf(fp, "%d\n", wall_node[i]);
	}
	fclose(fp);
	/********************************************************************************************************************************************/
	
	hh = 1;
	for (k=1; k<KMAX; k++)
	{
		for (j=1; j<JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				if (i < IMAX)
				{
					CD[hh].connect[0] = ordered_node[i+1][j][k].global;
					CD[hh].connect[1] = ordered_node[i+1][j][k+1].global;
					CD[hh].connect[2] = ordered_node[i][j][k+1].global;
					CD[hh].connect[3] = ordered_node[i][j][k].global;
					CD[hh].connect[4] = ordered_node[i+1][j+1][k].global;
					CD[hh].connect[5] = ordered_node[i+1][j+1][k+1].global;
					CD[hh].connect[6] = ordered_node[i][j+1][k+1].global;
					CD[hh].connect[7] = ordered_node[i][j+1][k].global;
				}
				if (i == IMAX)
				{
					CD[hh].connect[0] = ordered_node[1][j][k].global;
					CD[hh].connect[1] = ordered_node[1][j][k+1].global;
					CD[hh].connect[2] = ordered_node[i][j][k+1].global;
					CD[hh].connect[3] = ordered_node[i][j][k].global;
					CD[hh].connect[4] = ordered_node[1][j+1][k].global;
					CD[hh].connect[5] = ordered_node[1][j+1][k+1].global;
					CD[hh].connect[6] = ordered_node[i][j+1][k+1].global;
					CD[hh].connect[7] = ordered_node[i][j+1][k].global;
				}
				hh++;
			}
		}
	}
	
	i = 0;
	sprintf(filename,"nodefile_%d.dat",i);
	fp= fopen(filename,"w");
	fprintf(fp,"TITLE = \"Node file\"\n");
	fprintf(fp, "VARIABLES = \"X\", \"Y\", \"Z\",\n");
	fprintf(fp, "ZONE NODES=%d, ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n",NUMNP,NELEM);
	for(i=1;i<=NUMNP;i++)
	{
		fprintf(fp,"%f\t%f\t%f\n", node[i].x, node[i].y, node[i].z);
	}
	fprintf(fp,"\n\n\n");
	
	for (k=1; k<KMAX; k++)
	{
		for (j=1; j<JMAX; j++)
		{
			for (i=1; i<=IMAX; i++)
			{
				if (i < IMAX)
				{
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ordered_node[i+1][j][k].global, ordered_node[i+1][j][k+1].global, ordered_node[i][j][k+1].global, ordered_node[i][j][k].global, ordered_node[i+1][j+1][k].global, ordered_node[i+1][j+1][k+1].global, ordered_node[i][j+1][k+1].global, ordered_node[i][j+1][k].global);
				}
				if (i == IMAX)
				{
					fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",ordered_node[1][j][k].global, ordered_node[1][j][k+1].global, ordered_node[i][j][k+1].global, ordered_node[i][j][k].global, ordered_node[1][j+1][k].global, ordered_node[1][j+1][k+1].global, ordered_node[i][j+1][k+1].global, ordered_node[i][j+1][k].global);
				}
			}
		}
	}	
	fclose(fp);			
	
	printf("Number of nodal points = %d\nNumber of Elements = %d\n", NUMNP, NELEM);
		
	/************Obtaining number of boundary nodes*******************************************/
	no_of_nodes = 0;
	no_of_nodes = inl + out + wal + bou;
	nodes = no_of_nodes;

	temp_node = (TMP_NODE*)malloc(((no_of_nodes*6)+(NUMNP+10))*sizeof(TMP_NODE));
	nodes = (no_of_nodes*15)+(NUMNP+1);	
	
	if (node == NULL || temp_node == NULL || CD == NULL )
	{
		printf("LINE %d out of memory\n",__LINE__);
	}
	
	corn_node = 0;
	wal_node = wal;	
	inl_node = inl;
	out_node = out;
	bou_node = bou;

	all_bou_node = inl+wal+bou+out+corn_node;

	/***************************************************************AIRFOIL END**********************************************************************************/
}	