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

void node_location()
{
	//int corn_node;

	/*****************	Finding location of node ***************************************************/
	corn_node = 0;
	for (i=1; i<=NUMNP; i++)
	{
		node[i].loc = 0;
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 1;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 2;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 3;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 4;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 5;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 6;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 7;
			node[i].val = node[i].loc;
			node[i].corner_ID = 1;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 8;
			node[i].val = node[i].loc;
			node[i].corner_ID = 2;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 9;
			node[i].val = node[i].loc;
			node[i].corner_ID = 3;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 10;
			node[i].val = node[i].loc;
			node[i].corner_ID = 4;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 11;
			node[i].val = node[i].loc;
			node[i].corner_ID = 5;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 12;
			node[i].val = node[i].loc;
			node[i].corner_ID = 6;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 13;
			node[i].val = node[i].loc;
			node[i].corner_ID = 7;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 14;
			node[i].val = node[i].loc;
			node[i].corner_ID = 8;
			corner_node[corn_node] = i;
			corn_node++;
		}
		if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 15;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 16;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 17;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 18;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 19;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 20;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 21;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 22;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 23;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 24;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 25;
			node[i].val = node[i].loc;
		}
		if(node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 26;
			node[i].val = node[i].loc;
		}
		
	/*	if((node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0 && (CD[node[i].e[1]].connect[5] != CD[node[i].e[2]].connect[4]) && (CD[node[i].e[1]].connect[6] != CD[node[i].e[2]].connect[7])))
		{
			node[i].loc = 27;
		}
		if((node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0 && (CD[node[i].e[2]].connect[0] != CD[node[i].e[6]].connect[3]) && (CD[node[i].e[3]].connect[4] != CD[node[i].e[7]].connect[7])))
		{
			node[i].loc = 28;
		} */
	}		
	
	for(i = 1; i <= NUMNP; i++)
	{
		 if (node[i].ID == 3 && node[i].loc == 0 && node[i].n_n[3] != 0 && ((node[i].n_n[4] != 0 && node[node[i].n_n[4]].ID == 0) || (node[i].n_n[5] != 0 && node[node[i].n_n[5]].ID == 0)) && node[i].loc < 33)
		{
			node[i].loc = 27;
			node[i].ID = 10;
			if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
			{
				node[i].n_n[1] = 0;
			}
		}
		if (node[i].ID == 3 && node[i].n_n[3] != 0 && ((node[i].n_n[0] !=0 && node[i].n_n[2] != 0)) && (node[node[i].n_n[4]].ID == 3 || node[node[i].n_n[5]].ID == 3) && node[i].loc < 33)
		{
			node[i].loc = 28;
			node[i].ID = 10;
			if (node[node[i].n_n[3]].loc == 0)
			{
				node[i].n_n[1] = 0;
			}
		}	 	
		
		if (node[i].ID == 3 && node[node[i].n_n[5]].loc == 0 && node[i].n_n[5] != 0)  
		{
			node[i].loc = 6;
		}
		
		if (node[i].ID == 3 && node[node[i].n_n[4]].loc == 0 && node[i].n_n[4] != 0)
		{
			node[i].loc = 5;
		}
	}
	
	
	
	for (i = 1; i<= NUMNP; i++)
	{
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 33;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 34;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 35;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 36;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 37;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 38;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 39;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 40;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 41;
		}
		if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 42;
		}
		if (node[i].e[0] == 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 43;
		}
		if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 44;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 45;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 46;
		}
		if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 47;
		}
		if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 48;
		}	
		
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 49;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] != 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 50;
		}
		if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 51;
		}
		if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] != 0)
		{
			node[i].loc = 52;
		}

		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 53; 
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 54;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 55;
		}
		if (node[i].e[0] != 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 56;
		}		

		/*special case for cavity*/
		if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] != 0 && node[i].e[7] == 0)
		{
			node[i].loc = 6;
		}
		if (node[i].e[0] == 0 && node[i].e[1] != 0 && node[i].e[2] != 0 && node[i].e[3] == 0 && node[i].e[4] == 0 && node[i].e[5] != 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 6;
		}
		if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] != 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] == 0)
		{
			node[i].loc = 5;
		}
		if (node[i].e[0] != 0 && node[i].e[1] == 0 && node[i].e[2] == 0 && node[i].e[3] != 0 && node[i].e[4] == 0 && node[i].e[5] == 0 && node[i].e[6] == 0 && node[i].e[7] != 0)
		{
			node[i].loc = 5;
		} 				
	}
	/**********************************************************************************************/
	/***********************************node connectivity******************************************/
//			#pragma omp parallel private(j,element, n1) shared(node, CD)
//			{	
//				#pragma omp for
		for(i = 1; i <= NUMNP; i++)
		{
			if (node[i].e[0] != 0)
			{
				node[i].n_n[0] = CD[node[i].e[0]].connect[6];
				node[i].n_n[3] = CD[node[i].e[0]].connect[4];
				node[i].n_n[4] = CD[node[i].e[0]].connect[1];
			}
			if (node[i].e[1] != 0)
			{
				node[i].n_n[0] = CD[node[i].e[1]].connect[2];
				node[i].n_n[3] = CD[node[i].e[1]].connect[0];
				node[i].n_n[5] = CD[node[i].e[1]].connect[5];
			}
			if (node[i].e[2] != 0)
			{
				node[i].n_n[0] = CD[node[i].e[2]].connect[3];
				node[i].n_n[1] = CD[node[i].e[2]].connect[1];
				node[i].n_n[5] = CD[node[i].e[2]].connect[4];
			}
			if (node[i].e[3] != 0)
			{
				node[i].n_n[0] = CD[node[i].e[3]].connect[7];
				node[i].n_n[1] = CD[node[i].e[3]].connect[5];
				node[i].n_n[4] = CD[node[i].e[3]].connect[0];
			}
			if (node[i].e[4] != 0)
			{
				node[i].n_n[2] = CD[node[i].e[4]].connect[5];
				node[i].n_n[3] = CD[node[i].e[4]].connect[7];
				node[i].n_n[4] = CD[node[i].e[4]].connect[2];
			}
			if (node[i].e[5] != 0)
			{
				node[i].n_n[2] = CD[node[i].e[5]].connect[1];
				node[i].n_n[3] = CD[node[i].e[5]].connect[3];
				node[i].n_n[5] = CD[node[i].e[5]].connect[6];
			}
			if (node[i].e[6] != 0)
			{
				node[i].n_n[2] = CD[node[i].e[6]].connect[0];
				node[i].n_n[1] = CD[node[i].e[6]].connect[2];
				node[i].n_n[5] = CD[node[i].e[6]].connect[7];
			}
			if (node[i].e[7] != 0)
			{
				node[i].n_n[2] = CD[node[i].e[7]].connect[4];
				node[i].n_n[1] = CD[node[i].e[7]].connect[6];
				node[i].n_n[4] = CD[node[i].e[7]].connect[3];
			}
		/*	if (node[i].loc >=27 && node[i].loc <=32)
			{
				node[i].n_n[1] = 0;
			} */
		}
		
		/****** special case extended Loc....carefull**********************/
//				#pragma omp for 
		/* for(i = 1; i <= NUMNP; i++)
		{
			 if (node[i].ID == 3 && node[i].loc == 0 && node[i].n_n[3] != 0 && ((node[i].n_n[4] != 0 && node[node[i].n_n[4]].ID == 0) || (node[i].n_n[5] != 0 && node[node[i].n_n[5]].ID == 0)) && node[i].loc < 33)
			{
				node[i].loc = 27;
				node[i].ID = 10;
				if (node[node[i].n_n[3]].loc == 0 || node[node[i].n_n[3]].loc == 0)
				{
					node[i].n_n[1] = 0;
				}
			}
			if (node[i].ID == 3 && node[i].n_n[3] != 0 && ((node[i].n_n[0] !=0 && node[i].n_n[2] != 0)) && (node[node[i].n_n[4]].ID == 3 || node[node[i].n_n[5]].ID == 3) && node[i].loc < 33)
			{
				node[i].loc = 28;
				node[i].ID = 10;
				if (node[node[i].n_n[3]].loc == 0)
				{
					node[i].n_n[1] = 0;
				}
			}	 	
			
			if (node[i].ID == 3 && node[node[i].n_n[5]].loc == 0 && node[i].n_n[5] != 0)  
			{
				node[i].loc = 6;
			}
			
			if (node[i].ID == 3 && node[node[i].n_n[4]].loc == 0 && node[i].n_n[4] != 0)
			{
				node[i].loc = 5;
			}
		} */
//			}
	

	for(i=1; i<=NUMNP; i++)
	{
	
		if (node[i].loc == 0 && node[node[i].n_n[0]].loc == 1)
		{
			node[i].val = 119;
		}
		if (node[i].loc == 0 && node[node[i].n_n[3]].loc == 4)
		{
			node[i].val = 117;
		}
		if (node[i].loc == 0 && node[node[i].n_n[1]].loc == 2)
		{
			node[i].val = 115;
		}
		if (node[i].loc == 0 && node[node[i].n_n[5]].loc == 5)
		{
			node[i].val = 116;
		}
		if (node[i].loc == 0 && node[node[i].n_n[4]].loc == 6)
		{
			node[i].val = 118;
		}
		if (node[i].loc == 0 && node[node[i].n_n[2]].loc == 3)
		{
			node[i].val = 120;
		}
	}
	
	for(i=1; i<=NUMNP; i++)
	{
		if (node[i].loc == 0)
		{			
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[5]].loc == 5 )
			{ 
				node[i].val =107;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2)
			{
				node[i].val =108;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[4]].loc == 6)
			{
				node[i].val =109;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4)
			{
				node[i].val =110;
			}
			
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[5]].loc == 5)
			{
				node[i].val =111;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2)
			{
				node[i].val =112;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[4]].loc == 6)
			{
				node[i].val =113;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4)
			{
				node[i].val =114;
			}
				
			
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[4]].loc == 6 )
			{
				node[i].val=99;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[5]].loc == 5 )
			{
				node[i].val=100;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[5]].loc == 5 )
			{
				node[i].val=101;
			}
			if (node[node[i].n_n[0]].loc == 1 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[4]].loc == 6 )
			{
				node[i].val=102;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[4]].loc == 6 )
			{
				node[i].val=103;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[3]].loc == 4 && node[node[i].n_n[5]].loc == 5 )
			{
				node[i].val=104;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[5]].loc == 5 )
			{
				node[i].val=105;
			}
			if (node[node[i].n_n[2]].loc == 3 && node[node[i].n_n[1]].loc == 2 && node[node[i].n_n[4]].loc == 6 )
			{
				node[i].val=106;
			}
		}
		
		if (node[i].loc != 0)
		{	
			if (node[node[i].n_n[2]].loc == 15 && node[node[i].n_n[3]].loc == 23 )
			{ 
				node[i].val =43;
			}
			if (node[node[i].n_n[2]].loc == 15 && node[node[i].n_n[1]].loc == 26)
			{
				node[i].val =50;
			}
			if (node[node[i].n_n[1]].loc == 26 && node[node[i].n_n[0]].loc == 19)
			{
				node[i].val =66;
			}
			if (node[node[i].n_n[0]].loc == 19 && node[node[i].n_n[3]].loc == 23)
			{
				node[i].val =59;
			}
				
			if (node[node[i].n_n[2]].loc == 16 && node[node[i].n_n[4]].loc == 23)
			{
				node[i].val =44;
			}
			if (node[node[i].n_n[4]].loc == 23 && node[node[i].n_n[0]].loc == 20)
			{
				node[i].val =60;
			}
			if (node[node[i].n_n[5]].loc == 24 && node[node[i].n_n[0]].loc == 20)
			{
				node[i].val =61;
			}
			if (node[node[i].n_n[5]].loc == 24 && node[node[i].n_n[2]].loc == 16)
			{
				node[i].val =45;
			}
			
			if (node[node[i].n_n[2]].loc == 17 && node[node[i].n_n[3]].loc == 24)
			{
				node[i].val =46;
			}
			if (node[node[i].n_n[3]].loc == 24 && node[node[i].n_n[0]].loc == 21)
			{
				node[i].val =62;
			}
			if (node[node[i].n_n[0]].loc == 21 && node[node[i].n_n[1]].loc == 25)
			{
				node[i].val =63;
			}
			if (node[node[i].n_n[2]].loc == 17 && node[node[i].n_n[1]].loc == 25)
			{
				node[i].val =47;
			}
				
			if (node[node[i].n_n[5]].loc == 25 && node[node[i].n_n[0]].loc == 22)
			{
				node[i].val =64;
			}
			if (node[node[i].n_n[4]].loc == 26 && node[node[i].n_n[0]].loc == 22)
			{
				node[i].val =65;
			}
			if (node[node[i].n_n[4]].loc == 26 && node[node[i].n_n[2]].loc == 18)
			{
				node[i].val =49;
			}
			if (node[node[i].n_n[5]].loc == 25 && node[node[i].n_n[2]].loc == 18)
			{
				node[i].val =48;
			}
				
			if (node[node[i].n_n[4]].loc == 19 && node[node[i].n_n[3]].loc == 20)
			{
				node[i].val =67;
			}
			if (node[node[i].n_n[3]].loc == 20 && node[node[i].n_n[5]].loc == 21)
			{
				node[i].val =68;
			}
			if (node[node[i].n_n[5]].loc == 21 && node[node[i].n_n[1]].loc == 22)
			{
				node[i].val =69;
			}
			if (node[node[i].n_n[1]].loc == 22 && node[node[i].n_n[4]].loc == 19)
			{
				node[i].val =70;
			}
			
			if (node[node[i].n_n[4]].loc == 15 && node[node[i].n_n[3]].loc == 16)
			{
				node[i].val =71;
			}
			if (node[node[i].n_n[3]].loc == 16 && node[node[i].n_n[5]].loc == 17)
			{
				node[i].val =72;
			}
			if (node[node[i].n_n[5]].loc == 17 && node[node[i].n_n[1]].loc == 18)
			{
				node[i].val =73;
			}
			if (node[node[i].n_n[4]].loc == 15 && node[node[i].n_n[1]].loc == 18)
			{
				node[i].val =74;
			}
			
			
			if (node[node[i].n_n[2]].loc == 16 && node[i].loc == 4 && node[node[i].n_n[4]].loc == 4 && node[node[i].n_n[5]].loc == 4)
			{
				node[i].val = 86;
			}
			if (node[node[i].n_n[2]].loc == 16 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
			{
				node[i].val = 77;
			}

			if (node[node[i].n_n[2]].loc == 17 && node[i].loc == 26 && node[node[i].n_n[1]].loc == 5 && node[node[i].n_n[3]].loc == 5)
			{
				node[i].val = 85;
			}
			if (node[node[i].n_n[5]].loc == 17 && node[i].loc == 3 && node[node[i].n_n[1]].loc == 3 && node[node[i].n_n[3]].loc == 3)
			{
				node[i].val = 93;
			}

			if (node[node[i].n_n[2]].loc == 18 && node[i].loc == 2 && node[node[i].n_n[4]].loc == 2 && node[node[i].n_n[5]].loc == 2)
			{
				node[i].val = 84;
			}					
			if (node[node[i].n_n[1]].loc == 18 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
			{
				node[i].val = 81;
			}					

			if (node[node[i].n_n[4]].loc == 19 && node[i].loc == 1 && node[node[i].n_n[1]].loc == 1 && node[node[i].n_n[3]].loc == 1)
			{
				node[i].val = 95;
			}
			if (node[node[i].n_n[0]].loc == 19 && node[i].loc == 6 && node[node[i].n_n[1]].loc == 6 && node[node[i].n_n[3]].loc == 6)
			{
				node[i].val = 87;
			}

			if (node[node[i].n_n[3]].loc == 20 && node[i].loc == 1 && node[node[i].n_n[4]].loc == 1 && node[node[i].n_n[5]].loc == 1)
			{
				node[i].val = 75;
			}
			if (node[node[i].n_n[0]].loc == 20 && node[i].loc == 4 && node[node[i].n_n[4]].loc == 4 && node[node[i].n_n[5]].loc == 4)
			{
				node[i].val = 90;
			}

			if (node[node[i].n_n[5]].loc == 21 && node[i].loc == 1 && node[node[i].n_n[1]].loc == 1 && node[node[i].n_n[3]].loc == 1)
			{
				node[i].val = 91;
			}
			if (node[node[i].n_n[0]].loc == 21 && node[i].loc == 5 && node[node[i].n_n[1]].loc == 5 && node[node[i].n_n[3]].loc == 5)
			{
				node[i].val = 89;
			}

			if (node[node[i].n_n[1]].loc == 22 && node[i].loc == 1 && node[node[i].n_n[4]].loc == 1 && node[node[i].n_n[5]].loc == 1)
			{
				node[i].val = 79;
			}
			if (node[node[i].n_n[0]].loc == 22 && node[i].loc == 3 && node[node[i].n_n[4]].loc == 3 && node[node[i].n_n[5]].loc == 3)
			{
				node[i].val = 88;
			}

			if (node[node[i].n_n[4]].loc == 23 && node[i].loc == 4 && node[node[i].n_n[0]].loc == 4 && node[node[i].n_n[2]].loc == 4)
			{
				node[i].val = 96;
			}
			if (node[node[i].n_n[3]].loc == 23 && node[i].loc == 6 && node[node[i].n_n[0]].loc == 6 && node[node[i].n_n[2]].loc == 6)
			{
				node[i].val = 78;
			}

			if (node[node[i].n_n[5]].loc == 24 && node[i].loc == 4 && node[node[i].n_n[0]].loc == 4 && node[node[i].n_n[2]].loc == 4)
			{
				node[i].val = 94;
			}
			if (node[node[i].n_n[3]].loc == 24 && node[i].loc == 5 && node[node[i].n_n[0]].loc == 5 && node[node[i].n_n[2]].loc == 5)
			{
				node[i].val = 76;
			}

			if (node[node[i].n_n[1]].loc == 25 && node[i].loc == 5 && node[node[i].n_n[0]].loc == 5 && node[node[i].n_n[2]].loc == 5)
			{
				node[i].val = 80;
			}
			if (node[node[i].n_n[5]].loc == 25 && node[i].loc == 2 && node[node[i].n_n[0]].loc == 2 && node[node[i].n_n[2]].loc == 2)
			{
				node[i].val = 92;
			}

			if (node[node[i].n_n[4]].loc == 26 && node[i].loc == 2 && node[node[i].n_n[0]].loc == 2 && node[node[i].n_n[2]].loc == 2)
			{
				node[i].val = 98;
			}
			if (node[node[i].n_n[1]].loc == 26 && node[i].loc == 6 && node[node[i].n_n[0]].loc == 6 && node[node[i].n_n[2]].loc == 6)
			{
				node[i].val = 82;
			}	
				
			if (node[node[i].n_n[3]].loc == 12 && node[i].loc == 17)
			{
				node[i].val = 38;
			}
			if (node[node[i].n_n[5]].loc == 12 && node[i].loc == 16)
			{
				node[i].val = 37;
			}
			
			if (node[node[i].n_n[4]].loc == 11 && node[i].loc == 16)
			{
				node[i].val = 36;
			}
			if (node[node[i].n_n[3]].loc == 11 && node[i].loc == 15)
			{
				node[i].val = 35;
			}				
				
			if (node[node[i].n_n[1]].loc == 14 && node[i].loc == 15)
			{
				node[i].val = 42;
			}
			if (node[node[i].n_n[4]].loc == 14 && node[i].loc == 18)
			{
				node[i].val = 41;
			}
				
			if (node[node[i].n_n[5]].loc == 13 && node[i].loc == 18)
			{
				node[i].val = 40;
			}
			if (node[node[i].n_n[3]].loc == 13 && node[i].loc == 17)
			{
				node[i].val = 39;
			}
				
			if (node[node[i].n_n[3]].loc == 8 && node[i].loc == 21)
			{
				node[i].val = 54;
			}
			if (node[node[i].n_n[5]].loc == 8 && node[i].loc == 20)
			{
				node[i].val = 53;
			}
				
			if (node[node[i].n_n[4]].loc == 7 && node[i].loc == 20)
			{
				node[i].val = 52;
			}
			if (node[node[i].n_n[3]].loc == 7 && node[i].loc == 19)
			{
				node[i].val = 51;
			}
				
			if (node[node[i].n_n[1]].loc == 10 && node[i].loc == 19)
			{
				node[i].val = 58;
			}
			if (node[node[i].n_n[4]].loc == 10 && node[i].loc == 22)
			{
				node[i].val = 57;
			}
			
			if (node[node[i].n_n[1]].loc == 9 && node[i].loc == 21)
			{
				node[i].val = 55;
			}
			if (node[node[i].n_n[5]].loc == 9 && node[i].loc == 22)
			{
				node[i].val = 56;
			}
				
			if (node[i].n_n[3] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[0]].loc == 7)
			{
				node[i].val = 27;
			}
			if (node[i].n_n[3] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[0]].loc == 8)
			{
				node[i].val = 28;
			}
			if (node[i].n_n[1] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[0]].loc == 9)
			{
				node[i].val = 29;
			}
			if (node[i].n_n[1] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[0]].loc == 10)
			{
				node[i].val = 30;
			}
				
			if (node[i].n_n[3] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[2]].loc == 11)
			{
				node[i].val = 31;
			}
			if (node[i].n_n[3] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[2]].loc == 12)
			{
				node[i].val = 32;
			}
			if (node[i].n_n[1] == 0 && node[i].n_n[5] == 0 && node[node[i].n_n[2]].loc == 13)
			{
				node[i].val = 33;
			}
			if (node[i].n_n[1] == 0 && node[i].n_n[4] == 0 && node[node[i].n_n[2]].loc == 14)
			{
				node[i].val = 34;
			}
		}												
	}
	/***** val near singular points **************/
	for (i=1; i<= NUMNP; i++)
	{			
		if (node[i].loc == 27 || node[i].loc == 28 || node[i].loc == 29 || node[i].loc == 30 || node[i].loc == 31 || node[i].loc == 32) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[4]].val = 116;
		}
		if (node[i].loc == 27 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[1]].val = 117;
			node[node[i].n_n[3]].val = 115;
		}
		if (node[i].loc == 28 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[0]].val = 120;
			node[node[i].n_n[2]].val = 119;
		}
		if (node[i].loc == 29 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[0]].val = 120;
			node[node[i].n_n[4]].val = 116;
			node[node[i].n_n[1]].val = 115;
		}
		if (node[i].loc == 30 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[0]].val = 120;
			node[node[i].n_n[1]].val = 117;
			node[node[i].n_n[4]].val = 116;
		}
		if (node[i].loc == 31 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[2]].val = 119;
			node[node[i].n_n[1]].val = 117;
			node[node[i].n_n[4]].val = 116;
		}
		if (node[i].loc == 32 ) 
		{
			node[i].corner_ID = 0;
			node[node[i].n_n[2]].val = 119;
			node[node[i].n_n[3]].val = 115;
			node[node[i].n_n[4]].val = 116;
		}
		if (node[node[i].n_n[4]].loc == 27 && node[i].n_n[1] == 0)
		{
			node[i].val = 96;
		}
		if (node[node[i].n_n[4]].loc == 27 && node[i].n_n[3] == 0)
		{
			node[i].val = 98;
		}
		if (node[node[i].n_n[4]].loc == 28 && node[i].n_n[0] == 0)
		{
			node[i].val = 95;
		}
		if (node[node[i].n_n[4]].loc == 28 && node[i].n_n[2] == 0)
		{
			node[i].val = 97;
		}
		
		if (node[node[i].n_n[4]].loc == 29 && node[i].loc == 20)
		{
			node[i].val = 52;
		}
		if (node[node[i].n_n[4]].loc == 30 && node[i].loc == 22)
		{
			node[i].val = 57;
		}
		if (node[node[i].n_n[4]].loc == 31 && node[i].loc == 18)
		{
			node[i].val = 41;
		}
		if (node[node[i].n_n[4]].loc == 32 && node[i].loc == 16)
		{
			node[i].val = 36;
		}
		
		if (node[node[i].n_n[4]].loc == 29 && node[i].loc == 18)
		{
			node[i].val = 41;
		}
		if (node[node[i].n_n[4]].loc == 30 && node[i].loc == 16)
		{
			node[i].val = 36;
		}
		if (node[node[i].n_n[4]].loc == 31 && node[i].loc == 20)
		{
			node[i].val = 52;
		}
		if (node[node[i].n_n[4]].loc == 32 && node[i].loc == 22)
		{
			node[i].val = 57;
		}
		
		if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[0]].loc == 20)
		{
			node[i].val = 60;
		}
		if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[0]].loc == 22)
		{
			node[i].val = 65;
		}
		if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[2]].loc == 16)
		{
			node[i].val = 44;
		}
		if (node[node[i].n_n[4]].loc == 27 && node[node[i].n_n[2]].loc == 18)
		{
			node[i].val = 49;
		}
		
		if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[3]].loc == 20)
		{
			node[i].val = 67;
		}
		if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[1]].loc == 22)
		{
			node[i].val = 70;
		}
		if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[3]].loc == 16)
		{
			node[i].val = 71;
		}
		if (node[node[i].n_n[4]].loc == 28 && node[node[i].n_n[1]].loc == 18)
		{
			node[i].val = 74;
		}
		if (node[node[i].n_n[3]].loc == 27 && node[node[i].n_n[0]].loc == 28)
		{
			node[i].val = 59;
		}
		if (node[node[i].n_n[1]].loc == 27 && node[node[i].n_n[0]].loc == 28)
		{
			node[i].val = 66;
		}
		if (node[node[i].n_n[3]].loc == 27 && node[node[i].n_n[2]].loc == 28)
		{
			node[i].val = 43;
		}
		if (node[node[i].n_n[1]].loc == 27 && node[node[i].n_n[2]].loc == 28)
		{
			node[i].val = 50;
		}
		
		if (node[node[i].n_n[3]].loc == 29 && node[i].loc == 19)
		{
			node[i].val = 51;
		}
		if (node[node[i].n_n[1]].loc == 30 && node[i].loc == 19)
		{
			node[i].val = 58;
		}
		if (node[node[i].n_n[3]].loc == 32 && node[i].loc == 15)
		{
			node[i].val = 35;
		}
		if (node[node[i].n_n[1]].loc == 31 && node[i].loc == 15)
		{
			node[i].val = 42;
		}
		if (node[node[i].n_n[0]].loc == 29)
		{
			node[i].val = 27;
		}
		if (node[node[i].n_n[0]].loc == 30)
		{
			node[i].val = 30;
		}
		if (node[node[i].n_n[2]].loc == 31)
		{
			node[i].val = 34;
		}
		if (node[node[i].n_n[2]].loc == 32)
		{
			node[i].val = 31;
		}
	}
	
	all_bou_node = inl+wal+bou+out+corn_node;
	
	for (i=1; i<= NUMNP; i++)
	{
		if (node[i].ID == 10)
		{
			node[i].n_n[1] = 0;
		}
		
	}

}
