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

void subdomain_neighbour()
{
	int k;
	/***************GLOBAL AND LOCAL NUMBERING OF NODES AND ITS NEIGHBOURS*************************/
	/**********DOMAIN NEIGHBOUR AND SUBDOMAIN BOUNDARY POINTS PROCS IDENTIFICATION*****************/

	h = sd_node + 1;
	gh = sd_node;
	k = 0;
	for (i = 1; i <= sd_node; i++)
	{
		for (k = 0; k <= 5; k++)
		{
			if (d_node[d_node[node[i].global].n_n[k]].proc != myrank && d_node[node[i].global].n_n[k] != 0 && d_node[d_node[node[i].global].n_n[k]].proc > 0 && d_node[d_node[node[i].global].n_n[k]].proc < size)
			{
				if ((k == 0 || k == 1) && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k + 2] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k + 2] = i;
					}
					if (k == 0)
					{
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					}
					if (k == 1)
					{
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					}
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if ((k == 2 || k == 3) && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k - 2] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k - 2] = i;
					}
					if (k == 2)
					{
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					}
					if (k == 3)
					{
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if (k == 4 && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k + 1] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k + 1] = i;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if (k == 5 && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k - 1] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k - 1] = i;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}
			}


			if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[node[i].global].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc < size)
			{
				if ((k == 0 || k == 1) && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 2] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 2] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if ((k == 2 || k == 3) && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 2] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 2] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if (k == 4 && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 1] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 1] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if (k == 5 && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 1] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 1] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}


				if (k == 0 || k == 1)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 2] = h;
					}
				}

				if (k == 2 || k == 3)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 2] = h;
					}
				}

				if (k == 4)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 1] = h;
					}
				}

				if (k == 5)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 1] = h;
					}
				}
			}

			if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc < size)
			{
				if (k == 0 || k == 1)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 2] = node[node[i].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 2] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 2 || k == 3)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 2] = node[node[i].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 2] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 4)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 1] = node[node[i].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 1] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 5)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[i].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 1] = node[node[i].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 1] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					}
				}


				if (k == 0 || k == 1)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 2] = h;
					}
				}

				if (k == 2 || k == 3)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 2] = h;
					}
				}

				if (k == 4)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 1] = h;
					}
				}

				if (k == 5)
				{
					if (d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 1] = h;
					}
				}
			}
			if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc < size)
			{
				if (k == 0 || k == 1)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 2] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 2] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 2 || k == 3)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 2] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 2] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 4)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 1] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 1] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 5)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local == 0)
					{
						node[node[node[node[i].n_n[k]].n_n[k]].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 1] = node[node[node[i].n_n[k]].n_n[k]].n_n[k];
						node[h].val = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k];
						d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 1] = d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].local;
					}
				}

				if (k == 0 || k == 1)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 2] = h;
					}
				}

				if (k == 2 || k == 3)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 2] = h;
					}
				}

				if (k == 4)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k + 1] = h;
					}
				}

				if (k == 5)
				{
					if (d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].proc == myrank)
					{
						node[h].n_n[k] = d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].n_n[k]].n_n[k]].local].n_n[k - 1] = h;
					}
				}
			}
		}
	}

	nodes = h;
	sd_gh = gh - sd_node;
	all_bou_node = temp;

	for (i = h - sd_gh; i <= nodes; i++)
	{
		for (k = 0; k <= 5; k++)
		{
			if (d_node[d_node[node[i].global].n_n[k]].proc != myrank && d_node[node[i].global].n_n[k] != 0 && d_node[d_node[node[i].global].n_n[k]].proc > 0 && d_node[d_node[node[i].global].n_n[k]].proc < size)
			{
				if ((k == 0 || k == 1) && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k + 2] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k + 2] = i;
					}
					if (k == 0)
					{
						node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					}
					if (k == 1)
					{
						node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					}
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if ((k == 2 || k == 3) && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k - 2] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k - 2] = i;
					}
					if (k == 2)
					{
						node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					}
					if (k == 3)
					{
						node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if (k == 4 && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k + 1] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k + 1] = i;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[5] = d_node[d_node[node[i].global].n_n[5]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}

				if (k == 5 && node[i].n_n[k] == 0)
				{
					if (d_node[d_node[node[i].global].n_n[k]].local == 0)
					{
						node[i].n_n[k] = h;
						node[h].n_n[k - 1] = i;
						node[h].val = d_node[d_node[node[i].global].n_n[k]].val;
						node[h].loc = d_node[d_node[node[i].global].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[node[i].global].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[node[i].global].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[node[i].global].n_n[k]].proc;
						sd_gh_node[h].global = d_node[node[i].global].n_n[k];
						node[h].global = d_node[node[i].global].n_n[k];
						d_node[d_node[node[i].global].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[node[i].global].n_n[k]].local != 0)
					{
						node[i].n_n[k] = d_node[d_node[node[i].global].n_n[k]].local;
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k - 1] = i;
					}
					node[i].n_n[0] = d_node[d_node[node[i].global].n_n[0]].local;
					node[i].n_n[1] = d_node[d_node[node[i].global].n_n[1]].local;
					node[i].n_n[2] = d_node[d_node[node[i].global].n_n[2]].local;
					node[i].n_n[3] = d_node[d_node[node[i].global].n_n[3]].local;
					node[i].n_n[4] = d_node[d_node[node[i].global].n_n[4]].local;
					node[i].val = d_node[node[i].global].val;
					node[i].loc = d_node[node[i].global].loc;
					node[i].corner_ID = d_node[node[i].global].corner_ID;
					node[i].ID = d_node[node[i].global].ID;
				}
			}


			if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc != myrank && d_node[d_node[node[i].global].n_n[k]].n_n[k] != 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc > 0 && d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc < size)
			{
				if ((k == 0 || k == 1) && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 2] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 2] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if ((k == 2 || k == 3) && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 2] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 2] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if (k == 4 && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k + 1] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 1] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}

				if (k == 5 && node[node[i].n_n[k]].n_n[k] == 0)
				{
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local == 0)
					{
						node[node[i].n_n[k]].n_n[k] = h;
						node[h].n_n[k - 1] = node[i].n_n[k];
						node[h].val = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].val;
						node[h].loc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].loc;
						node[h].corner_ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].corner_ID;
						node[h].ID = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].ID;
						sd_gh_node[h].proc = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc;
						sd_gh_node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						node[h].global = d_node[d_node[node[i].global].n_n[k]].n_n[k];
						d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local = h;
						gh++;
						h++;
					}
					if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local != 0)
					{
						node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
						node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 1] = d_node[d_node[node[i].global].n_n[k]].local;
					}
				}
			}

			if (k == 0 || k == 1)
			{
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
				{
					node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 2] = d_node[d_node[node[i].global].n_n[k]].local;
				}
			}

			if (k == 2 || k == 3)
			{
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
				{
					node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 2] = d_node[d_node[node[i].global].n_n[k]].local;
				}
			}

			if (k == 4)
			{
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
				{
					node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k + 1] = d_node[d_node[node[i].global].n_n[k]].local;
				}
			}

			if (k == 5)
			{
				if (d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].proc == myrank)
				{
					node[d_node[d_node[node[i].global].n_n[k]].local].n_n[k] = d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local;
					node[d_node[d_node[d_node[node[i].global].n_n[k]].n_n[k]].local].n_n[k - 1] = d_node[d_node[node[i].global].n_n[k]].local;
				}
			}

		}
	}

	nodes = h;
	sd_gh = gh - sd_node;
	all_bou_node = temp;






}
