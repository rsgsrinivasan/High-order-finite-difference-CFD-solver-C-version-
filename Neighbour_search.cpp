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

void node_neighbour()
{
	int lenght1, lenght2, lenght, lenghts2;

	for (element = 1; element <= NELEM; element++)
	{
		temp = 0;
		if (fabs(node[CD[element].connect[0]].x - node[CD[element].connect[4]].x) >= fabs(node[CD[element].connect[0]].y - node[CD[element].connect[4]].y) && fabs(node[CD[element].connect[0]].x - node[CD[element].connect[4]].x) >= fabs(node[CD[element].connect[0]].z - node[CD[element].connect[4]].z))
		{
			value = 1;
		}
		if (fabs(node[CD[element].connect[0]].y - node[CD[element].connect[4]].y) >= fabs(node[CD[element].connect[0]].x - node[CD[element].connect[4]].x) && fabs(node[CD[element].connect[0]].y - node[CD[element].connect[4]].y) >= fabs(node[CD[element].connect[0]].z - node[CD[element].connect[4]].z))
		{
			value = 2;
		}
		if (fabs(node[CD[element].connect[0]].z - node[CD[element].connect[4]].z) >= fabs(node[CD[element].connect[0]].x - node[CD[element].connect[4]].x) && fabs(node[CD[element].connect[0]].z - node[CD[element].connect[4]].z) >= fabs(node[CD[element].connect[0]].y - node[CD[element].connect[4]].y))
		{
			value = 3;
		}

		switch (value)
		{
		case 1:

			max_angle = 0;
			lenght1 = 0;
			lenght2 = 0;
			lenght = 0;
			lenghts2 = 0;
			for (i = 0; i<4; i++)
			{
				if (i<3)
				{
					denominator = (node[CD[element].connect[i + 1]].z - node[CD[element].connect[i]].z);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i + 1]].y - node[CD[element].connect[i]].y) / denominator)));
					}
				}
				else
				{
					denominator = (node[CD[element].connect[i]].z - node[CD[element].connect[0]].z);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i]].y - node[CD[element].connect[0]].y) / denominator)));
					}
				}
				if (angle[i] >max_angle)
				{
					max_angle = angle[i];
					lenght1 = i;
				}
			}
			if (lenght1 == 0)
			{
				lenght2 = 2;
			}
			else if (lenght1 == 1)
			{
				lenght2 = 3;
			}
			else if (lenght1 == 2)
			{
				lenght2 = 0;
			}
			else
			{
				lenght2 = 1;
			}
			if (lenght1 == 0)
			{
				if (node[CD[element].connect[lenght1]].y > node[CD[element].connect[lenght1 + 1]].y)
				{
					if (node[CD[element].connect[0]].z > node[CD[element].connect[lenght2 + 1]].z)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[0]].z > node[CD[element].connect[lenght2 + 1]].z)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 2)
			{
				if (node[CD[element].connect[lenght1]].y < node[CD[element].connect[lenght1 + 1]].y)
				{
					if (node[CD[element].connect[lenght1]].z < node[CD[element].connect[1]].z)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].z < node[CD[element].connect[1]].z)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 1)
			{
				lenght = 2;
				lenghts2 = 0;

				if (node[CD[element].connect[lenght1]].y < node[CD[element].connect[lenght]].y)
				{
					if (node[CD[element].connect[lenght1]].z > node[CD[element].connect[lenghts2]].z)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						CD[element].connect[3] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						CD[element].connect[7] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
					else
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].z < node[CD[element].connect[lenghts2]].z)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}
			if (lenght1 == 3)
			{
				if (lenght1 == 3)
				{
					lenght = 0;
					lenghts2 = 2;
				}

				if (node[CD[element].connect[lenght1]].y > node[CD[element].connect[lenght]].y)
				{
					if (node[CD[element].connect[lenght1]].z > node[CD[element].connect[lenghts2]].z)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[0];
						CD[element].connect[0] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[4];
						CD[element].connect[4] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].z > node[CD[element].connect[lenghts2]].z)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}
			P0 = CD[element].connect[0];
			P1 = CD[element].connect[1];
			P2 = CD[element].connect[2];
			P3 = CD[element].connect[3];
			P4 = CD[element].connect[4];
			P5 = CD[element].connect[5];
			P6 = CD[element].connect[6];
			P7 = CD[element].connect[7];

			if (node[CD[element].connect[0]].x - node[CD[element].connect[4]].x <= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P2;
				CD[element].connect[1] = P6;
				CD[element].connect[2] = P5;
				CD[element].connect[3] = P1;
				CD[element].connect[4] = P3;
				CD[element].connect[5] = P7;
				CD[element].connect[6] = P4;
				CD[element].connect[7] = P0;

				temp = 1;
			}
			if (node[CD[element].connect[0]].x - node[CD[element].connect[4]].x >= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P6;
				CD[element].connect[1] = P2;
				CD[element].connect[2] = P1;
				CD[element].connect[3] = P5;
				CD[element].connect[4] = P7;
				CD[element].connect[5] = P3;
				CD[element].connect[6] = P0;
				CD[element].connect[7] = P4;

				temp = 1;
			}
			break;

		case 2:
			max_angle = 0;
			lenght1 = 0;
			lenght2 = 0;
			lenght = 0;
			lenghts2 = 0;
			for (i = 0; i<4; i++)
			{
				if (i<3)
				{
					denominator = (node[CD[element].connect[i + 1]].x - node[CD[element].connect[i]].x);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i + 1]].z - node[CD[element].connect[i]].z) / denominator)));
					}
				}
				else
				{
					denominator = (node[CD[element].connect[i]].x - node[CD[element].connect[0]].x);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i]].z - node[CD[element].connect[0]].z) / denominator)));
					}
				}
				if (angle[i] >max_angle)
				{
					max_angle = angle[i];
					lenght1 = i;
				}
			}
			if (lenght1 == 0)
			{
				lenght2 = 2;
			}
			else if (lenght1 == 1)
			{
				lenght2 = 3;
			}
			else if (lenght1 == 2)
			{
				lenght2 = 0;
			}
			else
			{
				lenght2 = 1;
			}
			if (lenght1 == 0)
			{
				if (node[CD[element].connect[lenght1]].z > node[CD[element].connect[lenght1 + 1]].z)
				{
					if (node[CD[element].connect[0]].x > node[CD[element].connect[lenght2 + 1]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[0]].x > node[CD[element].connect[lenght2 + 1]].x)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 2)
			{
				if (node[CD[element].connect[lenght1]].z < node[CD[element].connect[lenght1 + 1]].z)
				{
					if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[1]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[1]].x)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 1)
			{
				lenght = 2;
				lenghts2 = 0;

				if (node[CD[element].connect[lenght1]].z < node[CD[element].connect[lenght]].z)
				{
					if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						CD[element].connect[3] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						CD[element].connect[7] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
					else
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[lenghts2]].x)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}
			if (lenght1 == 3)
			{
				if (lenght1 == 3)
				{
					lenght = 0;
					lenghts2 = 2;
				}

				if (node[CD[element].connect[lenght1]].z > node[CD[element].connect[lenght]].z)
				{
					if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[0];
						CD[element].connect[0] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[4];
						CD[element].connect[4] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}

			P0 = CD[element].connect[0];
			P1 = CD[element].connect[1];
			P2 = CD[element].connect[2];
			P3 = CD[element].connect[3];
			P4 = CD[element].connect[4];
			P5 = CD[element].connect[5];
			P6 = CD[element].connect[6];
			P7 = CD[element].connect[7];

			if (node[CD[element].connect[0]].y - node[CD[element].connect[4]].y <= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P0;
				CD[element].connect[1] = P1;
				CD[element].connect[2] = P5;
				CD[element].connect[3] = P4;
				CD[element].connect[4] = P3;
				CD[element].connect[5] = P2;
				CD[element].connect[6] = P6;
				CD[element].connect[7] = P7;
				temp = 1;
			}
			if (node[CD[element].connect[0]].y - node[CD[element].connect[4]].y >= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P4;
				CD[element].connect[1] = P5;
				CD[element].connect[2] = P1;
				CD[element].connect[3] = P0;
				CD[element].connect[4] = P7;
				CD[element].connect[5] = P6;
				CD[element].connect[6] = P2;
				CD[element].connect[7] = P3;
				temp = 1;
			}
			break;

		case 3:
			max_angle = 0;
			lenght1 = 0;
			lenght2 = 0;
			lenght = 0;
			lenghts2 = 0;
			for (i = 0; i<4; i++)
			{
				if (i<3)
				{
					denominator = (node[CD[element].connect[i + 1]].x - node[CD[element].connect[i]].x);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i + 1]].y - node[CD[element].connect[i]].y) / denominator)));
					}
				}
				else
				{
					denominator = (node[CD[element].connect[i]].x - node[CD[element].connect[0]].x);
					if (denominator <1e-16 && denominator >= 0.0)
					{
						angle[i] = 90;
					}
					else
					{
						angle[i] = atan(fabs(((node[CD[element].connect[i]].y - node[CD[element].connect[0]].y) / denominator)));
					}
				}
				if (angle[i] >max_angle)
				{
					max_angle = angle[i];
					lenght1 = i;
				}
			}
			if (lenght1 == 0)
			{
				lenght2 = 2;
			}
			else if (lenght1 == 1)
			{
				lenght2 = 3;
			}
			else if (lenght1 == 2)
			{
				lenght2 = 0;
			}
			else
			{
				lenght2 = 1;
			}
			if (lenght1 == 0)
			{
				if (node[CD[element].connect[lenght1]].y > node[CD[element].connect[lenght1 + 1]].y)
				{
					if (node[CD[element].connect[0]].x > node[CD[element].connect[lenght2 + 1]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[0]].x > node[CD[element].connect[lenght2 + 1]].x)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 2)
			{
				if (node[CD[element].connect[lenght1]].y < node[CD[element].connect[lenght1 + 1]].y)
				{
					if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[1]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						shift2 = CD[element].connect[2];
						CD[element].connect[2] = CD[element].connect[1];
						CD[element].connect[1] = shift;
						CD[element].connect[3] = shift2;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						shift2 = CD[element].connect[6];
						CD[element].connect[6] = CD[element].connect[5];
						CD[element].connect[5] = shift;
						CD[element].connect[7] = shift2;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[1];
						CD[element].connect[1] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[5];
						CD[element].connect[5] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[1]].x)
				{
					shift = CD[element].connect[2];
					CD[element].connect[2] = CD[element].connect[0];
					CD[element].connect[0] = shift;

					shift = CD[element].connect[6];
					CD[element].connect[6] = CD[element].connect[4];
					CD[element].connect[4] = shift;
				}
				else
				{
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					shift2 = CD[element].connect[3];
					CD[element].connect[3] = shift;
					CD[element].connect[1] = CD[element].connect[2];
					CD[element].connect[2] = shift2;

					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					shift2 = CD[element].connect[7];
					CD[element].connect[7] = shift;
					CD[element].connect[5] = CD[element].connect[6];
					CD[element].connect[6] = shift2;
				}
			}

			if (lenght1 == 1)
			{
				lenght = 2;
				lenghts2 = 0;

				if (node[CD[element].connect[lenght1]].y < node[CD[element].connect[lenght]].y)
				{
					if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[3];
						CD[element].connect[3] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[7];
						CD[element].connect[7] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
					else
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x < node[CD[element].connect[lenghts2]].x)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}
			if (lenght1 == 3)
			{
				if (lenght1 == 3)
				{
					lenght = 0;
					lenghts2 = 2;
				}

				if (node[CD[element].connect[lenght1]].y > node[CD[element].connect[lenght]].y)
				{
					if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
					{
						shift = CD[element].connect[0];
						CD[element].connect[0] = CD[element].connect[2];
						CD[element].connect[2] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[3];
						CD[element].connect[3] = shift;

						shift = CD[element].connect[4];
						CD[element].connect[4] = CD[element].connect[6];
						CD[element].connect[6] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[7];
						CD[element].connect[7] = shift;
					}
					else
					{
						shift = CD[element].connect[3];
						CD[element].connect[3] = CD[element].connect[0];
						CD[element].connect[0] = shift;
						shift = CD[element].connect[1];
						CD[element].connect[1] = CD[element].connect[2];
						CD[element].connect[2] = shift;

						shift = CD[element].connect[7];
						CD[element].connect[7] = CD[element].connect[4];
						CD[element].connect[4] = shift;
						shift = CD[element].connect[5];
						CD[element].connect[5] = CD[element].connect[6];
						CD[element].connect[6] = shift;
					}
				}
				else if (node[CD[element].connect[lenght1]].x > node[CD[element].connect[lenghts2]].x)
				{
					shift = CD[element].connect[3];
					CD[element].connect[3] = CD[element].connect[2];
					CD[element].connect[2] = shift;
					shift = CD[element].connect[0];
					CD[element].connect[0] = CD[element].connect[1];
					CD[element].connect[1] = shift;

					shift = CD[element].connect[7];
					CD[element].connect[7] = CD[element].connect[6];
					CD[element].connect[6] = shift;
					shift = CD[element].connect[4];
					CD[element].connect[4] = CD[element].connect[5];
					CD[element].connect[5] = shift;
				}
			}

			P0 = CD[element].connect[0];
			P1 = CD[element].connect[1];
			P2 = CD[element].connect[2];
			P3 = CD[element].connect[3];
			P4 = CD[element].connect[4];
			P5 = CD[element].connect[5];
			P6 = CD[element].connect[6];
			P7 = CD[element].connect[7];

			if (node[CD[element].connect[0]].z - node[CD[element].connect[4]].z <= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P7;
				CD[element].connect[1] = P6;
				CD[element].connect[2] = P5;
				CD[element].connect[3] = P4;
				CD[element].connect[4] = P3;
				CD[element].connect[5] = P2;
				CD[element].connect[6] = P1;
				CD[element].connect[7] = P0;
				temp = 1;
			}
			if (node[CD[element].connect[0]].z - node[CD[element].connect[4]].z >= 0.0 && temp == 0)
			{
				CD[element].connect[0] = P3;
				CD[element].connect[1] = P2;
				CD[element].connect[2] = P1;
				CD[element].connect[3] = P0;
				CD[element].connect[4] = P7;
				CD[element].connect[5] = P6;
				CD[element].connect[6] = P5;
				CD[element].connect[7] = P4;
				temp = 1;
			}
			break;
		}
	}

	/**********************************************************************************************/
	/***********************************node connectivity******************************************/
	for (n1 = 1; n1 <= NUMNP; n1++)
	{
		for (element = 0; element <= 7; element++)
		{
			for (j = 0; j < 8; j++)
			{
				if (n1 == CD[temp_node[n1].e[element]].connect[j])
				{
					if (j == 0)
					{
						node[n1].e[2] = temp_node[n1].e[element];
					}
					else if (j == 1)
					{
						node[n1].e[1] = temp_node[n1].e[element];
					}
					else if (j == 2)
					{
						node[n1].e[5] = temp_node[n1].e[element];
					}
					else if (j == 3)
					{
						node[n1].e[6] = temp_node[n1].e[element];
					}
					else if (j == 4)
					{
						node[n1].e[3] = temp_node[n1].e[element];
					}
					else if (j == 5)
					{
						node[n1].e[0] = temp_node[n1].e[element];
					}
					else if (j == 6)
					{
						node[n1].e[4] = temp_node[n1].e[element];
					}
					else if (j == 7)
					{
						node[n1].e[7] = temp_node[n1].e[element];
					}
				}
			}
		}
	}

/*
	for (i = 1; i <= NUMNP; i++)
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
	} */
}
