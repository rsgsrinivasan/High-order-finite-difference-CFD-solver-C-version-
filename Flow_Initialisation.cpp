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

void intialise(int RK)
{
	int l, k, d_temp, d_temp6;
	double Cd, Cdr1, Cdr2, Cdr3, Cdc1, Cdc2, Cdc3, l_determinant, l_angle, l_angle2, Vt, Vn, Vt2, Vn2, Vn_sign;
	Cd = 0.07;
	Cdr1 = 0.0123; 
	Cdr2 = 0.0279;
	Cdr3 = 0.0264;
	Cdc1 = 0.0121;
	Cdc2 = 0.0138;
	Cdc3 = 0.0;
	
	if (initial == 0)
	{
		for (i=1; i<=sd_node; i++)
		{				
			u[RK][i] = 1.0;
			v[RK][i] = 0.0;
			w[RK][i] = 0.0;
			p[RK][i] = 1.0/(1.4*Mach*Mach);
			t[RK][i] = 1.0;
			rho[RK][i] = (1.4*Mach*Mach)*(p[RK][i]/t[RK][i]);
			e[RK][i] = p[RK][i]/(0.4*rho[RK][i]);
			a[RK][i] = sqrt(1.4*p[RK][i]/rho[RK][i]);																
		}		
	}
	
/****************************************************OUTLET********************************************************/		
	for (i=0; i<out_node; i++)
	{	
		if (node[outlet_node[i]].loc > 0 && node[outlet_node[i]].loc <= 6)
		{	
			if (node[outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[RK][outlet_node[i]]= u[RK][node[outlet_node[i]].n_n[k]];	
			v[RK][outlet_node[i]]= v[RK][node[outlet_node[i]].n_n[k]];	
			w[RK][outlet_node[i]]= w[RK][node[outlet_node[i]].n_n[k]];	
			p[RK][outlet_node[i]]= p[RK][node[outlet_node[i]].n_n[k]];
		/* 	if (back_pressure != 0.0)
			{
				if (node[outlet_node[i]].y < 412.0 || node[outlet_node[i]].y > 767.0)
				{
					p[RK][outlet_node[i]]= p[RK][node[outlet_node[i]].n_n[k]];	
				}
				if (node[outlet_node[i]].y > 412.0 && node[outlet_node[i]].y < 767.0)
				{
					p[RK][outlet_node[i]]= back_pressure*(1.0/(1.4*Mach*Mach));	
				}
			}	 */			
			t[RK][outlet_node[i]]= t[RK][node[outlet_node[i]].n_n[k]];	
			rho[RK][outlet_node[i]]= (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
			a[RK][outlet_node[i]] = sqrt(1.4*p[RK][outlet_node[i]]/rho[RK][outlet_node[i]]); 

			u[RK][node[outlet_node[i]].n_n[l]] = u[RK][outlet_node[i]];
			v[RK][node[outlet_node[i]].n_n[l]] = v[RK][outlet_node[i]];
			w[RK][node[outlet_node[i]].n_n[l]] = w[RK][outlet_node[i]];
			p[RK][node[outlet_node[i]].n_n[l]] = p[RK][outlet_node[i]];
			rho[RK][node[outlet_node[i]].n_n[l]] = rho[RK][outlet_node[i]];
			t[RK][node[outlet_node[i]].n_n[l]] = t[RK][outlet_node[i]];	
			a[RK][node[outlet_node[i]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[outlet_node[i]].n_n[l]] = e[RK][outlet_node[i]];
			
			u[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = u[RK][outlet_node[i]];
			v[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = v[RK][outlet_node[i]];
			w[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = w[RK][outlet_node[i]];
			p[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = p[RK][outlet_node[i]];
			rho[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = rho[RK][outlet_node[i]];
			t[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = t[RK][outlet_node[i]];	
			a[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = e[RK][outlet_node[i]];
	
			u[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][outlet_node[i]];
			v[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][outlet_node[i]];
			w[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][outlet_node[i]];
			p[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][outlet_node[i]];
			rho[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][outlet_node[i]];
			t[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][outlet_node[i]];	
			a[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][outlet_node[i]];
			
			 
		/* 	u[RK][node[outlet_node[i]].n_n[l]] = 2.0*u[RK][outlet_node[i]]-u[RK][node[outlet_node[i]].n_n[k]];
			v[RK][node[outlet_node[i]].n_n[l]] = 2.0*v[RK][outlet_node[i]]-v[RK][node[outlet_node[i]].n_n[k]];
			w[RK][node[outlet_node[i]].n_n[l]] = 2.0*w[RK][outlet_node[i]]-w[RK][node[outlet_node[i]].n_n[k]];
			p[RK][node[outlet_node[i]].n_n[l]] = 2.0*p[RK][outlet_node[i]]-p[RK][node[outlet_node[i]].n_n[k]];
			rho[RK][node[outlet_node[i]].n_n[l]] = 2.0*rho[RK][outlet_node[i]]-rho[RK][node[outlet_node[i]].n_n[k]];
			t[RK][node[outlet_node[i]].n_n[l]] = 2.0*t[RK][outlet_node[i]]-t[RK][node[outlet_node[i]].n_n[k]];	
			a[RK][node[outlet_node[i]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[outlet_node[i]].n_n[l]] = e[RK][outlet_node[i]];
			
			u[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*u[RK][node[outlet_node[i]].n_n[l]]-u[RK][outlet_node[i]];
			v[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*v[RK][node[outlet_node[i]].n_n[l]]-v[RK][outlet_node[i]];
			w[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*w[RK][node[outlet_node[i]].n_n[l]]-w[RK][outlet_node[i]];
			p[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*p[RK][node[outlet_node[i]].n_n[l]]-p[RK][outlet_node[i]];
			rho[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*rho[RK][node[outlet_node[i]].n_n[l]]-rho[RK][outlet_node[i]];
			t[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = 2.0*t[RK][node[outlet_node[i]].n_n[l]]-t[RK][outlet_node[i]];
			a[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]] = e[RK][outlet_node[i]];
	
			u[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*u[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-u[RK][node[outlet_node[i]].n_n[l]];
			v[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*v[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-v[RK][node[outlet_node[i]].n_n[l]];
			w[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*w[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-w[RK][node[outlet_node[i]].n_n[l]];
			p[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*p[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-p[RK][node[outlet_node[i]].n_n[l]];
			rho[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*rho[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-rho[RK][node[outlet_node[i]].n_n[l]];
			t[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = 2.0*t[RK][node[node[outlet_node[i]].n_n[l]].n_n[l]]-t[RK][node[outlet_node[i]].n_n[l]];
			a[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][outlet_node[i]];
			e[RK][node[node[node[outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][outlet_node[i]];
			 */
		}
		else if (node[outlet_node[i]].loc == 7 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);

		}
		else if (node[outlet_node[i]].loc == 8 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 9 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 10 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		
		else if (node[outlet_node[i]].loc == 11 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 12 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 13 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);			
		}
		else if (node[outlet_node[i]].loc == 14 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[node[outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);

		}
		else if (node[outlet_node[i]].loc == 15 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[5]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 16 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[1]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 17 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[4]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 18 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[3]].n_n[0]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		
		else if (node[outlet_node[i]].loc == 19 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[5]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 20 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[1]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 21 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[4]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 22 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[3]].n_n[2]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 23 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[1]].n_n[5]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 24 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[1]].n_n[4]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 25 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[3]].n_n[4]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
		else if (node[outlet_node[i]].loc == 26 )
		{
			u[RK][outlet_node[i]] = u[RK][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			v[RK][outlet_node[i]] = v[RK][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			w[RK][outlet_node[i]] = w[RK][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			p[RK][outlet_node[i]] = p[RK][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			t[RK][outlet_node[i]] = t[RK][node[node[outlet_node[i]].n_n[3]].n_n[5]];
			rho[RK][outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][outlet_node[i]]/t[RK][outlet_node[i]]);
			e[RK][outlet_node[i]] = p[RK][outlet_node[i]]/(0.4*rho[RK][outlet_node[i]]);
		}
	}

/****************************************************INLET********************************************************/
	for (i=0; i< inl_node; i++)
	{	
		if (node[inlet_node[i]].loc > 0 && node[inlet_node[i]].loc <= 6)
		{	
			if (node[inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			u[RK][inlet_node[i]]= 1.0;		
			v[RK][inlet_node[i]]= 0.0;	
			w[RK][inlet_node[i]]= 0.0;				
			p[RK][inlet_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[RK][inlet_node[i]]= 1.0;		
			rho[RK][inlet_node[i]]= (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
			a[RK][inlet_node[i]] = sqrt(1.4*p[RK][inlet_node[i]]/rho[RK][inlet_node[i]]); 

			u[RK][node[inlet_node[i]].n_n[l]] = u[RK][inlet_node[i]];
			v[RK][node[inlet_node[i]].n_n[l]] = v[RK][inlet_node[i]];
			w[RK][node[inlet_node[i]].n_n[l]] = w[RK][inlet_node[i]];
			p[RK][node[inlet_node[i]].n_n[l]] = p[RK][inlet_node[i]];
			rho[RK][node[inlet_node[i]].n_n[l]] = rho[RK][inlet_node[i]];
			t[RK][node[inlet_node[i]].n_n[l]] = t[RK][inlet_node[i]];	
			a[RK][node[inlet_node[i]].n_n[l]] = a[RK][inlet_node[i]];
			e[RK][node[inlet_node[i]].n_n[l]] = e[RK][inlet_node[i]];
			//mu[RK][node[inlet_node[i]].n_n[l]] = mu[RK][inlet_node[i]];
			
			u[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = u[RK][inlet_node[i]];
			v[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = v[RK][inlet_node[i]];
			w[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = w[RK][inlet_node[i]];
			p[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = p[RK][inlet_node[i]];
			rho[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = rho[RK][inlet_node[i]];
			t[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = t[RK][inlet_node[i]];	
			a[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = a[RK][inlet_node[i]];
			e[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = e[RK][inlet_node[i]];
			//mu[RK][node[node[inlet_node[i]].n_n[l]].n_n[l]] = mu[RK][inlet_node[i]];
	
			u[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][inlet_node[i]];
			v[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][inlet_node[i]];
			w[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][inlet_node[i]];
			p[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][inlet_node[i]];
			rho[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][inlet_node[i]];
			t[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][inlet_node[i]];	
			a[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][inlet_node[i]];
			e[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][inlet_node[i]];
			//mu[RK][node[node[node[inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][inlet_node[i]];
		}
		else if (node[inlet_node[i]].loc == 7 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 8 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 9 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 10 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		
		else if (node[inlet_node[i]].loc == 11 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 12 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 13 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);			
		}
		else if (node[inlet_node[i]].loc == 14 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[node[inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 15 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[5]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 16 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[1]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 17 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[4]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 18 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[3]].n_n[0]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		
		else if (node[inlet_node[i]].loc == 19 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[5]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 20 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[1]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 21 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[4]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 22 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[3]].n_n[2]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 23 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[1]].n_n[5]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 24 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[1]].n_n[4]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 25 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[3]].n_n[4]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
		else if (node[inlet_node[i]].loc == 26 )
		{
			u[RK][inlet_node[i]] = u[RK][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			v[RK][inlet_node[i]] = v[RK][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			w[RK][inlet_node[i]] = w[RK][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			p[RK][inlet_node[i]] = p[RK][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			t[RK][inlet_node[i]] = t[RK][node[node[inlet_node[i]].n_n[3]].n_n[5]];
			rho[RK][inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][inlet_node[i]]/t[RK][inlet_node[i]]);
			e[RK][inlet_node[i]] = p[RK][inlet_node[i]]/(0.4*rho[RK][inlet_node[i]]);
		}
	}
	
	/****************************************************BOUNDARY********************************************************/
	for (i=0; i< bou_node; i++)
	{	
		if (node[boundary_node[i]].loc > 0 && node[boundary_node[i]].loc <= 6)
		{	
			if (node[boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

/*			u[RK][boundary_node[i]]= 1.0;		
			v[RK][boundary_node[i]]= 0.0;	
			w[RK][boundary_node[i]]= 0.0;				
			p[RK][boundary_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[RK][boundary_node[i]]= 1.0;		
			rho[RK][boundary_node[i]]= (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
			a[RK][boundary_node[i]] = sqrt(1.4*p[RK][boundary_node[i]]/rho[RK][boundary_node[i]]); 
*/
			u[RK][boundary_node[i]] = u[RK][node[boundary_node[i]].n_n[k]];
			v[RK][boundary_node[i]] = v[RK][node[boundary_node[i]].n_n[k]];
			w[RK][boundary_node[i]] = w[RK][node[boundary_node[i]].n_n[k]];
			p[RK][boundary_node[i]] = p[RK][node[boundary_node[i]].n_n[k]];
			rho[RK][boundary_node[i]] = rho[RK][node[boundary_node[i]].n_n[k]];
			t[RK][boundary_node[i]] = t[RK][node[boundary_node[i]].n_n[k]];	
			a[RK][boundary_node[i]] = a[RK][node[boundary_node[i]].n_n[k]];
			e[RK][boundary_node[i]] = e[RK][node[boundary_node[i]].n_n[k]];
			//mu[RK][boundary_node[i]] = mu[RK][boundary_node[i]];
			
			u[RK][node[boundary_node[i]].n_n[l]] = u[RK][node[boundary_node[i]].n_n[k]];
			v[RK][node[boundary_node[i]].n_n[l]] = v[RK][node[boundary_node[i]].n_n[k]];
			w[RK][node[boundary_node[i]].n_n[l]] = w[RK][node[boundary_node[i]].n_n[k]];
			p[RK][node[boundary_node[i]].n_n[l]] = p[RK][node[boundary_node[i]].n_n[k]];
			rho[RK][node[boundary_node[i]].n_n[l]] = rho[RK][node[boundary_node[i]].n_n[k]];
			t[RK][node[boundary_node[i]].n_n[l]] = t[RK][node[boundary_node[i]].n_n[k]];	
			a[RK][node[boundary_node[i]].n_n[l]] = a[RK][node[boundary_node[i]].n_n[k]];
			e[RK][node[boundary_node[i]].n_n[l]] = e[RK][node[boundary_node[i]].n_n[k]];
			//mu[RK][node[boundary_node[i]].n_n[l]] = mu[RK][boundary_node[i]];
			
			u[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = u[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			v[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = v[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			w[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = w[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			p[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = p[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			rho[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = rho[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			t[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = t[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];	
			a[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = a[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			e[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = e[RK][node[node[boundary_node[i]].n_n[k]].n_n[k]];
			//mu[RK][node[node[boundary_node[i]].n_n[l]].n_n[l]] = mu[RK][boundary_node[i]];
	
			u[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];	
			a[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[node[boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[RK][node[node[node[boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][boundary_node[i]];
		}
		else if (node[boundary_node[i]].loc == 7 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);

		}
		else if (node[boundary_node[i]].loc == 8 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 9 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 10 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		
		else if (node[boundary_node[i]].loc == 11 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 12 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 13 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);			
		}
		else if (node[boundary_node[i]].loc == 14 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[node[boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);

		}
		else if (node[boundary_node[i]].loc == 15 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[5]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 16 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[1]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 17 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[4]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 18 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[3]].n_n[0]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		
		else if (node[boundary_node[i]].loc == 19 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[5]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 20 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[1]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 21 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[4]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 22 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[3]].n_n[2]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 23 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[1]].n_n[5]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 24 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[1]].n_n[4]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 25 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[3]].n_n[4]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}
		else if (node[boundary_node[i]].loc == 26 )
		{
			u[RK][boundary_node[i]] = u[RK][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			v[RK][boundary_node[i]] = v[RK][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			w[RK][boundary_node[i]] = w[RK][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			p[RK][boundary_node[i]] = p[RK][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			t[RK][boundary_node[i]] = t[RK][node[node[boundary_node[i]].n_n[3]].n_n[5]];
			rho[RK][boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][boundary_node[i]]/t[RK][boundary_node[i]]);
			e[RK][boundary_node[i]] = p[RK][boundary_node[i]]/(0.4*rho[RK][boundary_node[i]]);
		}		
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< wal_node; i++)
	{	
		if (node[wall_node[i]].loc > 0 && node[wall_node[i]].loc <= 6)
		{	
			if (node[wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[wall_node[i]].n_n[k]];
			//p[RK][wall_node[i]]= 1.0/(1.4*Mach*Mach);
			t[RK][wall_node[i]]= t[RK][node[wall_node[i]].n_n[k]];	
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);		
			a[RK][wall_node[i]] = sqrt(1.4*p[RK][wall_node[i]]/rho[RK][wall_node[i]]);	
		
			u[RK][node[wall_node[i]].n_n[l]] = (-1.0)*u[RK][node[wall_node[i]].n_n[k]];
			v[RK][node[wall_node[i]].n_n[l]] = (-1.0)*v[RK][node[wall_node[i]].n_n[k]];
			w[RK][node[wall_node[i]].n_n[l]] = (-1.0)*w[RK][node[wall_node[i]].n_n[k]];
			p[RK][node[wall_node[i]].n_n[l]] = p[RK][node[wall_node[i]].n_n[k]];
			rho[RK][node[wall_node[i]].n_n[l]] = rho[RK][node[wall_node[i]].n_n[k]];
			t[RK][node[wall_node[i]].n_n[l]] = t[RK][node[wall_node[i]].n_n[k]];
			a[RK][node[wall_node[i]].n_n[l]] = a[RK][node[wall_node[i]].n_n[k]];
			e[RK][node[wall_node[i]].n_n[l]] = e[RK][node[wall_node[i]].n_n[k]];
			//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
		
			u[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			v[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			w[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			p[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			rho[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			t[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			a[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			e[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[node[wall_node[i]].n_n[k]].n_n[k]];
			//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
		
			u[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			a[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[node[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
			
			/*******************************BLEED*********************************************/
		/* 	if (back_pressure != 0.0)
			{
				if ((node[wall_node[i]].x >= 1420.0 && node[wall_node[i]].x <= 1226.5+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1226.5+700.0 && node[wall_node[i]].x <= 1411.224+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1411.224+700.0 && node[wall_node[i]].x <= 1564.39+700.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1015.238+700.0 && node[wall_node[i]].x <= 1222.0+700.0 && node[wall_node[i]].loc == 1) || (node[wall_node[i]].x >= 1222.25+700.0 && node[wall_node[i]].x <= 1482.0+700.0 && node[wall_node[i]].loc == 1))
	//			if ((node[wall_node[i]].x >= 1226.5+700.0 && node[wall_node[i]].x <= 1411.224+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3) || (node[wall_node[i]].x >= 1411.224+700.0 && node[wall_node[i]].x <= 1564.39+700.0 && node[wall_node[i]].loc == 3))	
				{
					if (node[wall_node[i]].x >= 1420.0 && node[wall_node[i]].x <= 1226.5+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3)
					{
						Cd = Cdr1;
					} 
					
					if (node[wall_node[i]].x >= 1226.5+700.0 && node[wall_node[i]].x <= 1411.224+700.0 && node[wall_node[i]].y < 674.0 && node[wall_node[i]].loc == 3)
					{
						Cd = Cdr2;
					}
					
					if (node[wall_node[i]].x >= 1411.224+700.0 && node[wall_node[i]].x <= 1564.39+700.0 && node[wall_node[i]].loc == 3)
					{
						Cd = Cdr3;
					}
					
					if (node[wall_node[i]].x >= 1015.238+700.0 && node[wall_node[i]].x <= 1222.0+700.0 && node[wall_node[i]].loc == 1)
					{
						Cd = Cdc1;
					}
				
					if (node[wall_node[i]].x >= 1222.25+700.0 && node[wall_node[i]].x <= 1482.0+700.0 && node[wall_node[i]].loc == 1)
					{
						Cd = Cdc2;
					}
					
					if (node[wall_node[i]].loc == 3)
					{
						Vn_sign = -1.0;
					}
					
					if (node[wall_node[i]].loc == 1)
					{
						Vn_sign = 1.0;
					}
					
					rho[RK][wall_node[i]]= rho[RK][node[wall_node[i]].n_n[k]];
					l_angle = atan((node[node[wall_node[i]].n_n[1]].y-node[wall_node[i]].y)/(node[node[wall_node[i]].n_n[1]].x-node[wall_node[i]].x));
					l_angle2 = atan((node[node[node[wall_node[i]].n_n[k]].n_n[1]].y-node[node[wall_node[i]].n_n[k]].y)/(node[node[node[wall_node[i]].n_n[k]].n_n[1]].x-node[node[wall_node[i]].n_n[k]].x));
					l_determinant = cos(l_angle)*cos(l_angle)+sin(l_angle)*sin(l_angle);

					Vt = u[RK][node[wall_node[i]].n_n[k]]*cos(l_angle)+v[RK][node[wall_node[i]].n_n[k]]*sin(l_angle);
					Vn = Vn_sign*Cd*sqrt(1.4*p[RK][node[wall_node[i]].n_n[k]]/rho[RK][node[wall_node[i]].n_n[k]]);
					Vn2 = v[RK][node[wall_node[i]].n_n[k]]*cos(l_angle2)-u[RK][node[wall_node[i]].n_n[k]]*sin(l_angle2);
					
					p[RK][wall_node[i]] = p[RK][node[wall_node[i]].n_n[k]]-2.0*rho[RK][wall_node[i]]*Vn*(Vn-Vn2);				
					u[RK][wall_node[i]] = (cos(l_angle)/l_determinant)*Vt-(sin(l_angle)/l_determinant)*Vn;
					v[RK][wall_node[i]] = (sin(l_angle)/l_determinant)*Vt+(cos(l_angle)/l_determinant)*Vn;
					w[RK][wall_node[i]]= 0.0;
					t[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/rho[RK][wall_node[i]]);
				//	rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
					e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);		
					a[RK][wall_node[i]] = sqrt(1.4*p[RK][wall_node[i]]/rho[RK][wall_node[i]]);	
					
					u[RK][node[wall_node[i]].n_n[l]] = u[RK][wall_node[i]];
					v[RK][node[wall_node[i]].n_n[l]] = v[RK][wall_node[i]];
					w[RK][node[wall_node[i]].n_n[l]] = w[RK][wall_node[i]];
					p[RK][node[wall_node[i]].n_n[l]] = p[RK][wall_node[i]];
					rho[RK][node[wall_node[i]].n_n[l]] = rho[RK][wall_node[i]];
					t[RK][node[wall_node[i]].n_n[l]] = t[RK][wall_node[i]];	
					a[RK][node[wall_node[i]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[wall_node[i]].n_n[l]] = e[RK][wall_node[i]];
					//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
					
					u[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = u[RK][wall_node[i]];
					v[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = v[RK][wall_node[i]];
					w[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = w[RK][wall_node[i]];
					p[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[RK][wall_node[i]];
					rho[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[RK][wall_node[i]];
					t[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[RK][wall_node[i]];	
					a[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[RK][wall_node[i]];
					//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
			
					u[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][wall_node[i]];
					v[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][wall_node[i]];
					w[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][wall_node[i]];
					p[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][wall_node[i]];
					rho[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][wall_node[i]];
					t[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][wall_node[i]];	
					a[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][wall_node[i]];
					//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
					
				} 		
			} */
		}
		else if (node[wall_node[i]].loc == 7 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);

		}
		else if (node[wall_node[i]].loc == 8 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 9 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 10 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		
		else if (node[wall_node[i]].loc == 11 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 12 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 13 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);			
		}
		else if (node[wall_node[i]].loc == 14 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[node[wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 15 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[wall_node[i]].n_n[5]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[5]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 16 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[wall_node[i]].n_n[1]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[1]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 17 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]]= p[RK][node[node[wall_node[i]].n_n[4]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[4]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 18 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[3]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[3]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		
		else if (node[wall_node[i]].loc == 19 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[5]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[5]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 20 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[1]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[1]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 21 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[4]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[4]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 22 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[3]].n_n[2]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[3]].n_n[2]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 23 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[1]].n_n[5]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[1]].n_n[5]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 24 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[1]].n_n[4]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[1]].n_n[4]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 25 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[3]].n_n[4]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[3]].n_n[4]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 26 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[node[wall_node[i]].n_n[3]].n_n[5]];
			t[RK][wall_node[i]] = t[RK][node[node[wall_node[i]].n_n[3]].n_n[5]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 29 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[singular[wall_node[i]].n_n[3]].n_n[4]];
			t[RK][wall_node[i]] = t[RK][node[singular[wall_node[i]].n_n[3]].n_n[4]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		else if (node[wall_node[i]].loc == 30 )
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][node[singular[wall_node[i]].n_n[3]].n_n[5]];
			t[RK][wall_node[i]] = t[RK][node[singular[wall_node[i]].n_n[3]].n_n[5]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);
		}
		
		if (node[wall_node[i]].ID == 10 && node[wall_node[i]].loc != 29 && node[wall_node[i]].loc != 30)
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][singular[wall_node[i]].n_n[3]];
			t[RK][wall_node[i]] = t[RK][singular[wall_node[i]].n_n[3]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
				
				if (h<=1)
				{
					u[RK][node[wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[wall_node[i]].n_n[k]];
					v[RK][node[wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[wall_node[i]].n_n[k]];
					w[RK][node[wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[wall_node[i]].n_n[k]];
					p[RK][node[wall_node[i]].n_n[l]] = p[RK][singular[wall_node[i]].n_n[k]];
					rho[RK][node[wall_node[i]].n_n[l]] = rho[RK][singular[wall_node[i]].n_n[k]];
					t[RK][node[wall_node[i]].n_n[l]] = t[RK][singular[wall_node[i]].n_n[k]];	
					a[RK][node[wall_node[i]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[wall_node[i]].n_n[l]] = e[RK][singular[wall_node[i]].n_n[k]];
					//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
				
					u[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					v[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					w[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					p[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					rho[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					t[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];	
					a[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
		
					u[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];	
				}
				if (h == 2)
				{
					u[RK][node[wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[wall_node[i]].n_n[k]];
					v[RK][node[wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[wall_node[i]].n_n[k]];
					w[RK][node[wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[wall_node[i]].n_n[k]];
					p[RK][node[wall_node[i]].n_n[l]] = p[RK][singular[wall_node[i]].n_n[k]];
					rho[RK][node[wall_node[i]].n_n[l]] = rho[RK][singular[wall_node[i]].n_n[k]];
					t[RK][node[wall_node[i]].n_n[l]] = t[RK][singular[wall_node[i]].n_n[k]];	
					a[RK][node[wall_node[i]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[wall_node[i]].n_n[l]] = e[RK][singular[wall_node[i]].n_n[k]];
					//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
				
					u[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					v[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					w[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					p[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					rho[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					t[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];	
					a[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
					//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
		
					u[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
					e[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];	
				}
			}
		}
		
		if (node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 56)
		{
			u[RK][wall_node[i]]= 0.0;
			v[RK][wall_node[i]]= 0.0;
			w[RK][wall_node[i]]= 0.0;
			p[RK][wall_node[i]] = p[RK][singular[wall_node[i]].n_n[0]];
			t[RK][wall_node[i]] = t[RK][singular[wall_node[i]].n_n[0]];
			rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
			e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);		
			
			if ((node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 36) || (node[wall_node[i]].loc >= 41 && node[wall_node[i]].loc <= 44) || (node[wall_node[i]].loc >= 49 && node[wall_node[i]].loc <= 52))
			{
				d_temp6 = 2;
			}
			if ((node[wall_node[i]].loc >= 37 && node[wall_node[i]].loc <= 40) || (node[wall_node[i]].loc >= 45 && node[wall_node[i]].loc <= 48))
			{
				d_temp6 = 3;
			}
			if (node[wall_node[i]].loc >= 53 && node[wall_node[i]].loc <= 56)
			{
				d_temp6 = 1;
			}
		
			for (d_temp = 0; d_temp<d_temp6; d_temp++)
			{
				if (d_temp == 0 && ((node[wall_node[i]].loc >= 33 && node[wall_node[i]].loc <= 36) || (node[wall_node[i]].loc >= 37 && node[wall_node[i]].loc <= 40)))
				{
					l = 2; 
					k = 0;
				}
				if (d_temp == 0 && ((node[wall_node[i]].loc >= 41 && node[wall_node[i]].loc <= 44) || (node[wall_node[i]].loc >= 45 && node[wall_node[i]].loc <= 48)))
				{
					l = 0; 
					k = 2;
				}
				if (d_temp == 0 && (node[wall_node[i]].loc == 49 || node[wall_node[i]].loc == 52))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 0 && (node[wall_node[i]].loc == 50 || node[wall_node[i]].loc == 51))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 33 || node[wall_node[i]].loc == 43 || node[wall_node[i]].loc == 51 || node[wall_node[i]].loc == 52 || node[wall_node[i]].loc == 39 || node[wall_node[i]].loc == 40 || node[wall_node[i]].loc == 47 || node[wall_node[i]].loc == 48))
				{
					l = 3; 
					k = 1;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 35 || node[wall_node[i]].loc == 41 || node[wall_node[i]].loc == 49 || node[wall_node[i]].loc == 50 || node[wall_node[i]].loc == 37 || node[wall_node[i]].loc == 38 || node[wall_node[i]].loc == 45 || node[wall_node[i]].loc == 46  ))
				{
					l = 1; 
					k = 3;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 34 || node[wall_node[i]].loc == 44))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 1 && (node[wall_node[i]].loc == 36 || node[wall_node[i]].loc == 42))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 2 && (node[wall_node[i]].loc == 37 || node[wall_node[i]].loc == 40 || node[wall_node[i]].loc == 48 || node[wall_node[i]].loc == 45))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 2 && (node[wall_node[i]].loc == 38 || node[wall_node[i]].loc == 39 || node[wall_node[i]].loc == 46 || node[wall_node[i]].loc == 47))
				{
					l = 4; 
					k = 5;
				}
				
				u[RK][wall_node[i]]= 0.0;
				v[RK][wall_node[i]]= 0.0;
				w[RK][wall_node[i]]= 0.0;
				p[RK][wall_node[i]] = p[RK][singular[wall_node[i]].n_n[k]];
				t[RK][wall_node[i]] = t[RK][singular[wall_node[i]].n_n[k]];
				rho[RK][wall_node[i]] = (1.4*Mach*Mach)*(p[RK][wall_node[i]]/t[RK][wall_node[i]]);
				e[RK][wall_node[i]] = p[RK][wall_node[i]]/(0.4*rho[RK][wall_node[i]]);	
				
				
				u[RK][node[wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[wall_node[i]].n_n[k]];
				v[RK][node[wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[wall_node[i]].n_n[k]];
				w[RK][node[wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[wall_node[i]].n_n[k]];
				p[RK][node[wall_node[i]].n_n[l]] = p[RK][singular[wall_node[i]].n_n[k]];
				rho[RK][node[wall_node[i]].n_n[l]] = rho[RK][singular[wall_node[i]].n_n[k]];
				t[RK][node[wall_node[i]].n_n[l]] = t[RK][singular[wall_node[i]].n_n[k]];	
				a[RK][node[wall_node[i]].n_n[l]] = a[RK][wall_node[i]];
				e[RK][node[wall_node[i]].n_n[l]] = e[RK][singular[wall_node[i]].n_n[k]];
				//mu[RK][node[wall_node[i]].n_n[l]] = mu[RK][wall_node[i]];
			
				u[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				v[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				w[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				p[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				rho[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				t[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];	
				a[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
				e[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[wall_node[i]].n_n[k]].n_n[k]];
				//mu[RK][node[node[wall_node[i]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];
	
				u[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				v[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				w[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				p[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				rho[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				t[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
				a[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][wall_node[i]];
				e[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				//mu[RK][node[node[node[wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][wall_node[i]];	
				
			}
		}	
		
	}
	
	/****************************************************sd_outlet********************************************************/		
	for (i=0; i<sd_out_node; i++)
	{	
		if (node[sd_outlet_node[i]].loc > 0 && node[sd_outlet_node[i]].loc <= 6)
		{	
			if (node[sd_outlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_outlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_outlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_outlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_outlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_outlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[RK][sd_outlet_node[i]]= u[RK][node[sd_outlet_node[i]].n_n[k]];	
			v[RK][sd_outlet_node[i]]= v[RK][node[sd_outlet_node[i]].n_n[k]];	
			w[RK][sd_outlet_node[i]]= w[RK][node[sd_outlet_node[i]].n_n[k]];	
			p[RK][sd_outlet_node[i]]= p[RK][node[sd_outlet_node[i]].n_n[k]];	
		/* 	if (back_pressure != 0.0)
			{
				if (node[sd_outlet_node[i]].y < 412.0 && node[sd_outlet_node[i]].y > 767.0)
				{
					p[RK][sd_outlet_node[i]]= p[RK][node[sd_outlet_node[i]].n_n[k]];	
				}
				if (node[sd_outlet_node[i]].y > 412.0 && node[sd_outlet_node[i]].y < 767.0)
				{
					p[RK][sd_outlet_node[i]]= back_pressure*(1.0/(1.4*Mach*Mach));	
				}
			}	 */	
			t[RK][sd_outlet_node[i]]= t[RK][node[sd_outlet_node[i]].n_n[k]];	
			rho[RK][sd_outlet_node[i]]= (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
			a[RK][sd_outlet_node[i]] = sqrt(1.4*p[RK][sd_outlet_node[i]]/rho[RK][sd_outlet_node[i]]); 
			
			u[RK][node[sd_outlet_node[i]].n_n[l]] = u[RK][sd_outlet_node[i]];
			v[RK][node[sd_outlet_node[i]].n_n[l]] = v[RK][sd_outlet_node[i]];
			w[RK][node[sd_outlet_node[i]].n_n[l]] = w[RK][sd_outlet_node[i]];
			p[RK][node[sd_outlet_node[i]].n_n[l]] = p[RK][sd_outlet_node[i]];
			rho[RK][node[sd_outlet_node[i]].n_n[l]] = rho[RK][sd_outlet_node[i]];
			t[RK][node[sd_outlet_node[i]].n_n[l]] = t[RK][sd_outlet_node[i]];	
			a[RK][node[sd_outlet_node[i]].n_n[l]] = a[RK][sd_outlet_node[i]];
			e[RK][node[sd_outlet_node[i]].n_n[l]] = e[RK][sd_outlet_node[i]];
			//mu[RK][node[sd_outlet_node[i]].n_n[l]] = mu[RK][sd_outlet_node[i]];
			
			u[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = u[RK][sd_outlet_node[i]];
			v[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = v[RK][sd_outlet_node[i]];
			w[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = w[RK][sd_outlet_node[i]];
			p[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = p[RK][sd_outlet_node[i]];
			rho[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = rho[RK][sd_outlet_node[i]];
			t[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = t[RK][sd_outlet_node[i]];	
			a[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = a[RK][sd_outlet_node[i]];
			e[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = e[RK][sd_outlet_node[i]];
			//mu[RK][node[node[sd_outlet_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_outlet_node[i]];
	
			u[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][sd_outlet_node[i]];
			v[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][sd_outlet_node[i]];
			w[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][sd_outlet_node[i]];
			p[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][sd_outlet_node[i]];
			rho[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][sd_outlet_node[i]];
			t[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][sd_outlet_node[i]];	
			a[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_outlet_node[i]];
			e[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][sd_outlet_node[i]];
			//mu[RK][node[node[node[sd_outlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_outlet_node[i]];
		}
		else if (node[sd_outlet_node[i]].loc == 7 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);

		}
		else if (node[sd_outlet_node[i]].loc == 8 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 9 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 10 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		
		else if (node[sd_outlet_node[i]].loc == 11 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 12 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 13 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);			
		}
		else if (node[sd_outlet_node[i]].loc == 14 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[node[sd_outlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);

		}
		else if (node[sd_outlet_node[i]].loc == 15 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 16 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 17 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 18 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[0]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		
		else if (node[sd_outlet_node[i]].loc == 19 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[5]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 20 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 21 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[4]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 22 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[2]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 23 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[5]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 24 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[1]].n_n[4]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 25 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[4]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
		else if (node[sd_outlet_node[i]].loc == 26 )
		{
			u[RK][sd_outlet_node[i]] = u[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			v[RK][sd_outlet_node[i]] = v[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			w[RK][sd_outlet_node[i]] = w[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			p[RK][sd_outlet_node[i]] = p[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			t[RK][sd_outlet_node[i]] = t[RK][node[node[sd_outlet_node[i]].n_n[3]].n_n[5]];
			rho[RK][sd_outlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_outlet_node[i]]/t[RK][sd_outlet_node[i]]);
			e[RK][sd_outlet_node[i]] = p[RK][sd_outlet_node[i]]/(0.4*rho[RK][sd_outlet_node[i]]);
		}
	}

/****************************************************sd_inlet********************************************************/
	for (i=0; i< sd_inl_node; i++)
	{	
		if (node[sd_inlet_node[i]].loc > 0 && node[sd_inlet_node[i]].loc <= 6)
		{	
			if (node[sd_inlet_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_inlet_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_inlet_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_inlet_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_inlet_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_inlet_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}
		
			u[RK][sd_inlet_node[i]]= 1.0;		
			v[RK][sd_inlet_node[i]]= 0.0;	
			w[RK][sd_inlet_node[i]]= 0.0;				
			p[RK][sd_inlet_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[RK][sd_inlet_node[i]]= 1.0;		
			rho[RK][sd_inlet_node[i]]= (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
			a[RK][sd_inlet_node[i]] = sqrt(1.4*p[RK][sd_inlet_node[i]]/rho[RK][sd_inlet_node[i]]); 

			u[RK][node[sd_inlet_node[i]].n_n[l]] = u[RK][sd_inlet_node[i]];
			v[RK][node[sd_inlet_node[i]].n_n[l]] = v[RK][sd_inlet_node[i]];
			w[RK][node[sd_inlet_node[i]].n_n[l]] = w[RK][sd_inlet_node[i]];
			p[RK][node[sd_inlet_node[i]].n_n[l]] = p[RK][sd_inlet_node[i]];
			rho[RK][node[sd_inlet_node[i]].n_n[l]] = rho[RK][sd_inlet_node[i]];
			t[RK][node[sd_inlet_node[i]].n_n[l]] = t[RK][sd_inlet_node[i]];	
			a[RK][node[sd_inlet_node[i]].n_n[l]] = a[RK][sd_inlet_node[i]];
			e[RK][node[sd_inlet_node[i]].n_n[l]] = e[RK][sd_inlet_node[i]];
			//mu[RK][node[sd_inlet_node[i]].n_n[l]] = mu[RK][sd_inlet_node[i]];
			
			u[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = u[RK][sd_inlet_node[i]];
			v[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = v[RK][sd_inlet_node[i]];
			w[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = w[RK][sd_inlet_node[i]];
			p[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = p[RK][sd_inlet_node[i]];
			rho[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = rho[RK][sd_inlet_node[i]];
			t[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = t[RK][sd_inlet_node[i]];	
			a[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = a[RK][sd_inlet_node[i]];
			e[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = e[RK][sd_inlet_node[i]];
			//mu[RK][node[node[sd_inlet_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_inlet_node[i]];
	
			u[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][sd_inlet_node[i]];
			v[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][sd_inlet_node[i]];
			w[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][sd_inlet_node[i]];
			p[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][sd_inlet_node[i]];
			rho[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][sd_inlet_node[i]];
			t[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][sd_inlet_node[i]];	
			a[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_inlet_node[i]];
			e[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][sd_inlet_node[i]];
			//mu[RK][node[node[node[sd_inlet_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_inlet_node[i]];
		}
		else if (node[sd_inlet_node[i]].loc == 7 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 8 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 9 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 10 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		
		else if (node[sd_inlet_node[i]].loc == 11 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 12 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 13 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);			
		}
		else if (node[sd_inlet_node[i]].loc == 14 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[node[sd_inlet_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);

		}
		else if (node[sd_inlet_node[i]].loc == 15 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 16 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 17 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 18 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[0]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		
		else if (node[sd_inlet_node[i]].loc == 19 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[5]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 20 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 21 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[4]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 22 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[2]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 23 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[5]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 24 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[1]].n_n[4]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 25 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[4]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
		else if (node[sd_inlet_node[i]].loc == 26 )
		{
			u[RK][sd_inlet_node[i]] = u[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			v[RK][sd_inlet_node[i]] = v[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			w[RK][sd_inlet_node[i]] = w[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			p[RK][sd_inlet_node[i]] = p[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			t[RK][sd_inlet_node[i]] = t[RK][node[node[sd_inlet_node[i]].n_n[3]].n_n[5]];
			rho[RK][sd_inlet_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_inlet_node[i]]/t[RK][sd_inlet_node[i]]);
			e[RK][sd_inlet_node[i]] = p[RK][sd_inlet_node[i]]/(0.4*rho[RK][sd_inlet_node[i]]);
		}
	}
	
	/****************************************************sd_boundary********************************************************/
	for (i=0; i< sd_bou_node; i++)
	{	
		if (node[sd_boundary_node[i]].loc > 0 && node[sd_boundary_node[i]].loc <= 6)
		{	
			if (node[sd_boundary_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_boundary_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_boundary_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_boundary_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_boundary_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_boundary_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

/*			u[RK][sd_boundary_node[i]]= 1.0;		
			v[RK][sd_boundary_node[i]]= 0.0;	
			w[RK][sd_boundary_node[i]]= 0.0;
			p[RK][sd_boundary_node[i]]= 1.0/(1.4*Mach*Mach);		
			t[RK][sd_boundary_node[i]]= 1.0;		
			rho[RK][sd_boundary_node[i]]= (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
			a[RK][sd_boundary_node[i]] = sqrt(1.4*p[RK][sd_boundary_node[i]]/rho[RK][sd_boundary_node[i]]); 
*/
			u[RK][sd_boundary_node[i]] = u[RK][node[sd_boundary_node[i]].n_n[k]];
			v[RK][sd_boundary_node[i]] = v[RK][node[sd_boundary_node[i]].n_n[k]];
			w[RK][sd_boundary_node[i]] = w[RK][node[sd_boundary_node[i]].n_n[k]];
			p[RK][sd_boundary_node[i]] = p[RK][node[sd_boundary_node[i]].n_n[k]];
			rho[RK][sd_boundary_node[i]] = rho[RK][node[sd_boundary_node[i]].n_n[k]];
			t[RK][sd_boundary_node[i]] = t[RK][node[sd_boundary_node[i]].n_n[k]];	
			a[RK][sd_boundary_node[i]] = a[RK][node[sd_boundary_node[i]].n_n[k]];
			e[RK][sd_boundary_node[i]] = e[RK][node[sd_boundary_node[i]].n_n[k]];
			//mu[RK][sd_boundary_node[i]] = mu[RK][sd_boundary_node[i]];
			
			u[RK][node[sd_boundary_node[i]].n_n[l]] = u[RK][node[sd_boundary_node[i]].n_n[k]];
			v[RK][node[sd_boundary_node[i]].n_n[l]] = v[RK][node[sd_boundary_node[i]].n_n[k]];
			w[RK][node[sd_boundary_node[i]].n_n[l]] = w[RK][node[sd_boundary_node[i]].n_n[k]];
			p[RK][node[sd_boundary_node[i]].n_n[l]] = p[RK][node[sd_boundary_node[i]].n_n[k]];
			rho[RK][node[sd_boundary_node[i]].n_n[l]] = rho[RK][node[sd_boundary_node[i]].n_n[k]];
			t[RK][node[sd_boundary_node[i]].n_n[l]] = t[RK][node[sd_boundary_node[i]].n_n[k]];	
			a[RK][node[sd_boundary_node[i]].n_n[l]] = a[RK][node[sd_boundary_node[i]].n_n[k]];
			e[RK][node[sd_boundary_node[i]].n_n[l]] = e[RK][node[sd_boundary_node[i]].n_n[k]];
			//mu[RK][node[sd_boundary_node[i]].n_n[l]] = mu[RK][sd_boundary_node[i]];
			
			u[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = u[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			v[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = v[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			w[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = w[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			p[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = p[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			rho[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = rho[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			t[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = t[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];	
			a[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = a[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			e[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = e[RK][node[node[sd_boundary_node[i]].n_n[k]].n_n[k]];
			//mu[RK][node[node[sd_boundary_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_boundary_node[i]];
	
			u[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];	
			a[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[node[sd_boundary_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_boundary_node[i]];
			//mu[RK][node[node[node[sd_boundary_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_boundary_node[i]];
		}
		else if (node[sd_boundary_node[i]].loc == 7 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);

		}
		else if (node[sd_boundary_node[i]].loc == 8 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 9 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 10 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		
		else if (node[sd_boundary_node[i]].loc == 11 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 12 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 13 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);			
		}
		else if (node[sd_boundary_node[i]].loc == 14 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[node[sd_boundary_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);

		}
		else if (node[sd_boundary_node[i]].loc == 15 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 16 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 17 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 18 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[0]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		
		else if (node[sd_boundary_node[i]].loc == 19 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[5]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 20 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 21 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[4]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 22 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[2]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 23 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[5]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 24 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[1]].n_n[4]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 25 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[4]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}
		else if (node[sd_boundary_node[i]].loc == 26 )
		{
			u[RK][sd_boundary_node[i]] = u[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			v[RK][sd_boundary_node[i]] = v[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			w[RK][sd_boundary_node[i]] = w[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			p[RK][sd_boundary_node[i]] = p[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			t[RK][sd_boundary_node[i]] = t[RK][node[node[sd_boundary_node[i]].n_n[3]].n_n[5]];
			rho[RK][sd_boundary_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_boundary_node[i]]/t[RK][sd_boundary_node[i]]);
			e[RK][sd_boundary_node[i]] = p[RK][sd_boundary_node[i]]/(0.4*rho[RK][sd_boundary_node[i]]);
		}		
	}
	
	/*****************************************************************************************************************************************************************/
	
	for (i=0; i< sd_wal_node; i++)
	{	
		if (node[sd_wall_node[i]].loc > 0 && node[sd_wall_node[i]].loc <= 6)
		{	
			if (node[sd_wall_node[i]].loc == 6)
			{
				l = 4; /********no element on WEST**************/
				k = 5;											
			}
			
			if (node[sd_wall_node[i]].loc == 5)
			{
				l = 5; /********no element on WEST**************/
				k = 4;											
			}
			
			if (node[sd_wall_node[i]].loc == 4)
			{
				l = 3; /********no element on WEST**************/
				k = 1;											
			}
		
			if (node[sd_wall_node[i]].loc == 3)
			{
				l = 2; /********no element on SOUTH**************/
				k = 0;							
			}
		
			if (node[sd_wall_node[i]].loc == 2)
			{
				l = 1; /********no element on EAST**************/
				k = 3;								
			}
		
			if (node[sd_wall_node[i]].loc == 1) 
			{
				l = 0; /********no element on NORTH**************/
				k = 2;							
			}

			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[sd_wall_node[i]].n_n[k]];
			//p[RK][sd_wall_node[i]]= 1.0/(1.4*Mach*Mach);
			t[RK][sd_wall_node[i]]= t[RK][node[sd_wall_node[i]].n_n[k]];	
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);		
			a[RK][sd_wall_node[i]] = sqrt(1.4*p[RK][sd_wall_node[i]]/rho[RK][sd_wall_node[i]]);	
		
			u[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[RK][node[sd_wall_node[i]].n_n[k]];
			v[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[RK][node[sd_wall_node[i]].n_n[k]];
			w[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[RK][node[sd_wall_node[i]].n_n[k]];
			p[RK][node[sd_wall_node[i]].n_n[l]] = p[RK][node[sd_wall_node[i]].n_n[k]];
			rho[RK][node[sd_wall_node[i]].n_n[l]] = rho[RK][node[sd_wall_node[i]].n_n[k]];
			t[RK][node[sd_wall_node[i]].n_n[l]] = t[RK][node[sd_wall_node[i]].n_n[k]];
			a[RK][node[sd_wall_node[i]].n_n[l]] = a[RK][node[sd_wall_node[i]].n_n[k]];
			e[RK][node[sd_wall_node[i]].n_n[l]] = e[RK][node[sd_wall_node[i]].n_n[k]];
			//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
		
			u[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			v[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			w[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			p[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			rho[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			t[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			a[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			e[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[node[sd_wall_node[i]].n_n[k]].n_n[k]];
			//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
		
			u[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			v[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			w[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			p[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			rho[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			t[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			a[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			e[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[node[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
			//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
			
			
			/*******************************BLEED*********************************************/
		/* 	if (back_pressure != 0.0)
			{
				if ((node[sd_wall_node[i]].x >= 1420.0 && node[sd_wall_node[i]].x <= 1226.5+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1226.5+700.0 && node[sd_wall_node[i]].x <= 1411.224+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1411.224+700.0 && node[sd_wall_node[i]].x <= 1564.39+700.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1015.238+700.0 && node[sd_wall_node[i]].x <= 1222.0+700.0 && node[sd_wall_node[i]].loc == 1) || (node[sd_wall_node[i]].x >= 1222.25+700.0 && node[sd_wall_node[i]].x <= 1482.0+700.0 && node[sd_wall_node[i]].loc == 1))
	//			if ((node[sd_wall_node[i]].x >= 1226.5+700.0 && node[sd_wall_node[i]].x <= 1411.224+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3) || (node[sd_wall_node[i]].x >= 1411.224+700.0 && node[sd_wall_node[i]].x <= 1564.39+700.0 && node[sd_wall_node[i]].loc == 3))	
				{
					if (node[sd_wall_node[i]].x >= 1420.0 && node[sd_wall_node[i]].x <= 1226.5+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3)
					{
						Cd = Cdr1;
					} 
					
					if (node[sd_wall_node[i]].x >= 1226.5+700.0 && node[sd_wall_node[i]].x <= 1411.224+700.0 && node[sd_wall_node[i]].y < 674.0 && node[sd_wall_node[i]].loc == 3)
					{
						Cd = Cdr2;
					}
					
					if (node[sd_wall_node[i]].x >= 1411.224+700.0 && node[sd_wall_node[i]].x <= 1564.39+700.0 && node[sd_wall_node[i]].loc == 3)
					{
						Cd = Cdr3;
					}
					
					 if (node[sd_wall_node[i]].x >= 1015.238+700.0 && node[sd_wall_node[i]].x <= 1222.0+700.0 && node[sd_wall_node[i]].loc == 1)
					{
						Cd = Cdc1;
					} 
					
					if (node[sd_wall_node[i]].x >= 1222.25+700.0 && node[sd_wall_node[i]].x <= 1482.0+700.0 && node[sd_wall_node[i]].loc == 1)
					{
						Cd = Cdc2;
					}
					 
					if (node[sd_wall_node[i]].loc == 3)
					{
						Vn_sign = -1.0;
					}
					
					if (node[sd_wall_node[i]].loc == 1)
					{
						Vn_sign = 1.0;
					}
					
					
					
					rho[RK][sd_wall_node[i]]= rho[RK][node[sd_wall_node[i]].n_n[k]];
					l_angle = atan((node[node[sd_wall_node[i]].n_n[1]].y-node[sd_wall_node[i]].y)/(node[node[sd_wall_node[i]].n_n[1]].x-node[sd_wall_node[i]].x));
					l_angle2 = atan((node[node[node[sd_wall_node[i]].n_n[k]].n_n[1]].y-node[node[sd_wall_node[i]].n_n[k]].y)/(node[node[node[sd_wall_node[i]].n_n[k]].n_n[1]].x-node[node[sd_wall_node[i]].n_n[k]].x));
					l_determinant = cos(l_angle)*cos(l_angle)+sin(l_angle)*sin(l_angle);

					Vt = u[RK][node[sd_wall_node[i]].n_n[k]]*cos(l_angle)+v[RK][node[sd_wall_node[i]].n_n[k]]*sin(l_angle);
					Vn = Vn_sign*Cd*sqrt(1.4*p[RK][node[sd_wall_node[i]].n_n[k]]/rho[RK][node[sd_wall_node[i]].n_n[k]]);
					Vn2 = v[RK][node[sd_wall_node[i]].n_n[k]]*cos(l_angle2)-u[RK][node[sd_wall_node[i]].n_n[k]]*sin(l_angle2);
					
					p[RK][sd_wall_node[i]] = p[RK][node[sd_wall_node[i]].n_n[k]]-2.0*rho[RK][sd_wall_node[i]]*Vn*(Vn-Vn2);				
					u[RK][sd_wall_node[i]] = (cos(l_angle)/l_determinant)*Vt-(sin(l_angle)/l_determinant)*Vn;
					v[RK][sd_wall_node[i]] = (sin(l_angle)/l_determinant)*Vt+(cos(l_angle)/l_determinant)*Vn;
					w[RK][sd_wall_node[i]]= 0.0;
					t[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/rho[RK][sd_wall_node[i]]);
				//	rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
					e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);		
					a[RK][sd_wall_node[i]] = sqrt(1.4*p[RK][sd_wall_node[i]]/rho[RK][sd_wall_node[i]]);	
					
					u[RK][node[sd_wall_node[i]].n_n[l]] = u[RK][sd_wall_node[i]];
					v[RK][node[sd_wall_node[i]].n_n[l]] = v[RK][sd_wall_node[i]];
					w[RK][node[sd_wall_node[i]].n_n[l]] = w[RK][sd_wall_node[i]];
					p[RK][node[sd_wall_node[i]].n_n[l]] = p[RK][sd_wall_node[i]];
					rho[RK][node[sd_wall_node[i]].n_n[l]] = rho[RK][sd_wall_node[i]];
					t[RK][node[sd_wall_node[i]].n_n[l]] = t[RK][sd_wall_node[i]];	
					a[RK][node[sd_wall_node[i]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[sd_wall_node[i]].n_n[l]] = e[RK][sd_wall_node[i]];
					//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
					
					u[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = u[RK][sd_wall_node[i]];
					v[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = v[RK][sd_wall_node[i]];
					w[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = w[RK][sd_wall_node[i]];
					p[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[RK][sd_wall_node[i]];
					rho[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[RK][sd_wall_node[i]];
					t[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[RK][sd_wall_node[i]];	
					a[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[RK][sd_wall_node[i]];
					//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
			
					u[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = u[RK][sd_wall_node[i]];
					v[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = v[RK][sd_wall_node[i]];
					w[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = w[RK][sd_wall_node[i]];
					p[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][sd_wall_node[i]];
					rho[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][sd_wall_node[i]];
					t[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][sd_wall_node[i]];	
					a[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][sd_wall_node[i]];
					//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
					
				} 
			} */			
		}
		else if (node[sd_wall_node[i]].loc == 7 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);

		}
		else if (node[sd_wall_node[i]].loc == 8 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 9 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 10 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		
		else if (node[sd_wall_node[i]].loc == 11 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[5]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 12 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[1]].n_n[4]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 13 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[4]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);			
		}
		else if (node[sd_wall_node[i]].loc == 14 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[node[sd_wall_node[i]].n_n[3]].n_n[5]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 15 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[sd_wall_node[i]].n_n[5]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[5]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 16 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 17 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]]= p[RK][node[node[sd_wall_node[i]].n_n[4]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[4]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 18 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		
		else if (node[sd_wall_node[i]].loc == 19 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[5]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[5]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 20 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 21 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[4]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[4]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 22 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[2]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[2]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 23 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[5]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[5]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 24 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[4]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[1]].n_n[4]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 25 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 26 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 29 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[4]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		else if (node[sd_wall_node[i]].loc == 30 )
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			t[RK][sd_wall_node[i]] = t[RK][node[node[sd_wall_node[i]].n_n[3]].n_n[5]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
		}
		
		
		if (node[sd_wall_node[i]].ID == 10 && node[sd_wall_node[i]].loc != 29 && node[sd_wall_node[i]].loc != 30)
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][node[sd_wall_node[i]].n_n[3]];
			t[RK][sd_wall_node[i]] = t[RK][node[sd_wall_node[i]].n_n[3]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);		
			
			for (h=0;h<3;h++)
			{
				if (h==0)
				{
					l = 2; /********no element on SOUTH**************/
					k = 0;
				}
				if (h==1) 
				{
					l = 0; /********no element on NORTH**************/
					k = 2;
				}
				if (h==2)
				{
					l = 1; /********no element on EAST**************/
					k = 3;										
				}
		
				if (h<=1)
				{
					u[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[sd_wall_node[i]].n_n[k]];
					v[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[sd_wall_node[i]].n_n[k]];
					w[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[sd_wall_node[i]].n_n[k]];
					p[RK][node[sd_wall_node[i]].n_n[l]] = p[RK][singular[sd_wall_node[i]].n_n[k]];
					rho[RK][node[sd_wall_node[i]].n_n[l]] = rho[RK][singular[sd_wall_node[i]].n_n[k]];
					t[RK][node[sd_wall_node[i]].n_n[l]] = t[RK][singular[sd_wall_node[i]].n_n[k]];	
					a[RK][node[sd_wall_node[i]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[sd_wall_node[i]].n_n[l]] = e[RK][singular[sd_wall_node[i]].n_n[k]];
					//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
							
					u[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					v[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					w[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					p[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					rho[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					t[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];	
					a[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
		
					u[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];	
				}
				if (h == 2)
				{
					u[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[sd_wall_node[i]].n_n[k]];
					v[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[sd_wall_node[i]].n_n[k]];
					w[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[sd_wall_node[i]].n_n[k]];
					p[RK][node[sd_wall_node[i]].n_n[l]] = p[RK][singular[sd_wall_node[i]].n_n[k]];
					rho[RK][node[sd_wall_node[i]].n_n[l]] = rho[RK][singular[sd_wall_node[i]].n_n[k]];
					t[RK][node[sd_wall_node[i]].n_n[l]] = t[RK][singular[sd_wall_node[i]].n_n[k]];	
					a[RK][node[sd_wall_node[i]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[sd_wall_node[i]].n_n[l]] = e[RK][singular[sd_wall_node[i]].n_n[k]];
					//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
				
					u[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					v[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					w[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					p[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					rho[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					t[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];	
					a[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
					//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
		
					u[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					v[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					w[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					p[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					rho[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					t[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
					a[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
					e[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
					//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];	
				}
			}
		}
		
		if (node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 56)
		{
			u[RK][sd_wall_node[i]]= 0.0;
			v[RK][sd_wall_node[i]]= 0.0;
			w[RK][sd_wall_node[i]]= 0.0;
			p[RK][sd_wall_node[i]] = p[RK][singular[sd_wall_node[i]].n_n[0]];
			t[RK][sd_wall_node[i]] = t[RK][singular[sd_wall_node[i]].n_n[0]];
			rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
			e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);		
			
			if ((node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 36) || (node[sd_wall_node[i]].loc >= 41 && node[sd_wall_node[i]].loc <= 44) || (node[sd_wall_node[i]].loc >= 49 && node[sd_wall_node[i]].loc <= 52))
			{
				d_temp6 = 2;
			}
			if ((node[sd_wall_node[i]].loc >= 37 && node[sd_wall_node[i]].loc <= 40) || (node[sd_wall_node[i]].loc >= 45 && node[sd_wall_node[i]].loc <= 48))
			{
				d_temp6 = 3;
			}
			if (node[sd_wall_node[i]].loc >= 53 && node[sd_wall_node[i]].loc <= 56)
			{
				d_temp6 = 1;
			}
		
			for (d_temp = 0; d_temp<d_temp6; d_temp++)
			{
				if (d_temp == 0 && ((node[sd_wall_node[i]].loc >= 33 && node[sd_wall_node[i]].loc <= 36) || (node[sd_wall_node[i]].loc >= 37 && node[sd_wall_node[i]].loc <= 40)))
				{
					l = 2; 
					k = 0;
				}
				if (d_temp == 0 && ((node[sd_wall_node[i]].loc >= 41 && node[sd_wall_node[i]].loc <= 44) || (node[sd_wall_node[i]].loc >= 45 && node[sd_wall_node[i]].loc <= 48)))
				{
					l = 0; 
					k = 2;
				}
				if (d_temp == 0 && (node[sd_wall_node[i]].loc == 49 || node[sd_wall_node[i]].loc == 52))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 0 && (node[sd_wall_node[i]].loc == 50 || node[sd_wall_node[i]].loc == 51))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 33 || node[sd_wall_node[i]].loc == 43 || node[sd_wall_node[i]].loc == 51 || node[sd_wall_node[i]].loc == 52 || node[sd_wall_node[i]].loc == 39 || node[sd_wall_node[i]].loc == 40 || node[sd_wall_node[i]].loc == 47 || node[sd_wall_node[i]].loc == 48))
				{
					l = 3; 
					k = 1;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 35 || node[sd_wall_node[i]].loc == 41 || node[sd_wall_node[i]].loc == 49 || node[sd_wall_node[i]].loc == 50 || node[sd_wall_node[i]].loc == 37 || node[sd_wall_node[i]].loc == 38 || node[sd_wall_node[i]].loc == 45 || node[sd_wall_node[i]].loc == 46  ))
				{
					l = 1; 
					k = 3;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 34 || node[sd_wall_node[i]].loc == 44))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 1 && (node[sd_wall_node[i]].loc == 36 || node[sd_wall_node[i]].loc == 42))
				{
					l = 4; 
					k = 5;
				}
				if (d_temp == 2 && (node[sd_wall_node[i]].loc == 37 || node[sd_wall_node[i]].loc == 40 || node[sd_wall_node[i]].loc == 48 || node[sd_wall_node[i]].loc == 45))
				{
					l = 5; 
					k = 4;
				}
				if (d_temp == 2 && (node[sd_wall_node[i]].loc == 38 || node[sd_wall_node[i]].loc == 39 || node[sd_wall_node[i]].loc == 46 || node[sd_wall_node[i]].loc == 47))
				{
					l = 4; 
					k = 5;
				}
				u[RK][sd_wall_node[i]]= 0.0;
				v[RK][sd_wall_node[i]]= 0.0;
				w[RK][sd_wall_node[i]]= 0.0;
				p[RK][sd_wall_node[i]] = p[RK][singular[sd_wall_node[i]].n_n[k]];
				t[RK][sd_wall_node[i]] = t[RK][singular[sd_wall_node[i]].n_n[k]];
				rho[RK][sd_wall_node[i]] = (1.4*Mach*Mach)*(p[RK][sd_wall_node[i]]/t[RK][sd_wall_node[i]]);
				e[RK][sd_wall_node[i]] = p[RK][sd_wall_node[i]]/(0.4*rho[RK][sd_wall_node[i]]);
				
				u[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*u[RK][singular[sd_wall_node[i]].n_n[k]];
				v[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*v[RK][singular[sd_wall_node[i]].n_n[k]];
				w[RK][node[sd_wall_node[i]].n_n[l]] = (-1.0)*w[RK][singular[sd_wall_node[i]].n_n[k]];
				p[RK][node[sd_wall_node[i]].n_n[l]] = p[RK][singular[sd_wall_node[i]].n_n[k]];
				rho[RK][node[sd_wall_node[i]].n_n[l]] = rho[RK][singular[sd_wall_node[i]].n_n[k]];
				t[RK][node[sd_wall_node[i]].n_n[l]] = t[RK][singular[sd_wall_node[i]].n_n[k]];	
				a[RK][node[sd_wall_node[i]].n_n[l]] = a[RK][sd_wall_node[i]];
				e[RK][node[sd_wall_node[i]].n_n[l]] = e[RK][singular[sd_wall_node[i]].n_n[k]];
				//mu[RK][node[sd_wall_node[i]].n_n[l]] = mu[RK][sd_wall_node[i]];
			
				u[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				v[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				w[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				p[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = p[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				rho[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = rho[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				t[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = t[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];	
				a[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
				e[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = e[RK][node[singular[sd_wall_node[i]].n_n[k]].n_n[k]];
				//mu[RK][node[node[sd_wall_node[i]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];
	
				u[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*u[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				v[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*v[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				w[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = (-1.0)*w[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				p[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = p[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				rho[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = rho[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				t[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = t[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];	
				a[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = a[RK][sd_wall_node[i]];
				e[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = e[RK][node[node[singular[sd_wall_node[i]].n_n[k]].n_n[k]].n_n[k]];
				//mu[RK][node[node[node[sd_wall_node[i]].n_n[l]].n_n[l]].n_n[l]] = mu[RK][sd_wall_node[i]];	
				
			}
		}
		
		
	}

}

