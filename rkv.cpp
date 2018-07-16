/*
 * rkv.cpp
 *
 *  Created on: 11.07.2018
 *      Author: ls
 */
#include "rkv.h"

#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

/* function prototypes */
double dnx(double, double [], double [], int );
double rk4_dn1(double(*)(double, double [], double [], int ),
               double, double, double [], double [], int);















/*==============================================================
  System of first order differential equations for the RK solver

  For a system of n first-order ODEs
  x [] array - x values
  dx[] array - dx/dt values

  For a system of n/2 second order ODEs follow the agreement
  In:  x[] array
  # first n/2 elements are x
  # last  n/2 elements are dx/dt
  Out: dx[] array
  # first n/2 elements are first order derivatives  (or dx/dt)
  # last  n/2 elements are second order derivatives (or d2x/dt2)
  example: 2D projectile motion in (x,y) plane
  In           Out
  x[0] - x     dx[0] - x'
  x[1] - y     dx[1] - y'
  x[2] - x'    dx[2] - x"
  x[3] - y'    dx[3] - y"
==============================================================*/


// Muss noch geändert werden!!!!!!!!!!:

 double dnx(double t, double r[], double dr[], int n)
{
/* first order */
    dr[0] = r[2];
    dr[1] = r[3];
/* second order */
    //dr[2] = r[2]*t+(-1.0)*(G*M/(pow(r[0],2)+pow(r[1],2)))*cos((atan(r[1]/r[0])*rad)+(1/2)*(c_w/m)*A*rho*pow(r[2]*cos((atan(r[1]/r[0]))*rad), 2));				//zweite Ableitungen (r und phi)
    //dr[3] = r[3]*t+(-1.0)*(G*M/(pow(r[0],2)+pow(r[1],2)))*sin((atan(r[1]/r[0])*rad)+(1/2)*(c_w/m)*A*rho*pow(r[3]*sin((atan(r[1]/r[0]))*rad), 2));

    dr[2] = ((-1.0)*G*M/(pow(r[0],2)+pow(r[1],2)))*cos((atan(r[1]/r[0])*rad))-(1/2)*(c_w/m)*A*rho*pow(r[2]*cos((atan(r[1]/r[0]))*rad), 2);				//zweite Ableitungen (r und phi)
    dr[3] = ((-1.0)*G*M/(pow(r[0],2)+pow(r[1],2)))*sin((atan(r[1]/r[0])*rad))-(1/2)*(c_w/m)*A*rho*pow(r[3]*sin((atan(r[1]/r[0]))*rad), 2);

    //printf("%f ",atan(r[1]/r[0])*rad);

    return 0;
}



/*==========================================================
 rk4_dn1.cpp: Solution of a system of n first-order ODE
 method:      Runge-Kutta 4th-order
 written by: Alex Godunov
 last revision: 7 October 2009
----------------------------------------------------------
 call ...
 dnx(t,x[],dx[],n)- functions dx/dt   (supplied by a user)
 input ...
 ti    - initial time
 tf    - solution time
 xi[]  - initial values
 n     - number of first order equations
 output ...
 xf[]  - solutions
==========================================================*/
double rk4_dn1(double(dnx)(double, double [], double [], int),
               double ti, double tf, double ri[], double rf[], int n)
{
      double h, t, r[n], dr[n];
      double k1[n],k2[n],k3[n],k4[n];
      int j;

      h = tf-ti;
      t = ti;
//k1
      dnx(t, ri, dr, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k1[j] = h*dr[j];
          r[j]  = ri[j] + k1[j]/2.0;
        }
//k2
      dnx(t+h/2.0, r, dr, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k2[j] = h*dr[j];
          r[j]  = ri[j] + k2[j]/2.0;
        }
//k3
      dnx(t+h/2.0, r, dr, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k3[j] = h*dr[j];
          r[j]  = ri[j] + k3[j];
        }
//k4 and result
      dnx(t+h, r, dr, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k4[j] = h*dr[j];
          rf[j] = ri[j] + k1[j]/6.0+k2[j]/3.0+k3[j]/3.0+k4[j]/6.0;
        }
    return 0.0;

    // system ("gnuplot Plot.gp");
}




