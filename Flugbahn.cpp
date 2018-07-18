/*
 * Flugbahn.cpp
 *
 *  Created on: 16.07.2018
 *      Author: ls
 */

#include "Flugbahn.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace std;

//Deklaration der Funktionen
double dnx(double, double [], double [], int );
double rkv4(double(*)(double, double [], double [], int ),
               double, double, double [], double [], int);


/*==============================================================
  In           Out
  r[0] - x     dr[0] - x'
  r[1] - y     dr[1] - y'
  r[2] - x'    dr[2] - x"
  r[3] - y'    dr[3] - y"
==============================================================*/



 double dnx(double t, double r[], double dr[], int n)
{
//erste Ordnung
    dr[0] = r[2];
    dr[1] = r[3];
//zweite Ordnung
    dr[2] = (-1.0)*G*M/pow(sqrt((pow(r[0],2)+pow(r[1],2))),3)*r[0]-(1/2)*(c_w)*A*(rho0*exp(-(pow(r[0],2)+pow(r[1],2))/H0))*pow(r[2], 2);
    dr[3] = (-1.0)*G*M/pow(sqrt((pow(r[0],2)+pow(r[1],2))),3)*r[1]-(1/2)*(c_w)*A*(rho0*exp(-(pow(r[0],2)+pow(r[1],2))/H0))*pow(r[3], 2);

    return 0.0;
}



/*==========================================================
 method:      Runge-Kutta 4th-order
 written by: Alex Godunov
 last revision: 7 October 2009
----------------------------------------------------------
 call ...
 dnx(t,r[],dr[],n)- functions dx/dt
 input ...
 ti    - initial time
 tf    - solution time
 ri[]  - initial values
 n     - number of first order equations
 output ...
 rf[]  - solutions
==========================================================*/

double rkv4(double(dnx)(double, double [], double [], int),
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
//k4 und Ergebnis
      dnx(t+h, r, dr, n);
      for (j = 0; j<=n-1; j = j+1)
        {
          k4[j] = h*dr[j];
          rf[j] = ri[j] + k1[j]/6.0+k2[j]/3.0+k3[j]/3.0+k4[j]/6.0;
        }
    return 0.0;

}
