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

/* global variables */

const double G = 6.67*pow(10, -11);             // Gravitationskonstante
const double m = 455000;                		// Masse ISS
const double M = 5.97*pow(10, 24);				//Masse Erde
const double c_w = 0.6;							//Widerstandsbeiwert
const double A = 2667.78;						//Querschnittsfläche ISS
const double rho = 1.77*pow(10, -22);			//Dichte der Atmosphäre (Geschätzter Wert!!!!!!!!!)
const double R = 6371000;						//Radius Erde


const double rad = 3.1415926/180.0;  // radians


int main()
{
    const int n=4;                   // number of first-order equations
    double ti, tf, dt, tmax;
    double ri[n], rf[n];		//Array
    double v0, phi0;
    double phi;
    int i;						// , key;


/* output: file and formats */
    ofstream file;
    file.open ("table01c.dat");         // Ergebnisse werden in dieses File geschrieben
/* output format */
    file.precision(3);
    file.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);




/* initial information */
    ti = 0.0;                // initial value for variable t
    v0 = 7666.67;             	 // Geschwindigkeit (m/s)
    phi0 =  0.0;             // initial angle (degrees)
    ri[0] = 408000+R;        // Anfangsposition in x (m)
    ri[1] = 0.0;             // Anfangsposition in y (m)
    ri[2] = 0.0; 			  // Geschwindigkeit in x-Richtung (m.s)
    ri[3] = v0; 				  // Geschwindigkeit in y-Richtung  (m/s)

    dt = 10;             // step size for integration (s)
    tmax = 170000;          // integrate till tmax (s)
/* end of initial information */

    file << setw(12) << "t"  << setw(12) << "x"   << setw(12) << "y"
         << setw(12) << "x'" << setw(12) << "y'"  << endl;

/* integration of ODE */
    while (ti <= tmax)
    {

     file << setw(12) << ti   << setw(12) << ri[0] << setw(12) << ri[1]
          << setw(12) << ri[2]<< setw(12) << ri[3] << endl;

     //if (sqrt(pow(ri[0], 2)+pow(ri[1],2)) < R) break;

     tf = ti + dt;
/*=================================*/
     rk4_dn1(dnx, ti, tf, ri, rf, n);
/*=================================*/

// prepare for the next step
        ti = tf;
        for (i = 0; i<=n-1; i = i+1)
        {
           ri[i] = rf[i];
        }
   }

    // system ("read");



    return 0;
}

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
    dr[2] = (-G*M/(pow(r[0],2)+pow(r[2],2)))*cos((atan(r[1]/r[0])*rad))-(1/2)*(c_w/m)*A*rho*pow(r[3]*cos((atan(r[1]/r[0]))*rad), 2);				//zweite Ableitungen (r und phi)
    dr[3] = (-G*M/(pow(r[0],2)+pow(r[2],2)))*sin((atan(r[1]/r[0])*rad))-(1/2)*(c_w/m)*A*rho*pow(r[1]*sin((atan(r[1]/r[0]))*rad), 2);

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

    system ("gnuplot Plot.gp");
}




