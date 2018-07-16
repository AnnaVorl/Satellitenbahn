/*
 * rkv.h
 *
 *  Created on: 11.07.2018
 *      Author: ls
 */

#ifndef RKV_H_
#define RKV_H_

#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>


using namespace std;				//für Ausgabe

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




#endif /* RKV_H_ */
