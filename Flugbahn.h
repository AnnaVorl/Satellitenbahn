/*
 * Flugbahn.h
 *
 *  Created on: 16.07.2018
 *      Author: ls
 */

#ifndef FLUGBAHN_H_
#define FLUGBAHN_H_

#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>



using namespace std;

//Deklaration der beiden Funktionen

double dnx(double, double [], double [], int );
double rkv4(double(*)(double, double [], double [], int ),
               double, double, double [], double [], int);


//globale Variablen

const double G = 6.67*pow(10, -11);             // Gravitationskonstante in m^3/(kg*s^2)
const double m = 455000;                		// Masse ISS in kg
const double M = 5.97*pow(10, 24);				// Masse Erde in kg
const double c_w = 0.6;							// Widerstandsbeiwert (dimensionslos)
const double A = 2667.78;						// Querschnittsfläche ISS in m^2
const double rho0 = 1.3;						// Dichte der Atmosphäre auf der Erdoberfläche in kg/m^3
const double R = 6371000;						// Radius Erde in m
const double H0 = 7945.19;						// Parameter in der barometrischen Höhenformel in m



int main()
{
    const int n=4;                   // Anzahl DGLs 1. Ordnung
    double ti, tf, dt, tmax;
    double ri[n], rf[n];
    double v0;
    int i;


// output: file and formats
    ofstream file;
    file.open ("Ergebnisse.dat");         // Ergebnisse werden in dieses File geschrieben
// output format
    file.precision(3);
    file.setf(ios::scientific | ios::showpoint);
    cout.precision(6);
    cout.setf(ios::fixed | ios::showpoint);

//Anfangssinformationen und Startwerte
    ti = 0.0;                // Anfangsinformation
    v0 = 7666.67;            // Geschwindigkeit in m/s
    ri[0] = 408000+R;        // x-Anfangsposition in m
    ri[1] = 0.0;             // y-Anfangsposition m
    ri[2] = 0.0; 			 // Geschwindigkeit in x-Richtung in m/s
    ri[3] = v0; 			 // Geschwindigkeit in y-Richtung  in m/s

    dt = 10;             	// Schrittweite für Integration in s
    tmax = 17000000;          // Integration bis tmax in s



    file << setw(12) << "t"  << setw(12) << "x"   << setw(12) << "y"
         << setw(12) << "x'" << setw(12) << "y'"  << endl;

//Integration
    while (ti <= tmax)
    {

     file << setw(12) << ti   << setw(12) << ri[0] << setw(12) << ri[1]
          << setw(12) << ri[2]<< setw(12) << ri[3] << endl;


     tf = ti + dt;
/*=================================*/
     rkv4(dnx, ti, tf, ri, rf, n);
/*=================================*/

// Vorbereitung für den nächsten Schritt
        ti = tf;
        for (i = 0; i<=n-1; i = i+1)
        {
           ri[i] = rf[i];
        }
   }

    return 0;
}




#endif /* FLUGBAHN_H_ */
