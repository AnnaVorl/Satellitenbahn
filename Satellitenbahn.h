/*
 * Sinkflug.h
 *
 *  Created on: Jun 29, 2018
 *      Author: m1741169
 */

#ifndef SATELLITENBAHN_H_
#define SATELLITENBAHN_H_

#include <stdio.h>
#include <math.h>

class Sinkflug {

	double G=6.67*pow(10, -11);				//Gravitationskonstante
	double M=5.97*pow(10, 24);				//Masse der Erde
	double R=6371000;						//Erdradius
	double pi=3.14;
	double x;								//Zeit
	double y;								//Summe aus R und h
	double z;								//Flughoehe ueber der Erdoberflaeche


public:
	Sinkflug(double t, double r, double h); //Konstruktor
	Sinkflug();								//Defaultkonstruktor

	virtual ~Sinkflug();

											//setter Funktionen
	void setx (double t);					//Zeit
	void sety (double r);					//Radius der Flugbahn (r=R+h)
	void setz (double h);					//Höhe




	void Höhe();
};

#endif /* SATELLITENBAHN_H_ */
