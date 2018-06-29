/*
 * Sinkflug.h
 *
 *  Created on: Jun 29, 2018
 *      Author: m1741169
 */

#ifndef SINKFLUG_H_
#define SINKFLUG_H_

#include <stdio.h>
#include <math.h>

class Sinkflug {

	double G=6.67*pow(10, -11);				//Gravitationskonstante
	double M=5.97*pow(10, 24);				//Masse der Erde
	double R=6371000;						//Erdradius
	double pi=3.14;
	double t;								//Zeit
	double r;								//Summe aus R und h
	double h;								//Flughoehe ueber der Erdoberflaeche





public:
	Sinkflug();
	virtual ~Sinkflug();
};

#endif /* SINKFLUG_H_ */
