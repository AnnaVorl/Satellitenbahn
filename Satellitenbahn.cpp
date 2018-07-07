/*
 * Sinkflug.cpp
 *
 *  Created on: Jun 29, 2018
 *      Author: m1741169
 */

#include "Satellitenbahn.h"

#include <stdio.h>
#include <math.h>


Sinkflug::Sinkflug() {
	// TODO Auto-generated constructor stub
	x=0;
	y=0;
	z=0;

}

Sinkflug::Sinkflug(double t, double r, double h) {
	// TODO Auto-generated destructor stub

	x=t;
	y=r;
	z=h;

}

	void Sinkflug::setx (double t) {x=t;}
	void Sinkflug::sety (double r) {y=r;}
	void Sinkflug::setz (double h) {z=h;}

	void Sinkflug::Höhe() {
		for(x=0 ; z>=0 ; x++) {
			z=pow((G*M)/(4*pi*pi*x*x), (1/3))-R;
		}
		printf ("\n %f", z);
		return;
	}

