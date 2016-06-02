//
//  Particle.h
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#ifndef __SecondaryContact__Particle__
#define __SecondaryContact__Particle__

#include <stdio.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include "output.h"


//struct Output; //forward declaration

struct particle
{
	int popSize;
	double ratio;
	int time;
	
	
	particle()
	{
		popSize = 1000;
		ratio = 0.5;
		time = 200;

	}
	particle(int pS, double R, double T) : popSize(pS), ratio(R), time(T)	{
	}
	
	particle(const particle& other);
	particle& operator=(const particle& other);
	
	
	
	//void getFromPrior(int infer);
};

std::istream& operator >> (std::istream& is, particle& p);
std::ostream& operator << (std::ostream& os, const particle& p);

void output(const std::vector<particle>& P, int time);
void outputcout(const std::vector<particle>& P, int time);

void readParticles(int t, std::vector<particle>& parts);

double calculateMeanBlockSize(const std::vector<double>& v);
double calculateMeanBlockSize2(const Output& O);

#endif /* defined(__SecondaryContact__Particle__) */
