//
//  Particle.cpp
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include "Particle.h"
#include <boost/lexical_cast.hpp>
#include "output.h"
#include "randomc.h"

double getFreqFit(const std::vector<double>& E, const std::vector<double>& simV)
{
	double output = 0;
	std::vector<double> V = simV;
	
	std::sort(V.begin(),V.end());
	
	for(std::size_t i = 0; i < E.size(); ++i)
	{
		output += ((E[i] - V[i]) * (E[i] - V[i])); //squared distance
		//double a = E[i] - V[i];
		//if(a < 0.0) a *= -1.0;    //absolute distance
		//output += a;
	}
	
	return output;
}


double calculateMeanBlockSize(const std::vector<double>& v)
{
	double sum = 0;
	int count = 0;	
	std::vector<double> v2 = v;
	for(int i = 0; i < (int)v.size(); ++i)
	{
		for(int j = 0; j < v[i]; ++j)
		{
			sum += i;
			count++;
		}
	}
	double mean = 1.0 * sum / count;
	return mean;
}

double calculateMeanBlockSize2(const Output& O)
{
	double meanBSize = 0;
	double sumFreq = 0.0;
	for(int i = 0; i < (int)O.Blocks.back().size(); ++i)
	{
		double bSize = i;
		double freq = O.Blocks.back()[i];
		double add = bSize * freq;
		sumFreq += freq;
		
		meanBSize += add;
		//meanBSize += 1.0 * i * O.Blocks.back()[i];
	}
	return(meanBSize);
}







double getBlockFit(const std::vector<double>& E, const std::vector<double>& V)
{
	double output = 0;
	for(std::size_t i = 0; i < E.size(); ++i)
	{
		output += ((E[i] - V[i]) * (E[i] - V[i])); //squared distance
	}
	
	return output;
}



particle::particle(const particle& other)
{
	popSize = other.popSize;
	ratio = other.ratio;
	time = other.time;

}


particle& particle::operator=(const particle& other)
{
	if(this == &other) return *this;
	
	popSize = other.popSize;
	ratio = other.ratio;

	time = other.time;
	
	return *this;
}



std::istream& operator >> (std::istream& is, particle& p)
{
	is >> p.popSize;
	is >> p.ratio;
	is >> p.time;
	return is;
}

std::ostream& operator << (std::ostream& os, const particle& p)
{
	os << p.popSize << "\t";
	os << p.ratio   << "\t";
	os << p.time << "\t";
	return os;
}

void readParticles(int t, std::vector<particle>& parts)
{
	std::string tt = boost::lexical_cast<std::string>(t);
	std::string file_name = "particles_t=" + tt + ".txt";
	
	std::ifstream oldFile(file_name.c_str());
	while(!oldFile.eof())
	{
		particle temp;
		oldFile >> temp;
		parts.push_back(temp);
	}
	parts.pop_back();
	return;
}


