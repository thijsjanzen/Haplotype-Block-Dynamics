//
//  main.cpp
//  SMC_freqs
//
//  Created by Thijs Janzen on 29/01/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <numeric>

#include <boost/timer.hpp>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>

#include "GetParams.h"
#include "Particle.h"
#include "Fish.h"
#include "output.h"
#include "randomc.h"

int replicate = 0;
int checkStep = 1;

template <typename T>
double calculateMean(const std::vector<T>& v)
{
	double sum = std::accumulate(v.begin(), v.end(), 0.0);
	double mean = 1.0 * sum / v.size();
	return(mean);
}
template <typename T>
double calculateSD(const std::vector<T>& v)
{
	double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m =  sum / v.size();
	
	double accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
		accum += (d - m) * (d - m);
	});
	
	double stdev = sqrt(accum / (v.size()-1));
	return stdev;
}

std::vector<int> countJunctionTypes(const std::vector<Fish> Pop)
{
	std::vector<int> types(4,0);
	for(auto i = Pop.begin(); i != Pop.end(); ++i)
	{
		for(int j = 1; j < (int)(*i).chromosome1.size(); ++j)
		{
			int a = (*i).chromosome1[j];
			int b = (*i).chromosome1[j-1];
			if(a == 1 && b == 0) types[0]++;
			if(a == 0 && b == 1) types[1]++;
			if(a == 1 && b == 1) types[2]++;
			if(a == 0 && b == 0) types[3]++;

			a = (*i).chromosome2[j];
			b = (*i).chromosome2[j-1];
			if(a == 1 && b == 0) types[0]++;
			if(a == 0 && b == 1) types[1]++;
			if(a == 1 && b == 1) types[2]++;
			if(a == 0 && b == 0) types[3]++;
		}
	}
	return types;
}


Output doSimulation(int genSize, const particle& Params, int trackB, int RecomDist, int selfing, double numRecombinations)
{
	Output O;
	std::vector<Fish> Pop;
	
	std::vector<int> fusions(3,0);
	
	Fish parent1 = Fish(0,genSize);
	Fish parent2 = Fish(1,genSize);
	
	//test:
	Fish test = Fish(0,genSize);
	test.chromosome1 = parent1.chromosome1;
	test.chromosome2 = parent2.chromosome1;
	Fish offspringtest = mate(test,test, RecomDist, numRecombinations);
	

	
	
	for(int i = 0; i < Params.popSize; ++i)
	{
		Fish p1 = parent1;
		Fish p2 = parent1;
		if(uniform() < Params.ratio) {
			p1 = parent2;
		}
		if(uniform() < Params.ratio) {
			p2 = parent2;
		}
		
		Pop.push_back(mate(p1,p2, RecomDist, numRecombinations));
	}
	int t = 0;
	checkStep = 1;//Params.time / 20;
	O.countBlocksizes_fast(Pop,t);
	for(t = 1; t < Params.time; ++t)
	{
		if(trackB == 5)
		{
			O.countBlocksizes_fast(Pop,t);
			std::vector<int> types = countJunctionTypes(Pop);
			std::ofstream outFile("junctiontypes.txt",std::ios::app);
			outFile << t << "\t";
			for(int i = 0; i < (int)types.size(); ++i) outFile << types[i] << "\t";
			outFile << "\n";
			outFile.close();
		}
		
		if(trackB == 666)
		{
			if(t < 5)
			{
				O.countBlocksizes_fast(Pop,t);
			} else {
				if(t % 10 == 0) O.countBlocksizes_fast(Pop,t);
				
				int maxIndex = (int)O.avgBlocks.size();
				double a = O.avgBlocks[maxIndex-1];
				double b = O.avgBlocks[maxIndex-2];
				double c = O.avgBlocks[maxIndex-3];
				if(a == b && b == c)
				{
					break;
				}
			}
		}
		
		if(trackB == 777)
		{
			O.countBlocksizes_fast(Pop,t);
		}
		
		
		if(trackB == 2)
		{
			O.countBlocksizes(Pop,t);
			double stdev = calculateSD(O.numBlocks.back());
			std::ofstream outFile("fullRecord.txt",std::ios::app);
			outFile << Params << "\t" << t << "\t" << O.avgBlocks.back() << "\t" << stdev << "\n";
			outFile.close();
			std::cout << t << "\t" << O.avgBlocks.back() << "\t" << stdev << "\n";
			
			std::ofstream outFile2("allBlocks.txt",std::ios::app);
			outFile2 << Params << "\t" << genSize << "\t" << t << "\t";
			std::vector<int>::iterator start = O.numBlocks[O.numBlocks.size()-1].begin();
			std::vector<int>::iterator end =   O.numBlocks[O.numBlocks.size()-1].end();
			for(std::vector<int>::iterator i = start; i != end; ++i)
			{
				outFile2 << (*i) << "\t";
			}
			outFile2 << "\n";
			outFile2.close();
		}
		if(trackB  == 0)
		{
			if(t % checkStep == 0) O.countBlocksizes_fast(Pop,t);
		}
		
		
		


		
		std::vector<Fish> newGeneration;
		
		for(int i = 0; i < Params.popSize; ++i)
		{
			int index1 = random_number(Params.popSize);
			
			Fish kid;
			int index2 = random_number(Params.popSize);
			
			if(selfing == 0) while(index1 == index2) index2 = random_number(Params.popSize);

			
			if(trackB == 777) {
				if(t == 200) {
					if(uniform() < 0.5) {
						Fish p1 = Pop[index1];
						Fish p2 = Pop[index2];
						
						if(uniform() < 0.5) {
							p1 = parent1;
							if(uniform() < Params.ratio) {
								p1 = parent2;
							}
						} else {
							p2 = parent1;
							if(uniform() < Params.ratio) {
								p2 = parent2;
							}
						}
						
						kid = mate(p1,p2, RecomDist, numRecombinations);
					} else {
						kid = mate(Pop[index1],Pop[index2], RecomDist,numRecombinations);
					}
	
				} else {
					kid = mate(Pop[index1],Pop[index2], RecomDist,numRecombinations);
				}
			} else {
				kid = mate(Pop[index1],Pop[index2], RecomDist,numRecombinations);
			}
			
			newGeneration.push_back(kid);
		}
		
		Pop = newGeneration;
	}
	
	O.calcMeanFreqs(Pop);
	O.countBlocksizes_fast(Pop,Params.time);
	
	return O;
}

Output doSimulation_growth(int genSize, const particle& Params, int trackB, int RecomDist, int selfing, double numRecombinations)
{
	Output O;
	std::vector<Fish> Pop;
	

	Fish parent1 = Fish(0,genSize);
	Fish parent2 = Fish(1,genSize);
	
	for(int i = 0; i < Params.popSize; ++i)
	{
		Fish p1 = parent1;
		Fish p2 = parent1;
		if(uniform() < Params.ratio) {
			p1 = parent2;
		}
		if(uniform() < Params.ratio) {
			p2 = parent2;
		}
		
		Pop.push_back(mate(p1,p2, RecomDist, numRecombinations));
	}
	int t = 0;
	checkStep = 1;//Params.time / 20;

	O.countBlocksizes_fast(Pop,t);


	
	
	
	
	double popSize = 1.0 * Params.popSize;
	
	double r = 0.01;
	int K = Params.popSize * 10;
	
	
	for(t = 1; t < Params.time; ++t)
	{
		if(t % checkStep == 0) {
			O.countBlocksizes_fast(Pop,t);
			std::ofstream outFile("fullRecord.txt",std::ios::app);
			outFile << Params << "\t" << genSize << "\t" << popSize << "\t" << t << "\t" << O.avgBlocks.back() << "\n";
			outFile.close();
		}
		
		std::vector<Fish> newGeneration;
		
		if(t > 100) popSize += r * popSize * (1 - 1.0 * popSize / K);
		
		for(int i = 0; i < popSize; ++i)
		{
			int index1 = random_number((int)Pop.size());
			
			Fish kid;
			int index2 = random_number((int)Pop.size());
			
			if(selfing == 0) while(index1 == index2) index2 = random_number((int)Pop.size());
			
			kid = mate(Pop[index1],Pop[index2], RecomDist,numRecombinations);
			
			newGeneration.push_back(kid);
		}
		
		Pop = newGeneration;
	}
	
	O.calcMeanFreqs(Pop);
	O.countBlocksizes_fast(Pop,Params.time);
	
	return O;
}




bool file_exists(const std::string& name)
{
	std::ifstream f(name.c_str());
	if (f.good()) {
		f.close();
		return true;
	} else {
		f.close();
		return false;
	}
}


double calcBSize(const Output& O)
{
	double meanBSize = 0;
	int count = 0;
	for(int i = 0; i < (int)O.Blocks.back().size(); ++i)
	{
		if(O.Blocks.back()[i] > 0)
		{
			count++;
			meanBSize += O.Blocks.back()[i];
		}
	}
	//meanBSize *= 1.0 / count;
	return(meanBSize);
}



int main(int argc, const char * argv[]) {
	
	std::cout << "Simulating hybridization without secondary geneflow...\n";
	std::cout << argv[0] << "\n";
	
	char *dirsep = strrchr( argv[0], '/' );
	
	if( dirsep != NULL ) *dirsep = 0;
	
	std::cout << "Changing dir to: " << argv[0] << "\n";
	
	int changeDir = chdir(argv[0]);
	std::cout << "Changing Dir: " << changeDir << "\n";
	
	std::string cwd = getcwd(NULL,0);
	std::cout << "Working in" << cwd << "\n";
	
	std::cout << "newversion\n";
	GetParams P;
	P.readFromIni("config.ini");
	set_seed(P.seed);
	
	std::vector<double> empiricalData;
	std::vector<double> empBlocks;
	
	
	if(P.refit == 11)
	{
		std::vector<double> ratios = {0.5, 0.7, 0.9};
		for(int p = 0; p < (int)ratios.size(); ++p)
		{
			particle genParams(P.popSize,ratios[p],P.maxTime);
			
			for(int r = 0; r < P.replicates; ++r) {
				replicate = r;	
				std::cout << genParams << "\t" << r << "\n";
				Output O = doSimulation(P.genomeSize, genParams, P.allOutput, P.recomDist, P.selfing, P.numberRecombinations); //we use infer = 3, to make sure
				
				std::ofstream outFile("record.txt",std::ios::app);
				outFile << genParams << "\t" << r << "\t" << checkStep << "\t";
				
				for(int t = 0;t < (int)O.avgBlocks.size(); ++t)
				{
					outFile << O.avgBlocks[t] << "\t";
				}
				outFile << "\n";
				outFile.close();
			}
		}
		return 0;
	}
	
	if(P.refit == 10)
	{
		particle genParams(P.popSize,P.initRatio,P.maxTime);
		
		for(int r = 0; r < P.replicates; ++r) {
			replicate = r;
			std::cout << genParams << "\t" << r << "\n";
			Output O = doSimulation(P.genomeSize, genParams, P.allOutput, P.recomDist, P.selfing, P.numberRecombinations); //we use infer = 3, to make sure
			
			std::ofstream outFile("record.txt",std::ios::app);
			outFile << genParams << "\t" << r << "\t" << checkStep << "\t";
			
			for(int t = 0; t < (int)O.avgBlocks.size(); ++t)
			{
				outFile << O.avgBlocks[t] << "\t";
			}
			outFile << "\n";
			outFile.close();
		}
		return 0;

	}
	
	if(P.refit == 666)
	{
		particle genParams(P.popSize,P.initRatio,P.maxTime);
		
		for(int r = 0; r < P.replicates; ++r) {
			replicate = r;
			std::cout << genParams << "\t" << r << "\n";
			Output O = doSimulation(P.genomeSize, genParams, 666 , P.recomDist, P.selfing, P.numberRecombinations); //we use infer = 3, to make sure
			
			double K = O.avgBlocks.back();
			double R1 = (O.avgBlocks[1] / K) - (O.avgBlocks[0] / K);
			double R2 = (O.avgBlocks[2] / K) - (O.avgBlocks[1] / K);
			double R3 = (O.avgBlocks[3] / K) - (O.avgBlocks[2] / K);
			std::ofstream outFile("minimal.txt",std::ios::app);
			outFile << genParams << "\t" << P.genomeSize << "\t" << r << "\t" << K << "\t" << R1 << "\t" << R2 << "\t" << R3 << "\n";
			outFile.close();
		}
		return 0;

		
		
	}
	
	
	if(P.refit == 12)
	{
		//std::vector<int> N = {25,50,75,100,150,200,250,500,750,1000,2000,5000,10000};
		std::vector<int> N;
		for(int i = 140; i < 300; i+=10)
		//for(int i = 75; i < 100; i+=5)
			N.push_back(i);
		
		for(int n = 0; n < (int)N.size(); n++)
		{
			particle genParams(N[n],0.5,P.maxTime);
			
			for(int r = 0; r < P.replicates; ++r) {
				replicate = r;
				std::cout << genParams << "\t" << r << "\n";
				Output O = doSimulation(P.genomeSize, genParams, P.allOutput, P.recomDist, P.selfing,P.numberRecombinations); //we use infer = 3, to make sure
				std::ofstream outFile("record.txt",std::ios::app);
				outFile << genParams << "\t" << r << "\t" << checkStep << "\t";
				
				for(int t = 0; t < (int)O.avgBlocks.size(); ++t)
				{
					outFile << O.avgBlocks[t] << "\t";
				}
				outFile << "\n";
				outFile.close();

			}
				
		}
		return 0;
	}
	
	if(P.refit == 100)
	{
		particle genParams(P.popSize,P.initRatio,P.maxTime);
		
		for(int r = 0; r < P.replicates; ++r) {
			replicate = r;
			std::cout << genParams << "\t" << r << "\n";
			Output O = doSimulation_growth(P.genomeSize, genParams, P.allOutput, P.recomDist, P.selfing, P.numberRecombinations); //we use infer = 3, to make sure
			
			std::ofstream outFile("record.txt",std::ios::app);
			outFile << genParams << "\t" << r << "\t" << checkStep << "\t";
			
			for(int t = 0; t < (int)O.avgBlocks.size(); ++t)
			{
				outFile << O.avgBlocks[t] << "\t";
			}
			outFile << "\n";
			outFile.close();
		}
		return 0;
		
	}

	
	
	
	std::cout << "Done!\n";
	
	
	return 0;
}



