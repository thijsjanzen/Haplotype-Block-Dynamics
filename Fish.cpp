//
//  Fish.cpp
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include "Fish.h"
#include "output.h"
#include "randomc.h"
#include <cmath>
#include <fstream>

#include <algorithm>


int getRecomPos(int L, int recomDist)
{
	//we want to draw a number dependent on the recombination rate distribution
	
	//if the rate is uniform:
	if(recomDist == 0) {
		int index = random_number(L);
		while(index == 0 || index == L)
		{
			index = random_number(L);
		}
		return index;
	}
	
	if(recomDist == 2) //hotspots
	{
		static std::vector<int> rate = {14,23,54,62,67,79,94,104,112,123,131,175,181,186,201,204,208,219,224,232,242,261,266,276,299,310,312,316,322,331,360,368,371,372,373,374,377,380,383,396,400,412,416,431,442,447,454,455,468,478,490,491,507,530,546,572,581,584,610,621,624,625,638,646,650,655,657,659,679,694,697,698,706,707,715,717,722,732,751,754,760,763,766,780,789,795,799,800,815,818,826,846,859,885,894,906,917,940,944,976};
		
		static std::vector<int> rate2 = {13,21,27,35,45,55,61,66,86,97};
		
		if(L == 100) rate = rate2;
		
		static std::vector<double> rates;
		if(rates.empty())
		{
			int cnt = 0;
			for(int i = 0; i < L; ++i)
			{
				if(rate[cnt] == i) {
					rates.push_back(9);
					cnt++;
				} else {
					rates.push_back(1);
				}
			}
		}
		int pos = -100;
		while(pos < 0)
		{
			int index = random_number(L);
			while(index == 0 || index == L)
			{
				index = random_number(L);
			}
			
			double rate = rates[index];
			if(uniform() < rate / 9)
			{
				pos = index;
			}
		}
		//std::ofstream outFile("text.txt",std::ios::app);
		//outFile << pos << "\n";
		//outFile.close();
		return pos;
	}
	
	
	//if the rate is not uniform:
	int pos = -100;
	int centromere = (int)(0.5 * L);
	int maxDist = L - centromere;
	const static float factor = 1.5 * logf(10); //to ensure that at relative position 1, the rate is 10, and at position 0, the rate is 1: exp(log(10)) = 10 & exp(0) = 1;
	const static float maxProb = expf(factor);
	
	
	while(pos < 0)
	{
		int index = random_number(L);
		while(index == 0 || index == L)
		{
			index = random_number(L);
		}

		int distToCentromere = index - centromere; if(distToCentromere < 0) distToCentromere *= -1.0; //absolute distance
		float relPos = 1.0f * distToCentromere / maxDist; //in %
		
		float prob = expf(factor * relPos); //based on the stickleback recombination map, Roesti et. al. 2012;
		if(uniform() < (prob / maxProb))
		{
			pos = index;
		}
	}
	/*std::ofstream outFile("probs.txt",std::ios::app);
	outFile << pos << "\n";
	outFile.close();
	*/

	return pos;
}

void Recombine(std::vector<bool>& offspring, std::vector<bool> chromosome1, std::vector<bool> chromosome2, int recomDist )
{
	int recomPos = getRecomPos((int)chromosome1.size(), recomDist);
	
	offspring.insert(offspring.end(),chromosome1.begin(), chromosome1.begin() + recomPos);
	offspring.insert(offspring.end(),chromosome2.begin() + recomPos, chromosome2.end());
	return;
}

void Recombine(std::vector<bool>& offspring, std::vector<bool> chromosome1, std::vector<bool> chromosome2, int recomDist, double numberRecombinations )
{
	//int recomPos = getRecomPos((int)chromosome1.size(), recomDist);
	
	if(floor(numberRecombinations) != numberRecombinations) //we have a decimal value
	{
		double remain = numberRecombinations - floor(numberRecombinations);
		int add = 0;
		if(uniform() < remain) add = 1;
		
		numberRecombinations = floor(numberRecombinations) + add;
	}
	
	if(numberRecombinations == 0)
	{
		offspring = chromosome1;
		return;
	}
	
	
	
	
	
	std::vector<int> recomPos;
	int L = (int)chromosome1.size();
	for(int i = 0; i < numberRecombinations; ++i) {
		int pos = random_number(L);
		while(pos == 0 || pos == L)
		{
			pos = random_number(L);
		}
		recomPos.push_back(pos);
	}
	
	std::sort(recomPos.begin(), recomPos.end());
	auto last = std::unique(recomPos.begin(), recomPos.end());
	recomPos.erase(last,recomPos.end());
	
	while(recomPos.size() < numberRecombinations)
	{
		int remain = numberRecombinations - (int)recomPos.size();
		for(int i = 0; i < remain; ++i)
		{
			int pos = random_number(L);
			while(pos == 0 || pos == L)
			{
				pos = random_number(L);
			}
			recomPos.push_back(pos);
		}
		std::sort(recomPos.begin(), recomPos.end());
		auto last = std::unique(recomPos.begin(), recomPos.end());
		recomPos.erase(last,recomPos.end());
	}
	
	
	
	
	
	int order = 0;
	int start = 0;
	
	for(int i = 0; i < (int)recomPos.size(); ++i)
	{
		int end = recomPos[i];
		if(order == 0) {
			offspring.insert(offspring.end(),chromosome1.begin() + start, chromosome1.begin() + end);
			order = 1;
		} else {
			offspring.insert(offspring.end(),chromosome2.begin() + start, chromosome2.begin() + end);
			order = 0;
		}
		start = end;
	}
	
	if(order == 0) {
		offspring.insert(offspring.end(),chromosome1.begin() + start, chromosome1.end());
	} else {
		offspring.insert(offspring.end(),chromosome2.begin() + start, chromosome2.end());
	}
	
	return;
}



Fish mate(const Fish& A, const Fish& B, int recomDist, double numberRecombinations)
{

	Fish offspring;
	offspring.chromosome1.clear();
	offspring.chromosome2.clear(); //just to be sure.
	
	if(uniform() < 0.5) {
		Recombine(offspring.chromosome1, A.chromosome1, A.chromosome2, recomDist,numberRecombinations);
	} else {
		Recombine(offspring.chromosome1, A.chromosome2, A.chromosome1, recomDist,numberRecombinations);
	}
	
	if(uniform() < 0.5) {
		Recombine(offspring.chromosome2, B.chromosome1, B.chromosome2,recomDist,numberRecombinations);
	} else {
		Recombine(offspring.chromosome2, B.chromosome2, B.chromosome1,recomDist,numberRecombinations);
	}
	return offspring;
}
