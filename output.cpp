//
//  output.cpp
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include "output.h"
#include <fstream>
#include <iostream>
#include <boost/lexical_cast.hpp>

template <typename T>
void countBlocks(const std::vector<bool>& B, std::vector<T>& V, int& numBlocks)
{
	numBlocks = 0;
	bool val = B[0];
	int blockSize = 1;
	for(std::size_t i = 1; i < B.size(); ++i)
	{
		if(B[i] != val) {
			V[blockSize]++;
			val = B[i];
			blockSize = 1;
			numBlocks++;
		} else {
			blockSize++;
		}
	}
	if(blockSize > 1) { V[blockSize]++; numBlocks++;}
	return;
}



int countBlocks_fast(const std::vector<bool>& B)
{
	int numBlocks = 0;
	bool val = B[0];
	int blockSize = 1;
	for(std::size_t i = 1; i < B.size(); ++i)
	{
		if(B[i] != val) {
			val = B[i];
			blockSize = 1;
			numBlocks++;
		} else {
			blockSize++;
		}
	}
	if(blockSize > 1) {numBlocks++;}
	return numBlocks;
}






void Output::countBlocksizes(const std::vector<Fish>& Pop, int time)
{
	std::vector<double> output(Pop[0].chromosome1.size()+1,0.0);
	std::vector<int> trackNumBlocks;
	double averageNumBlocks = 0;
	for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i)
	{
		int numB = 0;
		countBlocks((*i).chromosome1,output,numB);
		averageNumBlocks += numB;
		trackNumBlocks.push_back(numB);
		numB = 0;
		countBlocks((*i).chromosome2,output,numB);
		averageNumBlocks += numB;
		trackNumBlocks.push_back(numB);
	}
	double mult = 1.0 / averageNumBlocks; //we need this one to later normalize our values in "output"
	
	averageNumBlocks = 1.0 * averageNumBlocks / (2 * Pop.size()); //the number of blocks averaged per individual
	avgBlocks.push_back(averageNumBlocks);
	numBlocks.push_back(trackNumBlocks);
	timesBlocks.push_back(time);
	
	//and now we normalize our output vector, such that each entry represents the probability of a block being that size.
	for(std::vector<double>::iterator it = output.begin(); it != output.end(); ++it)
	{
		(*it) *= mult;
	}
	
	Blocks.push_back(output);
	return;
}

void Output::countBlocksizes_fast(const std::vector<Fish>& Pop, int time)
{
	double averageNumBlocks = 0;
	for(std::vector<Fish>::const_iterator i = Pop.begin(); i != Pop.end(); ++i)
	{
		averageNumBlocks += countBlocks_fast((*i).chromosome1);
		averageNumBlocks += countBlocks_fast((*i).chromosome2);
	}
	
	averageNumBlocks = 1.0 * averageNumBlocks / (2 * Pop.size()); //the number of blocks averaged per individual
	avgBlocks.push_back(averageNumBlocks);
	return;
}



void Output::calcMeanFreqs(const std::vector<Fish>& V)
{
	if(V.empty())
	{
		std::cout << "ERROR!!! Can't calculate Frequencies, population empty!\n";
		return;
	}
	
	std::vector<double> output(V[0].chromosome1.size(),0.0);
	for(std::vector<Fish>::const_iterator i = V.begin(); i != V.end(); ++i)
	{
		for(std::size_t j = 0; j < (*i).chromosome1.size(); ++j)
		{
			output[j] += (*i).chromosome1[j] + (*i).chromosome2[j];
		}
	}
	
	for(std::vector<double>::iterator i = output.begin(); i != output.end(); ++i) {
		(*i) = (*i) * 1.0 / (2*V.size());  //divide by 2 times popsize, because diploid
	}
	
	meanFreqs = output;
}




