//
//  output.h
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#ifndef __SecondaryContact__output__
#define __SecondaryContact__output__

#include <stdio.h>
#include <vector>
#include "Fish.h"
//#include "Particle.h"

//struct particle;


struct Output
{
	std::vector<double> meanFreqs;
	std::vector< std::vector<double> > Blocks;
	std::vector<int> timesBlocks;
	std::vector<double> avgBlocks;
	std::vector< std::vector<int> > numBlocks;
	std::vector< std::vector<double> > freqs;


	void calcMeanFreqs(const std::vector<Fish>& V);
	void countBlocksizes(const std::vector<Fish>& Pop, int time);
	void countBlocksizes_fast(const std::vector<Fish>& Pop, int time);
	//void writeToFile(const particle& Params,int time);
	//void writeToFileBlocks(const particle &Params, int L);
};

//void countBlocks(const std::vector<bool>& B, std::vector<double>& V, int& numBlocks);
template <typename T>
void countBlocks(const std::vector<bool>& B, std::vector<T>& V, int& numBlocks);

#endif /* defined(__SecondaryContact__output__) */
