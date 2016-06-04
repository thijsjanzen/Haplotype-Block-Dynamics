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
#include <numeric>  // std::accumulate
#include <string>

#include <unistd.h>  // for mac specific run code

#include "./GetParams.h"
#include "./Fish.h"
#include "./randomc.h"

template <typename T>
double calculateMean(const std::vector<T>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = 1.0 * sum / v.size();
    return(mean);
}

struct particle {
    int popSize;
    double ratio;
    int time;
    int genomeSize;

    particle() {   // default constructor
        popSize = 1000;
        ratio = 0.5;
        time = 200;
        genomeSize = 100;
    }

    particle(int pS, double R, double T, int G):
      popSize(pS), ratio(R), time(T), genomeSize(G)	{
    }
};

std::ostream& operator << (std::ostream& os, const particle& p) {
    os << p.popSize << "\t";
    os << p.ratio   << "\t";
    os << p.time << "\t";
    os << p.genomeSize << "\t";
    return os;
}

// count the number of blocks in a chromosome
int countBlocks(const std::vector<bool>& B) {
    int numBlocks = 0;
    bool val = B[0];
    int blockSize = 1;
    for (std::size_t i = 1; i < B.size(); ++i) {
        if (B[i] != val) {
            val = B[i];
            blockSize = 1;
            numBlocks++;
        } else {
            blockSize++;
        }
    }
    if (blockSize > 1) {numBlocks++;}
    return numBlocks;
}

// calculate the average number of blocks in the population
double calcAvgBlocks(const std::vector<Fish>& Pop) {
    double averageNumBlocks = 0;
    for (auto i = Pop.begin(); i != Pop.end(); ++i) {
        averageNumBlocks += countBlocks((*i).chromosome1);
        averageNumBlocks += countBlocks((*i).chromosome2);
    }

    // the number of blocks averaged per (diploid) individual
    averageNumBlocks = 1.0 * averageNumBlocks / (2 * Pop.size());
    return averageNumBlocks;
}


std::vector<double> doSimulation(const particle& Params,
                                 int RecomDist,
                                 int selfing,
                                 double numRecombinations) {
    std::vector< double > averageNumberOfBlocks;
    std::vector<Fish> Pop;

    // reserve memory
    averageNumberOfBlocks.reserve(Params.time);
    Pop.reserve(Params.popSize);

    // these are the two ancestral species
    Fish parent1 = Fish(0, Params.genomeSize);
    Fish parent2 = Fish(1, Params.genomeSize);

    // initialize the hybrid population
    for (int i = 0; i < Params.popSize; ++i) {
        Fish p1 = parent1;
        Fish p2 = parent1;
        if (uniform() < Params.ratio) {
            p1 = parent2;
        }
        if (uniform() < Params.ratio) {
            p2 = parent2;
        }
        Pop.push_back(mate(p1, p2, RecomDist, numRecombinations));
    }

    // now run for T generations:
    for (int t = 1; t < Params.time; ++t) {
        averageNumberOfBlocks.push_back(calcAvgBlocks(Pop));

        std::vector<Fish> newGeneration;

        for (int i = 0; i < Params.popSize; ++i) {
            int index1 = random_number(Params.popSize);
            int index2 = random_number(Params.popSize);

            if (selfing == 0) {
                while (index1 == index2) {
                    index2 = random_number(Params.popSize);
                }
            }

            Fish kid = mate(Pop[index1],
                            Pop[index2],
                            RecomDist,
                            numRecombinations);

            newGeneration.push_back(kid);
        }

        Pop = newGeneration;
    }

    averageNumberOfBlocks.push_back(calcAvgBlocks(Pop));

    return averageNumberOfBlocks;
}

bool file_exists(const std::string& name)   {
    std::ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }
}

void macstart(const char * argv[]);  // forward declaration


int main(int argc, const char * argv[]) {
    std::cout << "Simulating hybridization without secondary geneflow...\n";
    // this function checks if the code is executed on OS X,
    // and then changes the working directory to
    // the directory of the executable.
    macstart(argv);

    GetParams P;
    P.readFromIni("config.ini");
    set_seed(P.seed);

    if (P.refit == 10)  {
        particle genParams(P.popSize, P.initRatio,
                           P.maxTime, P.genomeSize);

        for (int r = 0; r < P.replicates; ++r) {
            std::cout << genParams << "\t" << r << "\n";

            std::vector<double> avgBlocks =
                doSimulation(genParams,
                             P.recomDist,
                             P.selfing,
                             P.numberRecombinations);

            std::ofstream outFile("record.txt", std::ios::app);
            outFile << genParams << "\t" << r << "\t";

            for (std::size_t t = 0; t < avgBlocks.size(); ++t) {
                outFile << avgBlocks[t] << "\t";
            }
            outFile << "\n";
            outFile.close();
        }
        return 0;
    }

    std::cout << "Done!\n";
    return 0;
}

void macstart(const char * argv[])  {
    std::cout << "\n\n\n";
#ifdef __APPLE__
    char *dirsep = strrchr(argv[0], '/');
    if ( dirsep != NULL ) *dirsep = 0;
    int changeDir = chdir(argv[0]);
    std::cout << "Changing Dir: " << changeDir << "\n";
    std::string cwd = getcwd(NULL, 0);
    std::cout << cwd << "\n";
    std::cout << "Starting simulation\n";
#endif
}
