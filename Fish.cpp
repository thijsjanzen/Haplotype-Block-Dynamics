//
//  Fish.cpp
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#include "./Fish.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <vector>

#include "./randomc.h"

// get a recombination position in a finite chromosome,
// given L and the type of recombination distribution
int getRecomPos(int L, int recomDist) {
    // we want to draw a number dependent on the recombination rate distribution
    int pos = -100;

    // if the rate is uniform:
    if (recomDist == 0) {
        int index = random_number(L);
        // exclude the ends of the chromosome
        while (index == 0 || index == L) {
            index = random_number(L);
        }
        pos = index;
    }

    // recombination hotspots, used for figure Appendix 2
    if (recomDist == 2) {
        // these recombination spots were generated using R,
        // and are kept fixed across simulations
        // rate is for L = 1000, rate2 is for L = 100
        static std::vector<int> rate = {14, 23, 54, 62, 67, 79, 94, 104, 112, 123, 131, 175, 181, 186, 201, 204, 208, 219, 224, 232, 242, 261, 266, 276, 299, 310, 312, 316, 322, 331, 360, 368, 371, 372, 373, 374, 377, 380, 383, 396, 400, 412, 416, 431, 442, 447, 454, 455, 468, 478, 490, 491, 507, 530, 546, 572, 581, 584, 610, 621, 624, 625, 638, 646, 650, 655, 657, 659, 679, 694, 697, 698, 706, 707, 715, 717, 722, 732, 751, 754, 760, 763, 766, 780, 789, 795, 799, 800, 815, 818, 826, 846, 859, 885, 894, 906, 917, 940, 944, 976};

        static std::vector<int> rate2 = {13, 21, 27, 35, 45, 55, 61, 66, 86, 97};

        if (L == 100) rate = rate2;

        // the vector rates will hold the relative recombination rates
        static std::vector<double> rates;
        if (rates.empty()) {
            int cnt = 0;
            for (int i = 0; i < L; ++i) {
                if (rate[cnt] == i) {
                    rates.push_back(9);
                    cnt++;
                } else {
                    rates.push_back(1);
                }
            }
        }

        // pick a site proportional to it's recombination rate
        pos = -100;
        while (pos < 0) {
            int index = random_number(L);
            while (index == 0 || index == L) {
                index = random_number(L);
            }

            double rate = rates[index];
            if (uniform() < 1.0 * rate / 9) {
                pos = index;
            }
        }
    }

    // increased recombination rate
    // towards the peripheral ends (figure appendix 1)
    if (recomDist == 1) {
        int centromere = static_cast<int>(0.5 * L);
        int maxDist = L - centromere;

        // to ensure that at relative position 1,
        // the rate is 10, and at position 0, the rate is 1:
        // exp(log(10)) = 10 & exp(0) = 1;
        const static float factor = 1.0 * logf(10);
        const static float maxProb = expf(factor);

        pos = -100;
        while (pos < 0) {
            int index = random_number(L);
            // exclude the ends of the chromosome
            while (index == 0 || index == L) {
                index = random_number(L);
            }

            int distToCentromere = index - centromere;
            // absolute distance
            if (distToCentromere < 0) distToCentromere *= -1.0;
            float relPos = 1.0f * distToCentromere / maxDist;  // in %

            // based on the stickleback recombination map, Roesti et. al. 2012;
            float prob = expf(factor * relPos);
            if (uniform() < (prob / maxProb)) {
                pos = index;
            }
        }
    }
    return pos;
}


// function to recombine two parental chromosomes,
// and produce a new one ("offspring"), given a
// recombination distribution, and the number of recombinations
void Recombine(std::vector<bool>* offspring,
               std::vector<bool> chromosome1,
               std::vector<bool> chromosome2,
               int recomDist,
               double numberRecombinations )  {
    // we have a decimal value, the fractional part is interpreted as
    // a probability of one extra recombination, note that we do NOT
    // assume a poisson distributed number of recombinations
    // this has been proven inaccurate, and tends to overestimate
    // the number of meioses without recombination
    if (floor(numberRecombinations) != numberRecombinations) {
        double remain = numberRecombinations - floor(numberRecombinations);
        int add = 0;
        if (uniform() < remain) add = 1;

        numberRecombinations = floor(numberRecombinations) + add;
    }

    // if there are not recombinations, preliminary exit
    if (numberRecombinations == 0) {
        offspring = &chromosome1;
        return;
    }

    std::vector<int> recomPos;
    // store L, so we avoid repeated calls of the function .size()
    int L = static_cast<int>(chromosome1.size());

    while (recomPos.size() < numberRecombinations) {
        int pos = getRecomPos(L, recomDist);
        recomPos.push_back(pos);
        // sort them, in case they are not sorted yet
        // we need this to remove duplicates, and later
        // to apply crossover
        std::sort(recomPos.begin(), recomPos.end());
        // remove duplicate recombination sites
        auto last = std::unique(recomPos.begin(), recomPos.end());
        recomPos.erase(last, recomPos.end());
    }

    // used to track which chromosome was used
    // during the last recombination event
    int order = 0;
    int start = 0;

    for (std::size_t i = 0; i < recomPos.size(); ++i) {
        int end = recomPos[i];
        if (order == 0) {  // add the first chromosome
            offspring->insert(offspring->end(),
                             chromosome1.begin() + start,
                             chromosome1.begin() + end);
            order = 1;
        } else {   // add the second chromosome
            offspring->insert(offspring->end(),
                             chromosome2.begin() + start,
                             chromosome2.begin() + end);
            order = 0;
        }
        start = end;
    }

    // add chromosomal content after
    // the last recombination site:
    if (order == 0) {
        offspring->insert(offspring->end(),
                         chromosome1.begin() + start,
                         chromosome1.end());
    } else {
        offspring->insert(offspring->end(),
                         chromosome2.begin() + start,
                         chromosome2.end());
    }

    return;
}


// mating function
Fish mate(const Fish& A,
          const Fish& B,
          int recomDist,
          double numberRecombinations) {
    Fish offspring;
    offspring.chromosome1.clear();
    offspring.chromosome2.clear();  // just to be sure.

    // random order or in other words,
    // we randomly select 1 of 2 produced chromosomes during recombination
    if (uniform() < 0.5) {
        Recombine(&offspring.chromosome1,
                  A.chromosome1, A.chromosome2,
                  recomDist, numberRecombinations);
    } else {
        Recombine(&offspring.chromosome1,
                  A.chromosome2, A.chromosome1,
                  recomDist, numberRecombinations);
    }

    if (uniform() < 0.5) {
        Recombine(&offspring.chromosome2,
                  B.chromosome1, B.chromosome2,
                  recomDist, numberRecombinations);
    } else {
        Recombine(&offspring.chromosome2,
                  B.chromosome2, B.chromosome1,
                  recomDist, numberRecombinations);
    }
    return offspring;
}
