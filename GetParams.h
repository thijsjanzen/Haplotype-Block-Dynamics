//
//  GetParams.h
//  Adapted from Hanno Hildenbrandt
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//

#ifndef FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_
#define FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_

#include <sstream>
#include <string>
#include <vector>


class GetParams {
 public:
    GetParams();
    void readFromIni(const char * filename);

    template <typename T>
    void readNameValuePair(std::stringstream& ss,
                           std::string iniName, T& value);

    // GLOBAL VARIABLES
    int seed;
    int refit;
    int popSize;
    int genomeSize;
    double initRatio;
    int maxTime;

    double meanBlockSize;

    int recomDist;

    int replicates;

    int selfing;

    double numberRecombinations;
};

#endif  // FINITE_CHROMOSOME_ALWAYS_RECOM_GETPARAMS_H_
