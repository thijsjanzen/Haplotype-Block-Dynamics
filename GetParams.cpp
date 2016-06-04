//
//  GetParams.cpp
//  Adapted from Hanno Hildenbrandt, 2007
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2016 Thijs Janzen. All rights reserved.
//



#include "GetParams.h"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

GetParams::GetParams() {
    seed = 1;
    refit = 0;
    popSize = 1000;
    genomeSize = 1000;
    initRatio = 0.5;
    maxTime = 201;
    meanBlockSize = -10;

    replicates = 1;

    selfing = 0;
}

void GetParams::readFromIni(const char * filename ) {
    // locate file and tranfer text to stringstream
    std::ifstream ifs(filename);
    std::stringstream ss;

    // only for succesfully created ifstream: otherwise null-pointer?
    if (ifs) {
        // config.ini content is transferred to stringstream: easier to search?
        ss << ifs.rdbuf();
    } else {
        throw "Can't locate file";
    }

    while (ss.good()) {
        readNameValuePair(ss,  "seed", seed);
        readNameValuePair(ss,  "refit", refit);
        readNameValuePair(ss,  "popSize", popSize);
        readNameValuePair(ss,  "genomeSize", genomeSize);
        readNameValuePair(ss,  "initRatio", initRatio);
        readNameValuePair(ss,  "maxTime", maxTime);
        readNameValuePair(ss,  "recomDist", recomDist);
        readNameValuePair(ss,  "replicates", replicates);
        readNameValuePair(ss,  "selfing", selfing);
        readNameValuePair(ss,  "numRecombinations", numberRecombinations);
    }
}

template <typename T>
void GetParams::readNameValuePair(std::stringstream& ss,
                                  std::string iniName, T& value ) {
    std::string name;
    char sign;

    // >> copy ss content to string until white space is encountered
    ss >> name;
    if (name != iniName )
        throw "expect parameter";
    ss >> sign;  // copies ss content to character
    if (sign != '=' )
        throw "text format of ini file is not compatible";
    ss >> value;
    std::cout << iniName << ": " << value << std::endl;
}
