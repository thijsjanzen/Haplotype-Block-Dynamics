//
//  Fish.h
//  SecondaryContact
//
//  Created by Thijs Janzen on 20/02/15.
//  Copyright (c) 2015 Thijs Janzen. All rights reserved.
//

#ifndef FINITE_CHROMOSOME_ALWAYS_RECOM_FISH_H_
#define FINITE_CHROMOSOME_ALWAYS_RECOM_FISH_H_

#include <stdio.h>
#include <vector>


struct Fish  {
  std::vector<bool> chromosome1;
  std::vector<bool> chromosome2;

  Fish() {
  }

  // constructor that sets all genome elements to "initLoc"
  Fish(const bool initLoc, const int genomeSize) {
      for ( int i = 0; i < genomeSize; ++i ) {
        chromosome1.push_back(initLoc);
        chromosome2.push_back(initLoc);
    }
  }

  // copy constructor
  Fish(const std::vector<bool>& A, const std::vector<bool>& B) {
    chromosome1 = A;
    chromosome2 = B;
  }
};

// mating function
Fish mate(const Fish& A, const Fish& B, int recomDist,
          double numberRecombinations);

#endif  // FINITE_CHROMOSOME_ALWAYS_RECOM_FISH_H_
