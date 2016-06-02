#ifndef GETPARAMS_H_INCLUDED
#define GETPARAMS_H_INCLUDED

#include <sstream>
#include <string>
#include <vector>


class GetParams {
public:
	GetParams();
	void readFromIni( const char * filename );
	
	template <typename T>
	void readNameValuePair( std::stringstream& ss, std::string iniName, T& value );
	


	//GLOBAL VARIABLES
	int seed;
	int refit;
	int popSize;
	int genomeSize;
	double initRatio;
	int maxTime;
	int allOutput;
	int infer;
	
	
	double meanBlockSize;
	
	bool secondStep;
	
	int recomDist;
	
	int replicates;
	
	int selfing;
	
	double numberRecombinations;
	
};

#endif //GETPARAMS_H_INCLUDED