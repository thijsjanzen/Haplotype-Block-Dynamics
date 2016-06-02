#include "GetParams.h"
#include <fstream>
#include <iostream>

//This cpp file was adapted from
//H.Hildenbrandt 2007




GetParams::GetParams() {
	seed = 1;
	refit = 0;
	popSize = 1000;
	genomeSize = 1000;
	initRatio = 0.5;
	maxTime = 201;
	allOutput = 0;
	meanBlockSize = -10;
	
	infer = 1;
	replicates = 1;
	
	secondStep = false;
	selfing = 0;
}

void GetParams::readFromIni( const char * filename ) {
	
	//locate file and tranfer text to stringstream
	std::ifstream ifs( filename ); 
	std::stringstream ss;

	if( ifs ) { //only for succesfully created ifstream: otherwise null-pointer?
		ss << ifs.rdbuf(); //config.ini content is transferred to stringstream: easier to search?
	}
	else {
		throw "Can't locate file";
		//std::cout << "Can't locate file" << std::endl;
	}
	while( ss.good() ) {	
				readNameValuePair( ss,  "seed", seed);
				readNameValuePair( ss,  "refit", refit);
				readNameValuePair( ss,  "popSize", popSize);
				readNameValuePair( ss,  "genomeSize", genomeSize);
				readNameValuePair( ss,  "initRatio", initRatio);
				readNameValuePair( ss,  "maxTime", maxTime);
				readNameValuePair( ss,  "allOutput", allOutput);
				readNameValuePair( ss,  "recomDist",recomDist);
				readNameValuePair( ss,  "replicates",replicates);
				readNameValuePair( ss,  "selfing",selfing);
				readNameValuePair( ss,  "numRecombinations",numberRecombinations);
		
	}
}

template <typename T>
void GetParams::readNameValuePair( std::stringstream& ss, std::string iniName, T& value ) {
	std::string name;
	char sign;
	ss >> name; //>> copies ss content to string until white space is encountered
	if( name != iniName )
		throw "expect parameter";
	ss >> sign; //copies ss content to character 
	if( sign != '=' )
		throw "text format of ini file is not compatible";
    ss >> value;
	std::cout << iniName << ": " << value << std::endl;
//	std::cerr << iniName << ": " << value << std::endl;
}


