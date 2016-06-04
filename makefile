all: test

test : main.cpp mersenne.cpp GetParams.cpp Fish.cpp
	g++ main.cpp mersenne.cpp GetParams.cpp Fish.cpp -Wall -O3 -o run

