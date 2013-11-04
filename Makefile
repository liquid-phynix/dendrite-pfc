all: main

timing.o: timing.cpp timing.hpp
	g++ -c timing.cpp -o timing.o

main.o: main.cpp distributed_matrix.hpp global_transposition.hpp distributed_fft.hpp technicality.hpp params.hpp solution.hpp contour.hpp
	mpic++ -std=c++0x -Wall -c main.cpp -o main.o

main: main.o timing.o
	mpic++ -lm -lfftw3 main.o timing.o -o main
