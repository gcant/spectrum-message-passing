CC=g++-9
CFLAGS=-O3 -fopenmp -march=native -std=c++11 -I include/eigen3 

spectrum_main:
	$(CC) spectrum_main.cpp -o spectrum $(CFLAGS)
