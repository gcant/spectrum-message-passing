CC=g++
CFLAGS=-O3 -fopenmp -march=native -std=c++11 -I include/eigen3 

ifeq ($(shell uname -s), Darwin)
	CC=clang++
	CFLAGS=-Xpreprocessor -fopenmp -lomp -O3 -march=native -std=c++11 -I include/eigen3
endif

spectrum_main:
	$(CC) spectrum_main.cpp -o spectrum $(CFLAGS)

single_thread:
	$(CC) spectrum_main.cpp -o spectrum -O3 -march=native -std=c++11 -I include/eigen3 
