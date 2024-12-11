#!/bin/bash

#rm printer.o

g++ -g -o ope_plot ./source/parallel_ope_plotting.cpp -O3 -std=c++14 -fopenmp 

