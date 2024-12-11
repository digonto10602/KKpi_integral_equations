#!/bin/bash

#rm printer.o

g++ -g -o ope_printer ./source/ope_plotting.cpp -O3 -std=c++14 -I/usr/include/eigen3/ -fopenmp 

