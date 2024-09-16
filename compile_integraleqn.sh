#!/bin/bash

#rm printer.o

g++ -g -o printer ./source/printer_test.cpp -O3 -std=c++14 -I/usr/include/eigen3/ -fopenmp 

