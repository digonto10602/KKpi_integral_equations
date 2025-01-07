#!/bin/bash

#rm printer.o

nvcc ./source/printer_test.cpp -I/usr/include/eigen3/ -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/math_libs/include -I/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/include -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/math_libs/lib64 -L/opt/nvidia/hpc_sdk/Linux_x86_64/24.3/cuda/lib64 -O3 -std=c++14 -lcublas -lcusolver -lcudart -lgomp -Xcompiler -fopenmp -o printer_gpu
