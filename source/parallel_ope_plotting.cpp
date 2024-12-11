/* The purpose of this script is to run the python script to plot 
the ope files in parallel */

/* Executable of this code should be run from within the folder where
we have the ope data files, this code needs to inputs after the declaration 
of the exe, 1) number of data files to plot, 2) destination of the python script 
e.g. ./ope_plot 100 ../../source/ope_surface_plot.py 
this will then iterate through 0..100 for the files and execute the plotting 
scripts in parallel */

/* THIS CODE HAS BEEN INTEGRATED INSIDE OPE_PLOTTING.CPP
FILE, NO NEED FOR THIS CODE ANYMORE */

#include<bits/stdc++.h>
#include<omp.h>

int main(int argc, char* argv[])
{
    char* endptr;  
    int number_of_files = std::strtol(argv[1], &endptr, 10);
    std::string python_script = argv[2];
    int omp_ind = 0;
    //std::cout<<"running from "<<omp_ind<<std::endl; 
    
    #pragma omp parallel for 
    for(omp_ind=0; omp_ind<number_of_files; ++omp_ind)
    {
        std::string runcommand_str = "python  " + python_script + " " + std::to_string(omp_ind);
        
        char const *runcommand = runcommand_str.c_str();
        //std::cout<<"running from "<<omp_ind<<std::endl; 
        
        int runfile = std::system(runcommand);
    }

    return 0; 
}