#This scripts plots functions
import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def N_dependence_plots_1():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "degenerate_mphib_vs_N_1.dat"
    filename2 = "degenerate_mphib_vs_N_2.dat"
    filename3 = "mphib_above_phibthreshold_a2_vs_N_2.dat"
    filename4 = "mphib_above_phibthreshold_a2_vs_N_1.dat"

    (En1, s1, Re_rhoMphib1, Im_rhoMphib1, temp, temp, diff1, N1) = np.genfromtxt(filename1, unpack=True)
    (En2, s2, Re_rhoMphib2, Im_rhoMphib2, temp, temp, diff2, N2) = np.genfromtxt(filename2, unpack=True)
    (En3, s3, Re_rhoMphib3, Im_rhoMphib3, diff3, N3) = np.genfromtxt(filename3, unpack=True)
    (En4, s4, Re_rhoMphib4, Im_rhoMphib4, diff4, N4) = np.genfromtxt(filename4, unpack=True)

    fig, ax = plt.subplots(2, 2, figsize = (24, 12))

    ax[0,0].set_title("degenerate")
    ax[0,0].plot(N1, Re_rhoMphib1, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,0].plot(N1, Im_rhoMphib1, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,0].plot(N2, Re_rhoMphib2, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,0].plot(N2, Im_rhoMphib2, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[0,0].set_ylabel("$|\\rho_{\\varphi b} \\mathcal{M}_{\\varphi b}|$")

    ax[1,0].plot(N1, diff1, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].plot(N2, diff2, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].set_ylabel("$\\Delta \\rho_{\\varphi b}$")
    ax[1,0].set_xlabel("$N$")

    ax[0,1].set_title("$2+1$ system")
    ax[0,1].plot(N3, Re_rhoMphib3, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,1].plot(N3, Im_rhoMphib3, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,1].plot(N4, Re_rhoMphib4, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,1].plot(N4, Im_rhoMphib4, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    

    ax[1,1].plot(N3, diff3, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,1].plot(N4, diff4, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    #ax[1,1].set_ylabel("$difference$")
    ax[1,1].set_xlabel("$N$")

    #ax[0,0].set_xscale('log')
    #ax[0,0].set_yscale('log')
    #ax[0,1].set_xscale('log')
    #ax[0,1].set_yscale('log')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')

    ax[0,0].legend()
    ax[0,1].legend()
    #ax[0,0].legend()
    #ax[0,0].legend()
    plt.tight_layout()
    output = "difference_between_degenerate_and_nondegenerate.pdf"
    plt.savefig(output)
    plt.show() 

def N_dependence_plots_2():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "degenerate_mphib_vs_N_en_onefourth_1.dat"
    filename2 = "degenerate_mphib_vs_N_en_onefourth_2.dat"
    filename3 = "mphib_above_phibthreshold_a2_vs_N_en_onefourth_1.dat"
    filename4 = "mphib_above_phibthreshold_a2_vs_N_en_onefourth_2.dat"

    (En1, s1, Re_rhoMphib1, Im_rhoMphib1, temp, temp, diff1, N1) = np.genfromtxt(filename1, unpack=True)
    (En2, s2, Re_rhoMphib2, Im_rhoMphib2, temp, temp, diff2, N2) = np.genfromtxt(filename2, unpack=True)
    (En3, s3, Re_rhoMphib3, Im_rhoMphib3, diff3, N3) = np.genfromtxt(filename3, unpack=True)
    (En4, s4, Re_rhoMphib4, Im_rhoMphib4, diff4, N4) = np.genfromtxt(filename4, unpack=True)

    fig, ax = plt.subplots(2, 2, figsize = (24, 12))

    ax[0,0].set_title("degenerate")
    ax[0,0].plot(N1, Re_rhoMphib1, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,0].plot(N1, Im_rhoMphib1, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,0].plot(N2, Re_rhoMphib2, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,0].plot(N2, Im_rhoMphib2, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    
    ax[0,0].set_ylabel("$|\\rho_{\\varphi b} \\mathcal{M}_{\\varphi b}|$")

    ax[1,0].plot(N1, diff1, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].plot(N2, diff2, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,0].set_ylabel("$\\Delta \\rho_{\\varphi b}$")
    ax[1,0].set_xlabel("$N$")

    ax[0,1].set_title("$2+1$ system")
    ax[0,1].plot(N3, Re_rhoMphib3, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0,1].plot(N3, Im_rhoMphib3, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    ax[0,1].plot(N4, Re_rhoMphib4, marker='o',markerfacecolor="white", markersize=20, color="darkred",linestyle='none', markeredgewidth=2)
    ax[0,1].plot(N4, Im_rhoMphib4, marker='o',markerfacecolor="white", markersize=20, color="teal",linestyle='none', markeredgewidth=2)
    

    ax[1,1].plot(N3, diff3, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    ax[1,1].plot(N4, diff4, marker='o',markerfacecolor="white", markersize=20, color="darkorange",linestyle='none', markeredgewidth=2)
    #ax[1,1].set_ylabel("$difference$")
    ax[1,1].set_xlabel("$N$")

    #ax[0,0].set_xscale('log')
    #ax[0,0].set_yscale('log')
    #ax[0,1].set_xscale('log')
    #ax[0,1].set_yscale('log')
    ax[1,0].set_xscale('log')
    ax[1,0].set_yscale('log')
    ax[1,1].set_xscale('log')
    ax[1,1].set_yscale('log')

    ax[0,0].legend()
    ax[0,1].legend()
    #ax[0,0].legend()
    #ax[0,0].legend()
    plt.tight_layout()
    output = "difference_between_degenerate_and_nondegenerate_en_onefourth.pdf"
    plt.savefig(output)
    plt.show() 



#N_dependence_plots_1()   
N_dependence_plots_2()   
    