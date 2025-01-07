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

def sigb_m2_a01_dependence():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "m2_dependence_on_sigb.dat"
    filename2 = "a0_1_dependence_on_sigb.dat"

    (m2, Re_sigb1_1, Im_sigb1_1, Re_sigb2_1, Im_sigb2_1, thres1) = np.genfromtxt(filename1, unpack=True)
    (a01, Re_sigb1_2, Im_sigb1_2, Re_sigb2_2, Im_sigb2_2, thres2) = np.genfromtxt(filename2, unpack=True)

    fig, ax = plt.subplots(1, 2, figsize = (24, 10))

    #ax[0].set_title("$m_2$ dependence")
    ax[0].plot(m2, Re_sigb1_1, color="darkred",linestyle='-', linewidth=7, label="Re $\sigma_b^{+}$", zorder=4)
    ax[0].plot(m2, Im_sigb1_1, color="#DAA520",linestyle='-', linewidth=7, label="Im $\sigma_b^{+}$", zorder=4)
    ax[0].plot(m2, Re_sigb2_1, color="teal",linestyle='-', linewidth=7, label="Re $\sigma_b^{-}$", zorder=3)
    ax[0].plot(m2, Im_sigb2_1, color="lightblue",linestyle='-', linewidth=7, label="Im $\sigma_b^{-}$", zorder=3)
    ax[0].plot(m2, thres1, color="grey",linestyle='-', linewidth=7, label="$\sigma_{th}$", zorder=1)
    ax[0].plot(m2, np.zeros(len(thres1)), color="black", linestyle='-', linewidth=1, zorder=5)
    ax[0].set_xlabel("$m_2$")

    ax[1].set_ylim([-2,3])

    #ax[1].set_title("$a_0$ dependence")
    ax[1].plot(a01, Re_sigb1_2, color="darkred",linestyle='-', linewidth=7, label="Re $\sigma_b^{+}$", zorder=4)
    ax[1].plot(a01, Im_sigb1_2, color="#DAA520",linestyle='-', linewidth=7, label="Im $\sigma_b^{+}$", zorder=4)
    ax[1].plot(a01, Re_sigb2_2, color="teal",linestyle='-', linewidth=7, label="Re $\sigma_b^{-}$", zorder=3)
    ax[1].plot(a01, Im_sigb2_2, color="lightblue",linestyle='-', linewidth=7, label="Im $\sigma_b^{-}$", zorder=3)
    ax[1].plot(a01, thres2, color="grey",linestyle='-', linewidth=7, label="$\sigma_{th}$", zorder=1)
    ax[1].plot(a01, np.zeros(len(thres2)), color="black", linestyle='-', linewidth=1, zorder=5)
    ax[1].set_xlabel("$a_0$")

    #ax[1,0].set_xscale('log')
    #ax[1,0].set_yscale('log')
    #ax[1,1].set_xscale('log')
    #ax[1,1].set_yscale('log')

    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    output = "m2_a01_dependence_on_sigb.pdf"
    plt.savefig(output)
    plt.show() 

def g_m2_a01_dependence():
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "g_i_residue.dat"

    (m2, Re_g1, Im_g1, Re_g2, Im_g2) = np.genfromtxt(filename1, unpack=True)

    fig, ax = plt.subplots(1, 2, figsize = (24, 10))

    #ax[0].set_title("$m_2$ dependence")
    ax[0].set_ylim([-20,20])
    ax[1].set_ylim([-20,20])
    ax[0].set_xlim([0,1])
    ax[1].set_xlim([0,1])
    ax[0].plot(m2, Re_g1, marker="o", color="darkred",markerfacecolor="white", markersize=20, linestyle='none', markeredgewidth=2, label="Re $g^{(+)}$", zorder=4)
    ax[0].plot(m2, Im_g1, marker="o", color="#DAA520",markerfacecolor="white", markersize=20, linestyle='none', markeredgewidth=2, label="Im $g^{(+)}$", zorder=4)
    ax[1].plot(m2, Re_g2, marker="o", color="teal",markerfacecolor="white", markersize=20, linestyle='none', markeredgewidth=2, label="Re $g^{(-)}$", zorder=3)
    ax[1].plot(m2, Im_g2, marker="o", color="lightblue",markerfacecolor="white", markersize=20, linestyle='none', markeredgewidth=2, label="Im $g^{(-)}$", zorder=3)
    ax[0].axvline(x=0.5, linestyle='solid',color='black')
    ax[1].axvline(x=0.5, linestyle='solid',color='black')
    ax[0].set_xlabel("$m_2$")
    ax[1].set_xlabel("$m_2$")

    
    ax[0].legend()
    ax[1].legend()
    plt.tight_layout()
    output = "g_i_residue_m2_dependence.pdf"
    plt.savefig(output)
    plt.show() 



#N_dependence_plots_1()   
#sigb_m2_a01_dependence()
g_m2_a01_dependence()