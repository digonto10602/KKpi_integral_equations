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

def M2k_plotter_a2_m1m2(m1, m2):
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "M2k_a0_2.000000.dat"

    (sigk, ReM2k_1, ImM2k_1, ReM2k_2, ImM2k_2) = np.genfromtxt(filename1, unpack=True)

    fig, ax = plt.subplots(1, 2, figsize = (24, 8))

    threshold1 = pow((m1 + m2),2.0)
    threshold2 = pow((m1 - m2),2.0)

    ax[0].axvline(x=threshold2, color='black')
    ax[1].axvline(x=threshold1, color='black')


    ax[0].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[1].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[0].set_ylabel("$\mathcal{M}_2^{(m_1,m_2)}$")
    
    ax[0].set_xlim([0,0.03])
    ax[0].set_ylim([-100,100])
    ax[0].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[1].set_xlim([0.03,4.53])
    ax[1].set_ylim([-1000,1000])
    ax[1].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[1].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[0].legend()
    #ax[1].legend() 
    plt.tight_layout()
    output = "M2k_a2_m1m2.pdf"
    plt.savefig(output)
    #plt.show() 

def M2k_plotter_a6_m1m2(m1, m2):
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "M2k_a0_6.000000.dat"

    (sigk, ReM2k_1, ImM2k_1, ReM2k_2, ImM2k_2) = np.genfromtxt(filename1, unpack=True)

    fig, ax = plt.subplots(1, 2, figsize = (24, 8))

    threshold1 = pow((m1 + m2),2.0)
    threshold2 = pow((m1 - m2),2.0)

    ax[0].axvline(x=threshold2, color='black')
    ax[1].axvline(x=threshold1, color='black')


    ax[0].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[1].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[0].set_ylabel("$\mathcal{M}_2^{(m_1,m_2)}$")
    
    ax[0].set_xlim([0.005,0.015])
    ax[0].set_ylim([-100,100])
    ax[0].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[1].set_xlim([2,4.53])
    ax[1].set_ylim([-1000,1000])
    ax[1].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[1].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[0].legend()
    #ax[1].legend() 
    plt.tight_layout()
    output = "M2k_a6_m1m2.pdf"
    plt.savefig(output)
    #plt.show() 

def M2k_plotter_a16_m1m2(m1, m2):
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "M2k_a0_16.000000.dat"

    (sigk, ReM2k_1, ImM2k_1, ReM2k_2, ImM2k_2) = np.genfromtxt(filename1, unpack=True)

    fig, ax = plt.subplots(1, 2, figsize = (24, 8))

    threshold1 = pow((m1 + m2),2.0)
    threshold2 = pow((m1 - m2),2.0)

    ax[0].axvline(x=threshold2, color='black')
    ax[1].axvline(x=threshold1, color='black')


    ax[0].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[1].set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax[0].set_ylabel("$\mathcal{M}_2^{(m_1,m_2)}$")
    ax[0].set_xlim([0.005,0.015])
    ax[0].set_ylim([-100,100])
    ax[0].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[0].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[1].set_xlim([2,4.53])
    ax[1].set_ylim([-1000,1000])
    ax[1].plot(sigk, ReM2k_1, marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="real")
    ax[1].plot(sigk, ImM2k_1, marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="imag")
    
    ax[0].legend()
    #ax[1].legend() 

    plt.tight_layout()
    output = "M2k_a16_m1m2.pdf"
    plt.savefig(output)
    #plt.show() 

def M2k_plotter_a2_6_16_m1m1(m1, m2):
    plt.rcParams.update({'font.size': 22})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)

    filename1 = "M2k_a0_2.000000.dat"
    filename2 = "M2k_a0_6.000000.dat"
    filename3 = "M2k_a0_16.000000.dat"

    (sigk1, ReM2k_11, ImM2k_11, ReM2k_21, ImM2k_21) = np.genfromtxt(filename1, unpack=True)
    (sigk2, ReM2k_12, ImM2k_12, ReM2k_22, ImM2k_22) = np.genfromtxt(filename2, unpack=True)
    (sigk3, ReM2k_13, ImM2k_13, ReM2k_23, ImM2k_23) = np.genfromtxt(filename3, unpack=True)


    fig, ax = plt.subplots(figsize = (12, 8))

    threshold1 = pow((m1 + m2),2.0)
    threshold2 = pow((m1 - m2),2.0)

    ax.axvline(x=threshold2, color='black')
    ax.axvline(x=threshold1, color='black')


    ax.set_xlabel("$\sigma_k^{1}/m_1^2$")
    #ax.set_xlabel("$\sigma_k^{1}/m_1^2$")
    ax.set_ylabel("$\mathcal{M}_2^{(m_1,m_1)}$")
    ax.set_xlim([2.25,4.75])
    ax.set_ylim([0,2000])
    ax.plot(sigk1, np.abs(np.sqrt(ReM2k_21**2 + ImM2k_21**2)), marker='o',markerfacecolor="white", markersize=10, color="darkred",linestyle='none', markeredgewidth=2, label="$m_1 a_0 = 2$")
    ax.plot(sigk2, np.abs(np.sqrt(ReM2k_22**2 + ImM2k_22**2)), marker='o',markerfacecolor="white", markersize=10, color="teal",linestyle='none', markeredgewidth=2, label="$m_1 a_0 = 6$")
    ax.plot(sigk3, np.abs(np.sqrt(ReM2k_23**2 + ImM2k_23**2)), marker='o',markerfacecolor="white", markersize=10, color="darkorange",linestyle='none', markeredgewidth=2, label="$m_1 a_0 = 16$")
    
    ax.legend()
    

    plt.tight_layout()
    output = "M2k_a2_6_16_m1m1.pdf"
    plt.savefig(output)
    #plt.show() 


m1 = 1.0
m2 = 0.9

M2k_plotter_a2_m1m2(m1, m2)

M2k_plotter_a6_m1m2(m1, m2)

M2k_plotter_a16_m1m2(m1, m2)

M2k_plotter_a2_6_16_m1m1(m1, m1)
