#This script is used to plot the surface of OPE
#function for G11 G12 and G21. This script will 
#be called from a c++ code that will run this in 
#parallel to generate the surface plots using 
#parallelism. 

import numpy as np 
import matplotlib as mpl 
import matplotlib.pyplot as plt 
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate
import sys #important for using argv from input

def surface_plot_1(opefile, N):
    plt.rcParams.update({'font.size': 18})
    plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    plt.rc('text', usetex=True)
    
    #The number of points in each axis is hard coded to 250
    #N = 250
    
    (re_En, im_En,
     re_qb1, im_qb1,
     re_qb2, im_qb2, 
     re_p, im_p, 
     re_G11, im_G11, 
     re_G12, im_G12, 
     re_G21, im_G21 ) = np.genfromtxt(opefile, skip_header=1, unpack=True)
    
    list_G = [[re_G11, im_G11], [re_G12, im_G12], [re_G21, im_G21]] 
    #list_G = [[1,2],[3,4],[5,6]]
    
    re_En_val = re_En[0]
    im_En_val = im_En[0]
    re_qb1_val = re_qb1[0]
    im_qb1_val = im_qb1[0]
    re_qb2_val = re_qb2[0]
    im_qb2_val = im_qb2[0] 

    xi = np.linspace(re_p.min(), re_p.max(), N)
    yi = np.linspace(im_p.min(), im_p.max(), N)
    
    Xi, Yi = np.meshgrid(xi, yi) 

    fig, ax = plt.subplots(3, 2, figsize = (24, 24))

    title_str = "En = " + str(re_En_val) 
    fig.suptitle(title_str, fontsize=20)

    Gind_list = ["(1,1)","(1,2)","(2,1)"]

    for i in range(0, len(list_G), 1):
        for j in range(0, len(list_G[0]), 1):
            Garr = list_G[i][j]
            Z = scipy.interpolate.griddata((re_p, im_p), Garr, (Xi, Yi), method='linear')
            Gind = Gind_list[i]
            if(j==0):
                #Gind = "(" + str(i+1) + "," + str(j+1) + ")" 
                sub_title_str = "Re $G_S^{" + Gind + "}$"
            elif(j==1):
                #Gind = "(" + str(i+1) + "," + str(j+1) + ")" 
                sub_title_str = "Im $G_S^{" + Gind + "}$"

            ax[i,j].set_title(sub_title_str, fontsize=20)
            ax[i,j].set_xlim([re_p.min(), re_p.max()])
            ax[i,j].set_ylim([im_p.min(), im_p.max()])
            ax[i,j].tick_params(axis='both', which='major', labelsize=20)
            ax[i,j].tick_params(axis='both', which='minor', labelsize=20)
            ax[i,j].set_xlabel("Re $p$",fontsize=20)
            ax[i,j].set_ylabel("Im $p$",fontsize=20)
            ax[i,j].contour(Xi, Yi, Z, levels=np.linspace(-200, 200, num=1000), linewidths=0.5, colors='black', linestyles='solid')
            h0 = ax[i,j].contourf(Xi, Yi, Z, levels=np.linspace(-30, 30, num=30))
            divider = make_axes_locatable(ax[i,j])
            cax = divider.append_axes("right", "5%", pad="3%")
            cbar = fig.colorbar(h0, cax=cax, orientation="vertical")
            cax.xaxis.set_ticks_position("default")
            cax.xaxis.set_label_position("bottom")
            cax.tick_params(labelsize=18)
            #x = np.linspace(0, 2 * np.pi, 100)
            #y = np.sin(x)
            #ax[i,j].plot(x,y)
            
    plt.tight_layout()
    outfilename = str(opefile).split(".")
    output = outfilename[0] + ".png"
    plt.savefig(output)
    plt.close()         

#opefile = 'ope_file_80.dat'

if len(sys.argv) > 1:
    filename = sys.argv[1] + sys.argv[2] + ".dat"
    N = int(sys.argv[3])
else: 
    print("python didn't get proper arguments")
    exit() 

opefile = filename 
surface_plot_1(opefile, N) 