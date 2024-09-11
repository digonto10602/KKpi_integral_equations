import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os.path
from scipy.optimize import curve_fit
import scipy.interpolate

def sorter_func():
    filename1 = "bs_virtual_N_500.dat"
    filename2 = "bs_virtual_N_500_1.dat"
    filename3 = "bs_virtual_N_500_2.dat"
    filename4 = "bs_virtual_N_500_3.dat"
    filename5 = "bs_virtual_N_500_4.dat"
    filename6 = "bs_virtual_N_500_5.dat"

    (a1, s1, E1, N1, phibth1) = np.genfromtxt(filename1,unpack=True)    
    (a2, s2, E2, N2, phibth2) = np.genfromtxt(filename2,unpack=True)    
    (a3, s3, E3, N3, phibth3) = np.genfromtxt(filename3,unpack=True)    
    (a4, s4, E4, N4, phibth4) = np.genfromtxt(filename4,unpack=True)    
    (a5, s5, E5, N5, phibth5) = np.genfromtxt(filename5,unpack=True)    
    (a6, s6, E6, N6, phibth6) = np.genfromtxt(filename6,unpack=True)    

    total_a = len(a1) + len(a2) + len(a3) + len(a4) + len(a5) + len(a6)
    a_arr = []
    E_arr = []
    s_arr = []
    phibth_arr = []
    for i in range(0, len(a1),1):
        a_arr.append(a1[i])
        E_arr.append(E1[i])
        s_arr.append(s1[i])
        phibth_arr.append(phibth1[i])
    for i in range(0, len(a2),1):
        a_arr.append(a2[i])
        E_arr.append(E2[i])
        s_arr.append(s2[i])
        phibth_arr.append(phibth2[i])
    for i in range(0, len(a3),1):
        a_arr.append(a3[i])
        E_arr.append(E3[i])
        s_arr.append(s3[i])
        phibth_arr.append(phibth3[i])
    for i in range(0, len(a4),1):
        a_arr.append(a4[i])
        E_arr.append(E4[i])
        s_arr.append(s4[i])
        phibth_arr.append(phibth4[i])
    for i in range(0, len(a5),1):
        a_arr.append(a5[i])
        E_arr.append(E5[i])
        s_arr.append(s5[i])
        phibth_arr.append(phibth5[i])
    for i in range(0, len(a6),1):
        a_arr.append(a6[i])
        E_arr.append(E6[i])
        s_arr.append(s6[i])
        phibth_arr.append(phibth6[i])

    l = sorted(zip(a_arr, E_arr,s_arr,phibth_arr), key=lambda x:x[0])

    for i in range(0,len(a_arr),1):
        am, Eb, sb, pth = l[i]
        print(am,Eb,sb,np.sqrt(pth))


sorter_func()