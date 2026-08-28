# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:11:14 2020

@author: bense
"""

import matplotlib.pyplot as plt
import numpy as np
# import cv2
import skimage.io as skio
import os
os.chdir(r'D:\CorrugatedSheet\Buckling\ReplicationPackage\BasicFunctions')
import sys
sys.path.append(r"C:\Users\bense\AppData\Local\Programs\Python\Python37-32\Lib\site-packages")
print (sys.path)
from CreateData import *
from BasicAnalysis import *

#===============================PARAMS PLOT===================================#
mpl.rcParams['axes.labelsize'] =30
mpl.rcParams['axes.titlesize'] = 30
mpl.rcParams['xtick.labelsize'] =25
mpl.rcParams['ytick.labelsize'] = 25
mpl.rcParams['axes.linewidth']=0.75
mpl.rcParams['xtick.major.width']=0.75
mpl.rcParams['ytick.major.width']=0.75
mpl.rcParams['ytick.major.width']=0.75
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8

mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['legend.markerscale'] = 2
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.handlelength'] = 1.5
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = .2
mpl.rcParams['legend.labelspacing'] =  .5
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['figure.figsize'] =[ 5.12,3.94]
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#===============================PARAMS PLOT=========


#colors
#a= 000, b=001, c= 011,d=010,e=110,f=100
a = (217/255,218/255,218/255)
b = (1,1,0)
c = (1,192/255,0)
d = (150/255,255/255,50/255)
e = (196/255,0,0)
f = (0,250/255,200/255)

def indice(f0):
    df=np.diff(f0)
    jumps=np.where(abs(df)>0.0079)[0]
    dj=np.diff(jumps)
    j=np.where(dj>15)[0]
    idx=jumps[j]
    idx=np.append(idx,jumps[-1])
    return idx


def colorplot(f0,x0,t0):#for the accumulator plot
    idx= indice(f0)
    l=len(idx)
    colors=[a,b,c,d,e,f]
    l_ad=l-(len(colors)+2)
    if l_ad==0:
        colors.append(a)
    else:
        for i in range(l_ad):
            if i%2==0:
                colors.append(e)
            else:
                colors.append(f)
        colors.append(e)
        colors.append(f)
        colors.append(a)
    fig1,ax1=plt.subplots(1,1)
    fig2,ax3=plt.subplots(1,1)
    for i in range(len(idx)):
        if i==0:
            plt.figure(1)
            plt.plot(t0[:idx[0]],x0[:idx[0]],c=colors[0])
            plt.figure(2)
            plt.plot(t0[:idx[0]],f0[:idx[0]],c=colors[0])
        elif i==len(idx)-1:
            plt.figure(1)
            plt.plot(t0[idx[i-1]:idx[i]],x0[idx[i-1]:idx[i]],c=colors[i])
            plt.plot(t0[idx[-1]:],x0[idx[-1]:],c=colors[0])
            plt.figure(2)
            plt.plot(t0[idx[i-1]:idx[i]],f0[idx[i-1]:idx[i]],c=colors[i])
            plt.plot(t0[idx[-1]:],f0[idx[-1]:],c=colors[0])
        else:
            plt.figure(1)
            plt.plot(t0[idx[i-1]:idx[i]],x0[idx[i-1]:idx[i]],c=colors[i])
            plt.figure(2)
            plt.plot(t0[idx[i-1]:idx[i]],f0[idx[i-1]:idx[i]],c=colors[i])
    plt.figure(1)
    # ax1.set_aspect(1.0/ax1.get_data_ratio()*1)
    plt.xlabel(r'$t(s))$')
    plt.ylabel(r'$U(mm)$')
    ax2 = ax1.twinx()
    ax2.set_ylabel(r"$|dF/dt| (a.u)$")
    ax2.plot(t0[1:],np.abs(1000*np.diff(f0)))
    ax2.set_ylim(0,30)
    ax2.set_yticks([0,10,20])
    fig1.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()  
    
    plt.figure(2)
    ax3.set_aspect(1.0/ax3.get_data_ratio()*1)
    plt.xlabel(r'$t(s)$')
    plt.ylabel(r'$F(N)$')
    return 

def colorplotWrpm(f0,x0,t0):#for the RPM plot
    idx= indice(f0)
    l=len(idx)
    colors=[a,b,c,d,c,d,c,d,c,d,c,d,c,b,a]
    l_ad=l-(len(colors)+2)
    fig3,ax4=plt.subplots(1,1)
    fig4,ax6=plt.subplots(1,1)
    for i in range(len(idx)):
        if i==0:
            plt.figure(3)
            plt.plot(t0[:idx[0]],x0[:idx[0]],c=colors[0])
            plt.figure(4)
            plt.plot(t0[:idx[0]],f0[:idx[0]],c=colors[0])
        elif i==len(idx)-1:
            plt.figure(3)
            plt.plot(t0[idx[i-1]:idx[i]],x0[idx[i-1]:idx[i]],c=colors[i])
            plt.plot(t0[idx[-1]:],x0[idx[-1]:],c=colors[0])
            plt.figure(4)
            plt.plot(t0[idx[i-1]:idx[i]],f0[idx[i-1]:idx[i]],c=colors[i])
            plt.plot(t0[idx[-1]:],f0[idx[-1]:],c=colors[0])
        else:
            plt.figure(3)
            plt.plot(t0[idx[i-1]:idx[i]],x0[idx[i-1]:idx[i]],c=colors[i])
            plt.figure(4)
            plt.plot(t0[idx[i-1]:idx[i]],f0[idx[i-1]:idx[i]],c=colors[i])
    plt.figure(3)
    ax4.set_aspect(1.0/ax4.get_data_ratio()*1)
    plt.xlabel(r'$t(s))$')
    plt.ylabel(r'$U(mm)$')
    ax5 = ax4.twinx()
    ax5.set_ylabel(r"$|dF/dt| (a.u)$")
    ax5.plot(t0[1:],np.abs(1000*np.diff(f0)))
    ax5.set_ylim(0,30)
    ax5.set_yticks([0,7.5,15])
    fig3.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()  
    plt.figure(4)
    ax6.set_aspect(1.0/ax6.get_data_ratio()*1)
    plt.xlabel(r'$t(s)$')
    plt.ylabel(r'$F(N)$')
    return 

      
path=r'D:\CorrugatedSheet\Buckling\ReplicationPackage\Figure6\RPM_Break'
exp,filesList=prep_file(path)
plt.close('all')


f_acc,x_acc,t_acc = exp[0].create_var()
f_RPM,x_RPM,t_RPM = exp[1].create_var()



#fig6c
plt.close('all')
colorplot(f_acc,x_acc,t_acc)


#fig6a
colorplotWrpm(f_RPM,x_RPM,t_RPM)

