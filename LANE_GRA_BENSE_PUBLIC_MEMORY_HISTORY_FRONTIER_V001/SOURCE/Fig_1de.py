# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:11:14 2020

@author: bense
"""
#calculates the different between two subsequent images (normalised) and plots
#the sum of the pixels squared

import matplotlib.pyplot as plt
import numpy as np
import cv2
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

mpl.rcParams['lines.markersize'] = 10
mpl.rcParams['legend.markerscale'] = 10
mpl.rcParams['legend.fontsize'] = 2
mpl.rcParams['legend.handlelength'] = 1.5
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = .2
mpl.rcParams['legend.labelspacing'] =  .5
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['figure.figsize'] =[ 5.12,3.94]
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#===============================PARAMS PLOT=========


def pictures(path):
    filesList = []
    os.chdir(eval('path'))
    for root, dirs, files in os.walk(path):
        for filename in files:
            if filename.endswith('tiff'):
                filesList.append(filename)
    return filesList
def last_chars(x,nin=12,nout=-6):
    return(float(x[nin:nout]))
    
def classe(filesList,nin=12,nout=-6):
    filesList_sort=sorted(filesList, key = last_chars) 
    return(filesList_sort)
    
def dif(pics,fps=0.1):
    D=[]
    Dn=[]
    for i in range(len(pics)-1):
        im1=skio.imread(pics[i]).astype(np.double)
        im2=skio.imread(pics[i+1]).astype(np.double)
        delta=im1-im2
        D.append(np.sum(delta**2))
        im1N=(im1-np.mean(im1))/(np.sum((im1-np.mean(im1))**2))**0.5
        im2N=(im2-np.mean(im2))/(np.sum((im2-np.mean(im2))**2))**0.5
        deltaN=im1N-im2N
        Dn.append(np.sum(deltaN**2))
    Dn=np.asarray(Dn)
    D=np.asarray(D)
    t=fps*np.arange(0,len(pics)-1,1)
    return D,Dn,t

def CreateNPY(pics,chemin):
    DELTA=[]   
    os.chdir(chemin)
    for i in range(len(pics)-1):
        im1=skio.imread(pics[i]).astype(np.double)
        im2=skio.imread(pics[i+1]).astype(np.double)
        im1N=(im1-np.mean(im1))/(np.sum((im1-np.mean(im1))**2))**0.5
        im2N=(im2-np.mean(im2))/(np.sum((im2-np.mean(im2))**2))**0.5
        delta=im1N-im2N
        DELTA.append(delta)
    np.save('Diff_imagesnormalise.npy',DELTA)    
    return DELTA

def SaveImage(DELTA,chemin):
    vmin=np.min(DELTA)
    vmax=np.max(DELTA)
    os.chdir(chemin)
    for i in range(shape(DELTA)[0]):
        plt.imshow(DELTA[i],cmap='gray',aspect='auto',vmin=vmin, vmax=vmax)
        plt.savefig(str(i)+'.png')
    return


path=r'D:\CorrugatedSheet\Buckling\ReplicationPackage\Figure1\Pictures'        
#uncomment if Diff_imagesnormalise.npy does not exist
# pic=pictures(path)
# pics=classe(pic)
# D,Dn,temps=dif(pics)
# DELTAN= CreateNPY(pics,path)
os.chdir(r'D:\CorrugatedSheet\Buckling\ReplicationPackage\Figure1\Pictures')
DELTA=np.load('Diff_imagesnormalise.npy')

exp,filesList=prep_file(r'D:\CorrugatedSheet\Buckling\ReplicationPackage\Figure1\Data\1de')
t=exp[0].t
x=exp[0].x
f=exp[0].f
n=exp[0].n

en=int(len(DELTA)/2)

temps=0.1*np.arange(0,len(DELTA),1)

Dn_crop=[]
for i in range (len(DELTA)):
    Dn_crop.append(1000*np.sum(DELTA[i][171:485]**2))



##omega 3
xb1=217
xe1=253
d1=[]    
for i in range (0,len(DELTA),1):
    d1.append(1000*np.sum(DELTA[i][xb1:xe1]**2))
#    
##omega 2
xb2=320
xe2=340
d2=[]
for i in range (0,len(DELTA),1):
    d2.append(1000*np.sum(DELTA[i][xb2:xe2]**2))

##omega 1
xb3=424
xe3=463
d3=[]
for i in range (0,len(DELTA),1):
    d3.append(1000*np.sum(DELTA[i][xb3:xe3]**2))
    
######################
long=n[3]+int((len(x)-n[3])/2)
long_pictures=243+int((len(Dn_crop)-243)/2)

fig,ax=plt.subplots()
plt.plot(t[n[3]:long]-t[n[3]],x[n[3]:long],color=(236/255,28/255,36/255))
plt.plot(t[long:]-t[n[3]],x[long:],color=(66/255,101/255,117/255))#
plt.xlabel(r'$t\,(s)$')
plt.ylabel(r'$U\, (mm)$')
plt.title('fig. 1d')

fig,ax=plt.subplots()
plt.plot(t[n[3]:long]-t[n[3]],f[n[3]:long],color=(236/255,28/255,36/255))
plt.plot(t[long:]-t[n[3]],f[long:],color=(66/255,101/255,117/255))
plt.xlabel(r'$t, (s)$')
plt.ylabel(r'$F\, (N)$')
ratio = 1
xleft, xright = ax.get_xlim()
ax.set_aspect(1.0/ax.get_data_ratio()*ratio)


fig=plt.figure()
plt.plot(temps[239:long_pictures-5]-temps[239],d1[239:long_pictures-5],c=(236/255,28/255,36/255))
plt.plot(temps[long_pictures-5:-8]-temps[239],d1[long_pictures-5:-8],color=(66/255,101/255,117/255))
plt.ylabel(r'$\Delta_3^2$')
plt.xlim(xleft, xright)

fig=plt.figure()
plt.plot(temps[239:long_pictures-5]-temps[239],d2[239:long_pictures-5],c=(236/255,28/255,36/255))
plt.plot(temps[long_pictures-5:-8]-temps[239],d2[long_pictures-5:-8],c=(66/255,101/255,117/255))
plt.ylabel(r'$\Delta_2^2$')
plt.xlim(xleft, xright)

fig=plt.figure()
plt.plot(temps[239:long_pictures-5]-temps[239],d3[239:long_pictures-5],c=(236/255,28/255,36/255))
plt.plot(temps[long_pictures-5:-8]-temps[239],d3[long_pictures-5:-8],c=(66/255,101/255,117/255))
plt.ylabel(r'$\Delta_1^2$')
plt.xlim(xleft, xright)


