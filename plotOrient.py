#!/bin/env python
''' read and plot the results from the check orientation script '''
import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import glob
from obspy import UTCDateTime

#first read in the data
path = os.getcwd()
file = "testResults_IU_ADK"
station="ADK"
thetaNS=[]
thetaNSrad=[]
thetaEW=[]
thetaEWrad=[]
thetaNSresid=[]
thetaEWresid=[]
thetaNScorr=[]
thetaEWcorr=[]
with open(file,'r') as f:
    data=f.readline()
    mydat=f.read()
    lines=mydat.split('\n')
    header=lines[0].split(',')
    for ln in lines[1:-3]:
        lv=ln.split(',')
# break up data into arrays
        thetaNS.append(float(lv[4]))
        thetaNSrad.append(np.deg2rad(float(lv[4])))
        thetaNSresid.append(lv[5])
        thetaNScorr.append(lv[6])
        thetaEW.append(float(lv[7]))
        thetaEWrad.append(np.deg2rad(float(lv[7])))
        thetaEWresid.append(lv[8])
        thetaEWcorr.append(lv[9])
    summary=lines[-3].split(',')
    print(summary)
    metadata=lines[-2].split(',')
    print(metadata)

# now create a nice plot. 
    plt.figure(figsize=(11,8.5))
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    #ax.set_ylabel=('Correlation')
    plt.set_ylabel=('Correlation')
    ax.set_title(station+'                 ')
    NSorient=float(summary[1][1:-1])
    NSstr=str("%.2f" % NSorient)
    EWorient=float(summary[5][1:-1])
    EWstr=str("%.2f" % EWorient)
    print(NSorient, EWorient)
    label1="NS"
    label2="EW"
    #print(isinstance(metadata[4],str))
    #metadataRad=np.deg2rad(metadata[4])
    plt.ylim([0, 1.2])
    plt.arrow(np.deg2rad(float(metadata[4])),0, 0,1.09,fc='b', ec='b',head_width=0.05, head_length = 0.1, alpha=0.8)
    plt.plot(thetaNSrad,thetaNScorr,'bo', label=label1, alpha =0.35)
    plt.arrow(np.deg2rad(float(metadata[5])),0, 0,1.09,fc='g', ec='g',head_width=0.05, head_length = 0.1, alpha=0.8)
    plt.plot(thetaEWrad,thetaEWcorr,'go', label=label2, alpha=0.35)
    print(np.deg2rad(float(metadata[5])))
    plt.legend(bbox_to_anchor=(0.95, 0.85, 1.2, 0.102),loc=3,borderaxespad=0.)
    plotString = str("NS metadata, calculated: " + str(metadata[4]) +\
            ", " + NSstr + "\nEW metadata, calculated: " + str(metadata[5]) +\
            ", " + EWstr)
    plt.text(5*np.pi/6,1.4,plotString,fontsize=12)
    plt.show()
