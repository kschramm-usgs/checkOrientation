#!/bin/env python
''' read and plot the results from the check orientation script '''
import matplotlib.pyplot as plt
import numpy as np
import sys
import os 
import glob
from obspy import UTCDateTime

#first read in the data - look for files starting with Results*
#path = os.getcwd()
#network="IU"
#station="ANMO"
#refChan="00"
#testChan="10"
#file = "Results_" + network +"_"+ station +"_"+ refChan +"_"+ testChan
# get list of station files:
files=glob.glob('Results_*')
for file in files:
    statInfo = file.split('_')
    network=statInfo[1]
    station=statInfo[2]
    refChan=statInfo[3]
    testChan=statInfo[4]
    
    thetaNS=[]
    thetaNSrad=[]
    thetaEW=[]
    thetaEWrad=[]
    thetaNSresid=[]
    thetaEWresid=[]
    thetaNScorr=[]
    thetaEWcorr=[]
    epochAngles=[]
    epochAngles.append(metadata[4])
    with open(file,'r') as f:
        #data=f.readline()
        mydat=f.read()
        lines=mydat.split('\n')
        metadata=lines[0].split(',')
        print(metadata)
        header=lines[1].split(',')
        #print(header)
        start=lines[2].split(',')
        startDate=str(start[3])+'-'+str(start[2])
        end=lines[-2].split(',')
        endDate=str(end[3])+'-'+str(end[2])
        metaAngOld=metadata[4]
        epoch=1
        for ln in lines[2:-1]:
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
            metaAng = lv[10])
            # check for changes in epoch
            if metaAng /= metaAngOld:
                print('You have another epoch')
                epoch=epoch+1
                epochAngles.append(metaAng)
                 


        # moving this to the plotting routine...
        numdays = len(thetaNS)
        NSorient=np.average(thetaNS)
        EWorient=np.average(thetaEW)
        NSstr=str("%.2f" % NSorient)
        EWstr=str("%.2f" % EWorient)
        thetaNSstd=np.std(thetaNS)
        thetaEWstd=np.std(thetaEW)
        NSstdstr=str("%.2f" % thetaNSstd)
        EWstdstr=str("%.2f" % thetaEWstd)

    # now create a nice plot. 
        plt.figure(figsize=(11,8.5))
        ax = plt.subplot(111, projection='polar')
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        #ax.set_ylabel=('Correlation')
        plt.set_ylabel=('Correlation')
        ax.set_title(station+' ,Ref Chan: '+ refChan + ' ,Test Chan: ' \
                + testChan + ', for ' + startDate + ' to ' + endDate )
        print(NSorient, EWorient)
        label1="NS"
        label2="EW"
        #print(isinstance(metadata[4],str))
        #metadataRad=np.deg2rad(metadata[4])
        plt.ylim([0, 1.2])

        for ang in epochAngles:
            # need to change the metadata ot be ang.
            plt.arrow(np.deg2rad(float(metadata[4])),0, 0,1.09,fc='b', ec='b',head_width=0.05, head_length = 0.1, alpha=0.8)
            plt.arrow(np.deg2rad(float(metadata[5])),0, 0,1.09,fc='g', ec='g',head_width=0.05, head_length = 0.1, alpha=0.8)

        plt.plot(thetaNSrad,thetaNScorr,'bo', label=label1, alpha =0.35)
        plt.plot(thetaNSrad,thetaNSresid,'bx', label=label1, alpha =0.35)
        plt.plot(thetaEWrad,thetaEWcorr,'go', label=label2, alpha=0.35)
        print(np.deg2rad(float(metadata[5])))
        plt.legend(bbox_to_anchor=(0.95, 0.85, 1.2, 0.102),loc=3,borderaxespad=0.)
        plotString = str("NS metadata, calculated (std): " + str(metadata[4]) +\
                ", " + NSstr +"(" + NSstdstr + ")\nEW metadata, calculated (std): " + str(metadata[5]) +\
                ", " + EWstr +"(" + EWstdstr + ")\n" + str(numdays) + " days in calculation")
        plt.text(19*np.pi/20,1.4,plotString,fontsize=12)
        plt.show()
