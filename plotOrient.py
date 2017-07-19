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
file = "Results_IU_ADK"
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
    summary=lines[-1].split(',')
    metadata=lines[-1].split(',')

# now create a nice plot. 
    plt.figure()
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    plt.plot(thetaNSrad,thetaNScorr,'bo')
    plt.plot(thetaEWrad,thetaEWcorr,'go')
    plt.show()
## get the information for the lines
#        calcR = [1., 1.]
#        label1="Baz_calc = %.2f" % (ang)
#        label2="Baz2_calc = %.2f" % (ang2)
#        label3="Baz_meas = %.2f" % (StationAziExpec[2])
## actually plot the things
#        plt.plot(thetaNS,calcR,'blue',label=label1)
#        plt.plot(thetaEW,calcR,'cyan',label=label2)
#        plt.plot(expcTheta,calcR,'black',label=label3)
#        #plt.plot(theta,r,'red',label='Particle Motion')
#        #plt.text(32*np.pi/20,2.7,(printstr+'\n'+
#        #         printstr1+'\n'+printstr2))
#        #fileName =(os.getcwd() +'/'+ resDir +'/Azimuth_'+
#        #        station[0] +'_'+ station[1] +'_'+
#        #        str(eventTime) + '.png')
#        fileName='ADKtestPlot'
#        print fileName
#
#        plt.legend(bbox_to_anchor=(0.8, 0.85, 1., 0.102),loc=3,borderaxespad=0.)
#        plt.savefig(fileName,format='png')
#        plt.show()
#        plt.close()
