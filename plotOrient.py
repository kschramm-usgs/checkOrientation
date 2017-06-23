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
for file in glob.glob("Result_*"):
    f = open(file,'r') 
    for line in f:
        data = 






# now create a nice plot. 
            plt.figure()
            ax = plt.subplot(111, projection='polar')
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
# get the particle motion
            #theta = np.arctan2(SignalBHN.data,SignalBHE.data)
            #r = np.sqrt(SignalBHE.data*SignalBHE.data 
            #      + SignalBHN.data*SignalBHN.data)
# get the gradient on the particle motion
            #gradPM = 
# get the information for the lines
            calcR = [1., 1.]
            calcTheta = [np.radians(ang),np.radians(ang+180)]
            calcTheta2 = [np.radians(ang2),np.radians(ang2+180)]
            expcTheta = ([np.radians(StationAziExpec[2]), 
                          np.radians(StationAziExpec[2]+180)])
            label1="Baz_calc = %.2f" % (ang)
            label2="Baz2_calc = %.2f" % (ang2)
            label3="Baz_meas = %.2f" % (StationAziExpec[2])
# actually plot the things
            plt.plot(calcTheta,calcR,'blue',label=label1)
            plt.plot(calcTheta2,calcR,'cyan',label=label2)
            plt.plot(expcTheta,calcR,'black',label=label3)
            #plt.plot(theta,r,'red',label='Particle Motion')
            plt.plot(SignalBHE.data,SignalBHN.data, 'red',label='Particle Motion')
            plt.text(7*np.pi/4,2.5,str(station[1]+' '+ rotated + ' '+str(eventTime)),fontsize=14)
            printstr="linearity %.2f" % (line)
            printstr1="SNR, BHN %.2f" % (SNR_BHN)
            printstr2="SNR, BHE %.2f" % (SNR_BHE)
            plt.text(32*np.pi/20,2.7,(printstr+'\n'+
                     printstr1+'\n'+printstr2))
            fileName =(os.getcwd() +'/'+ resDir +'/Azimuth_'+
                    station[0] +'_'+ station[1] +'_'+
                    str(eventTime) + '.png')
            print fileName
            #plt.close()
            #plt.figure()
            #x = SignalBHE.data**2.
            #y = SignalBHN.data**2.
            #plt.plot(x,y,marker="o")
            #p = np.polyfit(x,y,1)
            #print p
            #ycalc = p[0]*x + p[1]
            #plt.plot(x,ycalc,marker="v")
            #pfitR=np.sqrt(x+ycalc**2.)
            #calcR=[max(pfitR), abs(min(pfitR))]
            #pfitTheta= np.arctan2(ycalc,SignalBHE.data)
            #aveTheta=np.average(pfitTheta)
            #plotAveTheta=[aveTheta, aveTheta+np.pi]
            #print np.degrees(aveTheta)

            #print "Max r "+str(max(r))
            #print "Min r "+str(min(r))
            #print "Max index r "+str(np.argmax(r))
            #print "Min index r "+str(np.argmin(r))
            ##print np.degrees(theta)
            #print "theta at r max "+str(np.degrees(theta[np.argmax(r)]))
            #print "theta at r min "+str(np.degrees(theta[np.argmin(r)]))
            #print "Max index theta "+str(np.argmax(theta))
            #print "Min index theta "+str(np.argmin(theta))
            #print "Max theta "+str(np.degrees(np.max(theta)))
            #print "Min theta "+str(np.degrees(np.min(theta)))

            #label4=("Baz_part = %.2f" % (np.degrees(aveTheta)))
            #plt.plot(plotAveTheta, calcR,'orange',label=label4)
            #plt.plot(pfitTheta, pfitR,'orange',label=label4)

            plt.legend(bbox_to_anchor=(0.8, 0.85, 1., 0.102),loc=3,borderaxespad=0.)
            plt.savefig(fileName,format='png')
            #plt.show()
            plt.close()
