#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import os
import glob
import sys
import numpy as np
from scipy.optimize import root
debug = True

class Rotation:
    def __init__(self, stref, sttest):
        self.stref=stref
        self.sttest=sttest

    # This function is used to rotate the data
    # This works on the North/South channel
    def rot1(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data1 = cosd*self.sttest[0].data -sind*self.sttest[1].data
        # Notice we are just checking the N/S we might want to add in the E/W
        resi = (abs(sum(data1*self.stref[0].data)/
                np.sqrt(sum(data1**2)*sum(self.stref[0].data**2)) -1.))
        return resi
        
        
   # This function is used to rotate the data
   # This one works on the East/West channel (stream[1])
    def rot2(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data2 = sind*self.sttest[0].data + cosd*self.sttest[1].data
        # Notice we are just checking the N/S we might want to add in the E/W
        resi = (abs(sum(data2*self.stref[1].data)/np.sqrt(sum(data2**2)*sum(self.stref[1].data**2)) -1.))
        return resi     
        
        
        
        
        

########################################################################
# start of the main program
if __name__ == "__main__":
    net = 'IU'
    station = "*"
# Here is our start and end time
    stime = UTCDateTime('2016-001T00:00:00.0')
    etime = UTCDateTime('2016-002T00:00:00.0')
    ctime = stime

# Grab all the stations
    stas = glob.glob('/msd/' + net + '_*'+station+'/' +  str(stime.year) )
    if debug:
        print stas

# reformat the stations
    newstas = []
    for sta in stas:
        newstas.append(sta.split('/')[2])

    if debug:
        print(newstas)
    stas = newstas
    del newstas

# Now we want to go one day at a time and one station at a time and 
# check the relative orientation

    for sta in stas:
        if debug:
            print('On sta: ' + sta)
        while ctime < etime:
        
            day = str(ctime.julday).zfill(3)
            if debug:
                print('On day: ' + day)
            try:
            # Use the following if True to avoid error handling when you are debugging
            #if True:
                # format the string
                string = '/msd/' + sta + '/' + str(ctime.year) + '/' + day + '/*LH*'
                # read in the data
                # Just grab one hour we might want to change this
                st = read(string, starttime=ctime, endtime=ctime+60.*60)
                if debug:
                    print(st)
                st.detrend('demean')
                st.merge()
                st.filter('bandpass',freqmin=1./8., freqmax=1./4.)
                st.taper(0.05)
                # okay time to process the relative orientation
                # We need to grab the different locations
                locs = []
                for tr in st:
                    locs.append(str(tr.stats.location))
                    locs = list(set(locs))
                    if debug:
                        print(locs)
                    # We now have all the location codes for the station
                    if len(locs) >= 2:
                    # We have at least two sensors so compare the azimuth
                    # First one will be the reference
                        refloc = locs.pop(0)
                        stref = st.select(location=refloc)
                        for loc in locs:
                            sttest = st.select(location=loc)
                            if debug:
                                print(stref)
                                print(sttest)
                            # This is really bad programming we need to fix this    
                            # it is really bad.  what the hell is he doing here?
                            # I got Austin to help and rewrote this as a class!
            
                            #rotdata is an object that stores the data and has
                            #the rotation method.
                            rotdata=Rotation(stref,sttest)
                            #root function - finds the roots of the rotation
                            #method. lm is the levenberg-marquardt method
                            result1 = root(rotdata.rot1, 0., method = 'lm')
                            #not sure what this line is doing.
                            theta1 = result1['x'][0]
                            #what is fun?
                            # This is the value of the residual function you are minimizing
                            resi1 = result1['fun']
                            result2 = root(rotdata.rot2, 0., method = 'lm')
                            #not sure what this line is doing. 
                            # This is grabbing the results from our minimization problem
                            theta2 = result2['x'][0]
                            #what is fun?
                            resi2 = result2['fun'] 

                            if debug:
                                print('Here is theta1: ' + str(theta1))
                                print('Here is the residual1: ' + str(resi1))
                                print('Here is theta2: ' + str(theta2))
                                print('Here is the residual2: ' + str(resi2))
                            if 'f' not in globals():
                                f=open('Results_' + sta, 'w')
                                f.write('ReferenceLoc, TestLoc, day, year, theta 1, residual 1, theta 2, residual 2\n')
                            f.write(refloc +', '+ loc +', '+ day +', ' + str(ctime.year) + ', ' + str(theta1) + ', ' + str(resi1) + ', '  + str(theta2) + ', ' + str(resi2) + '\n')
            except:
                print('Problem with ' + sta + ' on day ' + day)

        #out of the try/except loop
            ctime += 24.*60.*60.
        #out of the station loop - move on to next station...
        ctime = stime
        # We need to close the file and delete the variable otherwise it shows as it is in global
        if 'f' in globals():
            f.close()
            del f
