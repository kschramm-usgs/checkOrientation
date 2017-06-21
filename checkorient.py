#!/usr/bin/env python

from obspy.core import read, UTCDateTime
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
    def rot(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data1 = cosd*self.sttest[0].data -sind*self.sttest[1].data
        data2 = sind*self.sttest[1].data + cosd*self.sttest[1].data
        # Notice we are just checking the N/S we might want to add in the E/W
        resi = (abs(sum(data1*self.stref[0].data)/
                np.sqrt(sum(data1**2)*sum(self.stref[0].data**2)) -1.))
        return resi

########################################################################
# start of the main program
if __name__ == "__main__":
    net = 'IU'
    station = "ANTO"
# Here is our start and end time
    stime = UTCDateTime('2016-001T00:00:00.0')
    etime = UTCDateTime('2016-005T00:00:00.0')
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
            #try:
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
            #rotdata is an object that stores the data and has
            #the rotation method.
                        rotdata=Rotation(stref,sttest)
            #root function - finds the roots of the rotation
            #method. lm is the levenberg-marquardt method
                        print('here')
                        result = root(rotdata.rot, 0., method = 'lm')
            #not sure what this line is doing.
                        theta = result['x'][0]
            #what is fun?
                        resi = result['fun'] 
                        print('result: ',str(result))
                        if debug:
                            print('Here is theta: ' + str(theta))
                            print('Here is the residual: ' + str(resi))
                        if 'f' not in globals():
                            print('opening file')
                            f=open('Results_' + sta, 'w')
                            f.write('ReferenceLoc, TestLoc, day, year, theta, residual\n')
                        f.write(refloc +', '+ loc +', '+ day +', ' + str(ctime.year) + ', ' + str(theta) + ', ' + str(resi) + '\n')
                    
                    
                    #sys.exit()
                    
        #    except:
        #        print('Problem with ' + sta + ' on day ' + day)
            ctime += 24.*60.*60.
        ctime = stime
        if 'f' in globals():
            f.close()

