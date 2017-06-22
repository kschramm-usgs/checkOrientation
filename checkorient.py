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
    def rotNS(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data1 = cosd*self.sttest[0].data -sind*self.sttest[1].data
        resi = (abs(sum(data1*self.stref[0].data)/
                np.sqrt(sum(data1**2)*sum(self.stref[0].data**2)) -1.))
        return resi

# This function is used to rotate the data
    def rotEW(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(theta))
        sind=np.sin(np.deg2rad(theta))
        data2 = sind*self.sttest[1].data + cosd*self.sttest[1].data
        resi = (abs(sum(data2*self.stref[1].data)/
                np.sqrt(sum(data2**2)*sum(self.stref[1].data**2)) -1.))
        return resi
########################################################################
# start of the main program
if __name__ == "__main__":
    net = 'IU'
    station = "GUMO"
# Here is our start and end time
    stime = UTCDateTime('2016-001T00:00:00.0')
    etime = UTCDateTime('2016-006T00:00:00.0')
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
        # format the string
            string = '/msd/' + sta + '/' + str(ctime.year) + '/' + day + '/*LH*'
        # read in the data
        # Just grab one hour we might want to change this
            try:
                st = read(string, starttime=ctime, endtime=ctime+60.*60)
            except:
                print('no data for '+ string)
                #better increment....
                ctime += 24.*60.*60.
                continue
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
        #make sure we have 3 component data 
                if (stref.count() < 3) :
                    print('No 3 component data: '+ string)
                    #better increment....
                    ctime += 24.*60.*60.
                    continue
        # now the test stream
                for loc in locs:
                    sttest = st.select(location=loc)
        # make sure we have 3 component data 
                    if (sttest.count() < 3) :
                        print('No 3 component data: '+ string)
                        #better increment....
                        ctime += 24.*60.*60.
                        continue
                    if debug:
                        print(stref)
                        print(sttest)
                # now make sure that we have the same number of samples
                    if (stref[0].count() != stref[1].count()):
                        print('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif (stref[0].count() != stref[2].count()):
                        print('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif stref[0].count() != sttest[0].count():
                        print('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif stref[0].count() != sttest[1].count():
                        print('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif stref[0].count() != sttest[2].count():
                        print('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    # we seem to have passed all the tests, so see if there is a file that needs opening
                    fileName='./Results_'+sta
                    if not os.path.isfile(fileName):
                        print('opening file '+fileName)
                        f=open(fileName, 'w')
                        f.write('ReferenceLoc, TestLoc, day, year, comp, \
                                NS theta, NS residual, EW theta, EW residual\n')
        #rotdata is an object that stores the data and has
        #the rotation method.
                    rotdata=Rotation(stref,sttest)
        #root function - finds the roots of the rotation
        #method. lm is the levenberg-marquardt method
                    resultNS = root(rotdata.rotNS, 0., method = 'lm')
                    resultEW = root(rotdata.rotEW, 0., method = 'lm')
        #grab the results from the minimization problem
                    thetaNS = resultNS['x'][0]
                    thetaEW = resultEW['x'][0]
        #This is the value of the residual function you are minimizing
                    resiNS = resultNS['fun'] 
                    resiEW = resultEW['fun'] 
                    print('resultNS: ',str(resultNS))
                    print('resultEW: ',str(resultEW))
                    if debug:
                        print('Here is thetaNS: ' + str(thetaNS))
                        print('Here is the residualEW: ' + str(resiEW))
                        print('Here is thetaEW: ' + str(thetaEW))
                        print('Here is the residualEW: ' + str(resiEW))
                   
                    # write results to file.
                    f.write(refloc +', '+ loc +', '+ day +', ' +  \
                            str(ctime.year) + ', ' + str(thetaNS) + ', ' + \
                            str(resiNS) + ', ' + str(thetaEW) + ', ' + str(resiEW) + '\n')
                
            # in the while ctime .lt. etime - need to increment this by a day.
            ctime += 24.*60.*60.
        # done with that station, exit the while loop.
        ctime = stime

        # close the file when you exit the while loop

        if 'f' in globals():
            f.close()

        #    except:
        #        print('Problem with ' + sta + ' on day ' + day)
