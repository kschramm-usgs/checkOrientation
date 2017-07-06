#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import os
import glob
import sys
import math
import numpy as np
from time import gmtime, strftime
from scipy.optimize import root
from obspy.io.xseed import Parser
debug = True

class Rotation:
    def __init__(self, stref, sttest):
        self.stref=stref
        self.sttest=sttest

# This function is used to rotate the data
    def rotNS(self,theta):
        theta = theta % 360.
        cosd=np.cos(np.deg2rad(-theta))
        sind=np.sin(np.deg2rad(-theta))
        data1 = cosd*self.sttest[0].data + sind*self.sttest[1].data
        resi = (abs(sum(data1*self.stref[0].data)/
                np.sqrt(sum(data1**2)*sum(self.stref[0].data**2)) -1.))
        return resi

# This function is used to rotate the data
    def rotEW(self,theta):
        theta = theta % 360.
        cosd=np.cos(-np.deg2rad(theta)- math.pi/2.)
        sind=np.sin(-np.deg2rad(theta)- math.pi/2.)
        data2 = sind*self.sttest[1].data - cosd*self.sttest[1].data
        resi = (abs(sum(data2*self.stref[1].data)/
                np.sqrt(sum(data2**2)*sum(self.stref[1].data**2)) -1.))
        return resi
    
# This function gets the correlation of the data
    def corrNS(self, windowLen):
        windowLen = self.stref[0].data.length()/2
        corValue = (self.stref,self.sttest,windowLen)
        return corValue
        
def getsncl(tr):
    """ Return the sncl """
    nslc = (tr.id).split('.')
    return nslc[1], nslc[0], nslc[3], nslc[2]        
        

def getorientation(tr, sp):
    """ 
    A function to get the orientation of a station at a specific time from
    the metadata.
    """
    sta, net, chan, loc = getsncl(tr)
    evetime = tr.stats.starttime
    for cursta in sp.stations:
# As we scan through blockettes we need to find blockettes 50 and 52
        for blkt in cursta:
            if blkt.id == 50:
# Pull the station info for blockette 50
                stacall = blkt.station_call_letters.strip()
            if stacall == sta:
                if blkt.id == 52 and blkt.location_identifier == loc and blkt.channel_identifier == chan:
                    if type(blkt.end_date) is str:
                        curdoy = strftime("%j", gmtime())
                        curyear = strftime("%Y", gmtime())
                        curtime = UTCDateTime(curyear + "-" +
                                              curdoy + "T00:00:00.0")
                        if blkt.start_date <= evetime:
                            azimuth = blkt.azimuth
                    elif blkt.start_date <= evetime and blkt.end_date >= evetime:
                        azimuth = blkt.azimuth
    return azimuth


########################################################################
# start of the main program
if __name__ == "__main__":
    net = 'IU'
    station = "*"
# Here is our start and end time
    stime = UTCDateTime('2017-001T00:00:00.0')
    etime = UTCDateTime('2017-031T00:00:00.0')
    ctime = stime

    sp = Parser('/APPS/metadata/SEED/' + net + '.dataless')


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
        while ctime <= etime:
        
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
            # look for gaps or masked values:
            if (st.get_gaps()): 
                print('Data has gaps')
                ctime += 24.*60.*60.
                continue
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
                stref.sort(['channel'])
        #make sure we have 3 component data 
                if (stref.count() < 3) :
                    print('No 3 component data: '+ string)
                    #better increment....
                    ctime += 24.*60.*60.
                    continue
        # now the test stream
                for loc in locs:
                    sttest = st.select(location=loc)
                    sttest.sort(['channel'])
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
                # This is clunky and probably needs to be changed
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
                        print('Here is theta NS: ' + str(thetaNS))
                        print('Here is the residual EW: ' + str(resiNS))
                        print('Here is theta EW: ' + str(thetaEW))
                        print('Here is the residual EW: ' + str(resiEW))
                    # we seem to have passed all the tests, so see if there is a file that needs opening
                    fileName='Results_' + sta
                    # The previous logic is based on the existence not if python has a copy
                    if debug:
                        print('opening file '+fileName)
                    # this if not in globals doesn't seem to work. going back to the other test.
                    #if 'f' not in globals():
                    # write results to file.
                    if not os.path.isfile(fileName):
                        f=open(fileName, 'w')
                        f.write('ReferenceLoc, TestLoc, day, year, comp,' \
                                'NS theta, NS residual, EW theta, EW residual, metadata Ref LH1,' \
                                 'metadata Ref LH2, metadata Test LH1, metadata Test LH2\n')
                    # get metadata orientation values
                    Ref1[nd] = getorientation(stref[0], sp)
                    Ref2[nd] = getorientation(stref[1], sp)
                    Test1[nd] = getorientation(sttest[0],sp)
                    Test2[nd] = getorientation(sttest[1], sp)
                    # We have some results so lets also include metadata
                    
                    f.write(refloc +', '+ loc +', '+ day +', ' +  \
                            str(ctime.year) + ', ' + str(thetaNS) + ', ' + \
                            str(resiNS) + ', ' + str(thetaEW) + ', ' + str(resiEW) + \
                            ', ' + str(Ref1) + ', ' + str(Ref2) + ', ' + str(Test1) + ', ' + str(Test2) +  '\n')
                
        # in the while ctime .lt. etime - need to increment this by a day.
            ctime += 24.*60.*60.
        # keep track of the number of days we have...           
            nd+=1
    # calculate some statistics...
        Ref1Ave = Ref1Total/nd
        Ref2Ave = Ref2Total/nd
        Test1Ave = Test1Total/nd
        Test2Ave = Test2Total/nd
        Ref1Std = np.std
        f.write(
    # done with that station, exit the while loop, reset ctime and numdays
        nd = 0 
        ctime = stime
    # calculate the standard deviation and average    

        # close the file when you exit the while loop

        if 'f' in globals():
            f.close()

        #    except:
        #        print('Problem with ' + sta + ' on day ' + day)
