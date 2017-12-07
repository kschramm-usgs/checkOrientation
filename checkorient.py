#!/usr/bin/env python

from obspy.core import read, UTCDateTime, Stream
import os
import glob
import sys
import math
import numpy as np
from time import gmtime, strftime
from scipy.optimize import root
from obspy.io.xseed import Parser
from obspy.signal.cross_correlation import xcorr
import logging
debug = True
debug = False 

def notSameSamp(st1,st2):
    retVal = True
    if (st1.count() == st2.count()):
        retVal = False
    return retVal
  

def threeChannels(numTraces):
    retVal = False
    if (numTraces == 3):
        retVal = True
    return retVal 


class Rotation:
    def __init__(self, stref, sttest):
        self.stref=stref
        self.sttest=sttest

# This function is used to rotate the data
    def rotNS(self,theta):
        theta = theta % 360.
        thetaRad=np.deg2rad(theta)
        #print(theta, thetaRad)
        cosd=np.cos(np.deg2rad(-theta))
        sind=np.sin(np.deg2rad(-theta))
        data1 = cosd*self.sttest[0].data + sind*self.sttest[1].data
        resi = (abs(sum(data1*self.stref[0].data)/
                np.sqrt(sum(data1**2)*sum(self.stref[0].data**2)) -1.))
        return resi

# This function is used to rotate the data
    def rotEW(self,theta):
        theta = theta % 360.
        cosd=np.cos(-np.deg2rad(theta)- np.pi/2.)
        sind=np.sin(-np.deg2rad(theta)- np.pi/2.)
        # trying the next method.  if this doesn't work then I will have to think
        # harder about what adam is doing...
        data2 = sind*self.sttest[0].data - cosd*self.sttest[1].data
        # tried it this way - still coming out with the wrong angle...
        # data2 = sind*self.sttest[1].data - cosd*self.sttest[0].data
        # this was the original - not sure why both had 1 and the n/s uses 0 and 1
        # data2 = sind*self.sttest[1].data - cosd*self.sttest[1].data
        resi = (abs(sum(data2*self.stref[1].data)/
                np.sqrt(sum(data2**2)*sum(self.stref[1].data**2)) -1.))
        return resi
    
# This function gets the correlation of the data
    @staticmethod
    def corrNS(self):
        windowLen = self.stref[0].data.length()/2
        a,corValue = xcorr(self.stref[0].data,self.sttest[0].data,windowLen)
        return corValue

# This function gets the correlation of the data
    def corrEW(self):
        windowLen = self.stref[1].data.length()/2
        a,corValue = xcorr(self.stref[1].data,self.sttest[1].data,windowLen)
        return corValue
        
# this is stuff Adam Ringler added to get information about the station
def getsncl(tr):
    """ Return the sncl """
    nslc = (tr.id).split('.')
    return nslc[1], nslc[0], nslc[3], nslc[2]        

def rotatehorizontal(stream, angle1, angle2):
    """
    A function to rotate the horizontal components of a seismometer from 
    radial and transverse into E and North components.
    """
    debugRot = False
    if stream[0].stats.channel in set(['LHE', 'LHN', 'BHE', 'BHN']):
        stream.sort(['channel'], reverse=True)
        angle1, angle2 = angle2, angle1
    if debugRot:
        print(stream)
        print 'Angle1: ' + str(angle1) + ' Angle2: ' + str(angle2)
    theta_r1 = math.radians(angle1)
    theta_r2 = math.radians(angle2)
    swapSecond = False
    if (angle2 >= 180. and angle2 <= 360.) or angle2 == 0.:
        swapSecond = True 
    # if the components are swaped swap the matrix
    if theta_r1 > theta_r2 and swapSecond:
        if debugRot:
            print 'Swap the components: ' + str((360. - angle1) - angle2)
        stream.sort(['channel'], reverse=True)
        theta_r1, theta_r2 = theta_r2, theta_r1
        if debugRot:
            print(stream)
    # create new trace objects with same info as previous
    rotatedN = stream[0].copy()
    rotatedE = stream[1].copy()
    # assign rotated data
    rotatedN.data = stream[0].data*math.cos(-theta_r1) +\
        stream[1].data*math.sin(-theta_r1)
    rotatedE.data = -stream[1].data*math.cos(-theta_r2-math.pi/2.) +\
        stream[0].data*math.sin(-theta_r2-math.pi/2.)
    rotatedN.stats.channel = 'LHN'
    rotatedE.stats.channel = 'LHE'
    # return new streams object with rotated traces
    streamsR = Stream(traces=[rotatedN, rotatedE])
    return streamsR        

def getorientation(tr, sp):
    """ 
    A function to get the orientation of a station at a specific time from
    the metadata.
    """
    sta, net, chan, loc = getsncl(tr)
    evetime = tr.stats.starttime
    azimuth=None
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

# def


########################################################################
# start of the main program
if __name__ == "__main__":
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    net = 'IU'
    station  = "*"
    # Here is our start and end time
    stime = UTCDateTime('2016-001T00:00:00.0')
    etime = UTCDateTime('2016-366T00:00:00.0')
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
        print ("processing station: "+sta)
# initialize arrays
        thetaNS = [] 
        thetaEW = [] 
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
                st = read(string, starttime=ctime, endtime=ctime+60.*60.*24.)
                if debug:
                    print(st)
            except:
                print('no data for '+ string)
                #better increment....
                ctime += 24.*60.*60.
                continue
            if debug:
                print(st)
            # look for gaps or masked values:
            if (st.get_gaps()): 
                if debug:
                    print('Data has gaps')
                ctime += 24.*60.*60.
                continue
            st.detrend('demean')
            st.merge()
            st.filter('bandpass',freqmin=1./8., freqmax=1./4., zerophase=True, corners=4)
            st.taper(0.05)
        
        # okay time to process the relative orientation
        # We need to grab the different locations of the sensors
            locs = []
            for tr in st:
                locs.append(str(tr.stats.location))
            locs = list(set(locs))
            locs.sort()
            logging.debug('Channels available: '+ str(locs))
            # We now have all the location codes for the statio
            # do we have 2 or three locs?  how do we handle the few stations with 
            # more than 2 locs?
            if len(locs) >= 2:
                logging.info(str(len(locs)))
                if (debug):
                    print (len(locs))
                    print (locs)
       # We have at least two sensors so compare the azimuth
       # First one will be the reference
       # pop will take the first one out and the list gets smaller
                refloc = locs.pop(0)
                stref = st.select(location=refloc)
                stref.sort(['channel'])
       #make sure we have 3 component data 
                if (not threeChannels(stref.count())):
                    if debug:
                        print('No 3 component data: '+ string)
                    #better increment....
                    ctime += 24.*60.*60.
                    continue
                # rotate the reference to metadata
                if (notSameSamp(stref[0],stref[1])):
                    logging.debug('samples in ref data not the same')
                    ctime += 24.*60.*60.
                    continue

                Ref1 = getorientation(stref[0], sp)
                Ref2 = getorientation(stref[1], sp)
                stref = rotatehorizontal(stref,Ref1,Ref2)
       # now loop over the sensors  
                for loc in locs:
                    fileName='Results_' + sta + '_' + refloc + '_' + loc 
                    if (debug):
                        print(fileName)
                    sttest = st.select(location=loc)
                    sttest.sort(['channel'])
        # make sure we have 3 component data 
                    if (not threeChannels(sttest.count())):
                        if debug:
                            print('No 3 component data: '+ string)
                        #better increment....
                        ctime += 24.*60.*60.
                        continue
                    logging.debug('print out stream info')
                    logging.debug(stref)
                    logging.debug(sttest)
                # now make sure that we have the same number of samples
                    if (notSameSamp(stref[0],sttest[0])):
                        logging.debug('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif (notSameSamp(stref[0],sttest[0])):
                        logging.debug('samples not the same')
                        ctime += 24.*60.*60.
                        continue
                    elif (notSameSamp(stref[0],sttest[0])):
                        logging.debug('samples not the same')
                        ctime += 24.*60.*60.
                        continue

        # get metadata orientation values
                    Test1 = getorientation(sttest[0],sp)
                    logging.debug('loc1 orientation '+ str(Test1))
                    Test2 = getorientation(sttest[1], sp)
                    logging.debug('loc2 orientation '+ str(Test2))
                    if (Test1 == None or Test2 == None):
                        print("Cannot find azimuth for test data.") 
                        ctime += 24.*60.*60.
                        continue
                        
        #rotdata is an object that stores the data and has
        #the rotation method.
                    rotdata=Rotation(stref,sttest)
                    sttest = rotatehorizontal(sttest, Test1, Test2)
                    
        #root function - finds the roots of the rotation
        #method. lm is the levenberg-marquardt method
                    resultNS = root(rotdata.rotNS, 0., method = 'lm')
                    resultEW = root(rotdata.rotEW, 0., method = 'lm')
        #grab the results from the minimization problem
        #and check to make sure we are between 0 and 360
                    thetaNS.append(resultNS['x'][0])
                    if abs(resultNS['x'][0]) > 360:
                        thetaNS[-1] = (thetaNS[-1] % 360)
                    elif resultNS['x'][0] < 0:
                        thetaNS[-1] = 360 + thetaNS[-1]

                    thetaEW.append(resultEW['x'][0])
                    if abs(resultEW['x'][0]) > 360:
                        thetaEW[-1] = (thetaEW[-1] % 360)
                    elif resultEW['x'][0] < 0:
                        thetaEW[-1] = 360 + thetaEW[-1]


        #This is the value of the residual function you are minimizing
                    resiNS = resultNS['fun'] 
                    resiEW = resultEW['fun'] 
                    from scipy.stats import pearsonr
                    corrvalNS = pearsonr(stref[0].data,sttest[0].data)[0]
                    corrvalEW = pearsonr(stref[1].data,sttest[1].data)[0]
                    if (debug):
                        print(corrvalNS)
                    
        # if there is a low correlation, don't use that data.  (one station might be noisy, or we ran a calibration)
                    if (abs(corrvalNS) < 0.5) or (abs(corrvalEW) < 0.5):
                        ctime += 24.*60.*60.
                        continue
                        ctime += 24.*60.*60.
                        continue
                       
                    if debug:
                        print('resultNS: ',str(resultNS))
                        print('resultEW: ',str(resultEW))
                        print('Here is theta NS: ' + str(thetaNS))
                        print('Here is the residual EW: ' + str(resiNS))
                        print('Here is theta EW: ' + str(thetaEW))
                        print('Here is the residual EW: ' + str(resiEW))
                    # The previous logic is based on the existence not if python has a copy
                    if debug:
                        print('opening file '+fileName)

                    if not os.path.isfile(fileName):
                        print('writing header info')
                        f=open(fileName, 'w')
                        f.write('Metadata information for ' + sta + \
                                ': Ref theta, ' +  str(Ref1) + ', ' + str(Ref2) + \
                                ', Test theta, ' + str(Test1) +  ', ' + str(Test2) +  '\n')
                        f.write('ReferenceLoc, TestLoc, day, year, comp,'\
                                +'NS theta, NS residual, NS corr, EW theta,'\
                                +'EW residual, EW coor \n')
                        f.close()

                    # Write some results and include metadata
                    # the index [-1] will print the last value in the list
                    f=open(fileName, 'a')
                    f.write(refloc +', '+ loc +', '+ day +', ' +  \
                            str(ctime.year) + ', ' + str(thetaNS[-1]) + ', ' + \
                            str(resiNS) + ', ' + str(corrvalNS) + ', ' + \
                            str(thetaEW[-1]) + ', ' + str(resiEW) + \
                            ', ' + str(corrvalEW) + ', ' + str(Ref1) + ', ' + \
                            str(Test1) + '\n')
                    f.close()
                
        # in the while ctime .lt. etime - need to increment this by a day.
            ctime += 24.*60.*60.
    # done with that station, exit the while loop, reset ctime and numdays
        ctime = stime

        # close the file when you exit the while loop
        if 'f' in globals():
            f.close()

