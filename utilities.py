
'''
$Id: utilities.py

Some useful tools for NEDscrape.
'''

__author__ = "David M. French - frenchd24@gmail.com"

import sys
import math
import numpy as np

################################################################################

def isNumber(s):
    '''
    Determines if input is a number or not.

    Inputs
    ------
        s  : anything you like

    Returns
    -------
        boolean
    '''

    try:
        float(s)
        if np.isnan(float(s)):
            return False
        else:
            return True

    except Exception as e:
        return False


def isNull(s):
    # determine if something is a null value or not. returns a boolean

    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'

    returnVal = False
    try:
        l = float(s)
        if np.isnan(l):
            returnVal = True

        if l == nullFloat:
            returnVal = True

        if l == nullInt:
            returnVal = True

    except ValueError:
        # ValueError means float(s) failed because s is a string
        # so check if it is the null string
        if s == nullStr:
            returnVal =  True

    return returnVal


def notNull(s):
    # determine if something is a null value or not. returns a boolean

    nullFloat = -99.99
    nullInt = -99
    nullStr = 'x'

    returnVal = True
    try:
        l = float(s)
        if np.isnan(l):
            returnVal = False

        if l == nullFloat:
            returnVal = False

        if l == nullInt:
            returnVal = False

    except ValueError:
        # ValueError means float(s) failed because s is a string
        # so check if it is the null string
        if s == nullStr:
            returnVal =  False

    return returnVal



def bfind(k,s):
    '''
    This is a find function that returns true/false instead of an index.

    Inputs:
    ------
    k   :   str
        Input string being searched
    s   :   str
        String being searched for

    Returns:
    ------
    True/False  : boolean
        Was 's' found in 'k'?
    
    '''

    if str(k).find(s) == -1:
        return False
    else:
        return True



def calculatevcorr(ra,dec,velocity):
    rav = 186.7833
    decv = 12.9333
    vcorr = velocity + 300*(np.sin(np.radians(dec)) * np.sin(np.radians(decv)) + np.cos(np.radians(dec))*np.cos(np.radians(decv)) * np.cos(np.radians(ra-rav)))
    return vcorr



def trunc(x,n):
    # truncates a number, x, to n places without rounding
    if isNumber(x):
        x = float(x)
        slen = len('%.*f' %(n,x))
        return str(x)[:slen]
    else:
        return x



def calculateInclination(major,minor):
    # outputs inclination in degress
    inclination = np.arccos(float(minor)/float(major))*180/np.pi
    return inclination



def calculateFancyInclination(maj,min,q0):
    # calculate more advanced version of inclination

#     q0 = 0.2
#     if float(min) > float(maj):
#         maj,min = min,maj
#
    q = float(min)/float(maj)
    if q >1:
        q = 1/q

    # if the computed axis ratio is smaller than the minimum, q0, avoid the domain error
    # by setting inc = 90 immediately
    if q > q0:
#         print 'q, q0',q,', ',q0
        inc = np.arccos(np.sqrt((q**2 - q0**2)/(1-q0**2)))*(180.0/np.pi)
    else:
        inc = 90.

    return inc




def convertRAandDec(ra,dec,toWhat):
    # this function converts ra and dec in degrees to HH:MM:SS and DD:MM:SS

    if toWhat == 'sexagesimal':
        raHours,raMinutes = str(((float(ra)*24)/360)).split('.')
        raMinutes,raSeconds = str(float('0.'+raMinutes) *60).split('.')
        raSeconds = float('0.'+raSeconds) *60

        decDegrees,decMinutes = str(dec).split('.')
        decMinutes,decSeconds = str(float('0.'+decMinutes)*60).split('.')
        decSeconds = float('0.'+decSeconds)*60

        return (int(raHours),int(raMinutes),round(raSeconds,4)),(int(decDegrees),int(decMinutes),round(decSeconds,4))

    elif toWhat == 'degree':
        raHours,raMinutes,raSeconds = ra
        decDegrees,decMinutes,decSeconds = dec

        raAfterDecimal = ((float(raSeconds)/60)+float(raMinutes))/60
        decAfterDecimal = ((float(decSeconds)/60)+float(decMinutes))/60

        if float(raHours) < 0:
            raDecimal = -15*(abs(float(raHours)) + raAfterDecimal)
        else:
            raDecimal = 15*(float(raHours) + raAfterDecimal)
        # check again
        if str(raHours).find('-') !=-1 and raDecimal >0:
            raDecimal = -raDecimal

        if float(decDegrees) < 0:
            decDecimal = -(abs(float(decDegrees)) + decAfterDecimal)
        else:
            decDecimal = (float(decDegrees) + decAfterDecimal)
        # check again
        if str(decDegrees).find('-') !=-1 and float(decDecimal) >0:
            decDecimal = - decDecimal

    else:
        sys.stdout.write("3rd argument must match either 'sexagesimal' or 'degree'")
        raDecimal = False
        decDecimal = False

    return raDecimal,decDecimal



def calculateAngularSeparation(ra1,dec1,ra2,dec2):
    # return angular separation between ra1,dec1 and ra2, dec2 in radians

#     ra1 = float(ra1)*(math.pi/180)
#     dec1 = (90-float(dec1))*(math.pi/180)
#     ra2 = float(ra2)*(math.pi/180)
#     dec2 = (90- float(dec2))*(math.pi/180)

    ra1 = np.radians(float(ra1))
    dec1 = np.radians(float(dec1))
    ra2 = np.radians(float(ra2))
    dec2 = np.radians(float(dec2))

    inside = np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2)
    if inside <1.:
        angSep = np.arccos(inside)
    else:
        angSep = 0.

    return angSep


def calculateImpactParameter(ra1,dec1,ra2,dec2,dist):
    # return physical separation between ra1,dec1 and ra2, dec2 in kpc

    angSep = calculateAngularSeparation(ra1,dec1,ra2,dec2)
    sep = round(np.tan(angSep) * float(dist) * 1000,3)
    return sep



def calculateAngularSeparation_slow(ra1,dec1,ra2,dec2):
    # returns the angular separation between supplied coordinates in arcsec
    # uses the astropy module, which is slower

    from astropy.coordinates import SkyCoord
    from astropy import units as u

    c1 = SkyCoord(ra1,dec1,frame='icrs',unit='degree')
    c2 = SkyCoord(ra2,dec2,frame='icrs',unit='degree')
    sep = c1.separation(c2)
    return sep.radian



def calculateImpactParameter_slow(ra1,dec1,ra2,dec2,dist):
    # returns the linear distance in kpc between two objects
    # uses the Astropy module, which is slower

    angSep = calculateAngularSeparation_slow(ra1,dec1,ra2,dec2)
    sep = round(np.tan(angSep) * dist * 1000,3)
    return sep




def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return np.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return np.arccos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def calculateAzimuth(galaxyRA,galaxyDec,AGNra,AGNdec,galaxyDist,pa):
    # calculates the QSO-galaxy azimuth

    gRA = float(galaxyRA)
    gDec = float(galaxyDec)
    agnRA = float(AGNra)
    agnDec = float(AGNdec)
    galaxyDist = float(galaxyDist)
    pa = float(pa)

    # calculate angular separations in ra and dec to determine positions on chart w.r.t. target AGN

    try:
        dRA = calculateImpactParameter(gRA,agnDec,agnRA,agnDec,galaxyDist)
        dDec = calculateImpactParameter(agnRA,gDec,agnRA,agnDec,galaxyDist)

        if gRA < agnRA:
            dRA = -dRA

        if gDec < agnDec:
            dDec = -dDec

        agnVector = np.array([dRA,dDec])

        if pa <90.:
            galaxyVector = np.array([np.tan(pa*(np.pi/180)),1.])
        if pa == 90.:
            galaxyVector = np.array([1,0.])
        if pa >90:
            galaxyVector = np.array([np.tan((180-pa)*(np.pi/180)),-1.])

        azimuth = angle(agnVector,galaxyVector)* (180/np.pi)

        if azimuth >90:
            azimuth = 180-azimuth

    except Exception as e:
        sys.stdout.write("Could not calculate azimuth. Here's the message built into the exception: \n %s \n" %e)
        azimuth = 'x'

    if isNumber(azimuth):
        azimuth = round(azimuth,1)

    return azimuth




def calculateVirialRadius(d):
    # this routine calculates the virial radius of a galaxy using the diameter
    # approximation found in Wakker & Savage 2009 and Wakker et al 2015.

    logRvir = 0.69*np.log10(float(d)) + 1.24

    rVir = 10**logRvir

    return rVir



def calculateLinearDiameters(major,minor,distance):
    # input major and minor in arcsec, distance in Mpc
    # outputs major and minor in kpc

    newMajor = np.tan(np.radians(float(major)/3600.))*(distance*1000)
    newMinor = np.tan(np.radians(float(minor)/3600.))*(distance*1000)
    return (newMajor,newMinor)



def calculateLikelihood_basic(impact,d):
    # input is impact parameter and major axis diameter.
    #
    # basic version has no velocity component

    rVir = calculateVirialRadius(d)
    l = np.exp(-(float(impact)/rVir)**2)

    return l


def weightedMean(x,e):
    # x is the array of data
    # e is the array of errors (i.e. std dev for each x_i)
    #
    # returns the weighted mean and the error on the weighted mean

    topSum = 0
    bottomSum = 0

    for xi,ei in zip(x,e):
        top = float(xi)/(float(ei)**2)
        topSum = topSum + top

        bottom = 1./(float(ei)**2)
        bottomSum = bottomSum + bottom

    # the weighted mean
    xm = topSum/bottomSum

    topStd = 1
    bottomStd = 0
    # now the weighted Std
    for ei in e:

        bottom = 1./(float(ei)**2)
        bottomStd = bottomStd + bottom

    xstd = np.sqrt(topStd/bottomStd)

    return xm, xstd


def round_to_sig(x,sig=3):
    # round to sig number of significant figures
    return round(x, sig-int(np.floor(np.log10(abs(x))))-1)


def median_low(l):
    # returns the closest element in a list to the median, rounding down

    # E.g.,
    # list = [1, 2, 3, 4, 5]
    # median_low(list) = 3
    #
    # list = [1, 2, 3, 4]
    # median_low(list) = 2

    # sort the input list
    try:
        l.sort()
    except AttributeError:
        l = list(l)
        l.sort()

    except Exception as e:
        sys.stdout.write('Could not sort the following: {0}'.format(l))
        sys.stdout.write("The error built into the exception was: {0}".format(e))
        sys.exit()


    # make an array version of the input list
    l_a = np.array(l)

    # calculate the median
    med = np.median(l_a)

    # difference between the median and each value in the list
    diff = abs(np.subtract(med,l_a))

    # round both off at 6 decimal places to overcome rounding errors
    mindiff = round(min(diff),6)
    diff2 = np.around(diff, decimals=6)

    # find the index where the first instance of the value closest to the median appears
    indexMin = np.where(diff2==mindiff)
    indexMin = indexMin[0][0]

    return l[indexMin]
