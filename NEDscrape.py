'''
NEDscrape Navigator

A program to grab galaxy data from the NED server by means of xml-based VOTable
distributions and html parsing, and recompile specified data into a new csv file(s)

'''
__author__ = "David M. French - frenchd24@gmail.com"
__version__="2.1"
__date__ = "10/09/2020"

import getpass
import csv
import math
import optparse
import os
from queue import Queue
import string
import sys
import time
import threading
# import urllib2
from urllib.request import urlopen, urlretrieve
from urllib.parse import quote, quote_plus, unquote_plus
# from urllib.parse import unquote_plus
import warnings

from astropy.io.votable import parse, tree
import numpy as np
import pandas as pd

# these are the parameters for this code, should be in the same directory
# from . import params
import params

# some useful tools for this program
from utilities import *

def parse_commandline():
    # Parse the options given on the command-line.

    parser = optparse.OptionParser(usage=__doc__,\
                                   version="$Id: NEDscrape, v {} {}".format(__version__, __date__))
    parser.add_option("-f", "--filename", dest='filename', action='store', type='string',\
            help="The name of the file containing the galaxy names to be queried")
    parser.add_option("-o", "--output", dest='output',action='store',\
        help="Specify the directory in which to save newly generated csv table")
    parser.add_option("-n", "--outname", dest='outname',action='store',\
        help="Specify the output filename (etc.csv) for the newly generated csv table")
    parser.add_option("-v","--verbose",action="store_true",\
        default=False, help="Run verbosely")

    opts,args = parser.parse_args()

    # check if necessary input exists

    if opts.filename is None and opts.outname is None:
        print()
        print("NEDscrape.py  version: {}   {}".format(__version__, __date__))
        print()
        print("A program to grab galaxy data from the NED server.")
        print()
        print("Normal usage example: python NEDscrape.py -f names.txt -n outFile -o /Users/me/ -v")
        print()
        print(" -f is a text file of names of objects, one per line, for which data will be retrieved")
        print(" -n is the name of the csv file that will be created and the NED data written to")
        print(" -o is the full pathname describing where the file should be saved")
        print(" -v tells it to run verbosely")
        print()
        sys.exit(1)

    test = True
    if opts.filename is None:
        print("--filename is a required parameter")
        test = False

    if opts.outname is None:
        print("--outname is a required parameter")
        opts.outname = opts.filename + '_results'

    if opts.output is None:
        print("--output is a required parameter")
        # set output to current directory
        opts.output = os.getcwd() + '/'

    if test == False:
        sys.exit(1)

    # make an output directory if it doesn't exist yet
#     if not os.path.exists(opts.output): os.makedirs(opts.output)

    # show parameters
    if opts.verbose:
        print()
        print("********************** PARAMETERS ******************************")
        print('filename: ', opts.filename)
        print('output directory: ', opts.output)
        print('output filename: ',opts.outname)
        print()
        print("************************ RESULTS *******************************")

    return opts


class ThreadUrl(threading.Thread):
    """Threaded Url Grab"""
    def __init__(self, queue, outQueue, numberQueue):
        threading.Thread.__init__(self)
        self.queue = queue
        self.outQueue = outQueue
        self.numberQueue = numberQueue

    def run(self):
        while True:
            #grabs host from queue
            host = self.queue.get()

            #grabs urls of hosts and reads or parses them according to their type
            try:
                # url = urllib2.urlopen(host[0])
                # url = urlopen(host[0])
                url_obj, headers = urlretrieve(host[0])

            except Exception as e:
                sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s\n" %e)

            if host[1] ==1:
                warnings.simplefilter("ignore")
                votable1 = 'x'
                try:
                    # votable1 = parse(url, pedantic=False)
                    votable1 = parse(url_obj, pedantic=False)
                except Exception as e:
                    sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)

                self.outQueue.put(votable1)
                warnings.resetwarnings()

            else:
                html = 'x'
                try:
                    # html = url.readlines()
                    with open(url_obj) as f:
                        html = f.readlines()
                except Exception as e:
                    sys.stderr.write("\n Unable to read html. Here is the message built into the exception: \n %s \n" %e)

                self.outQueue.put(html)

            self.numberQueue.put(host[1])
            self.queue.task_done()
            # url.close()


class galaxyClass(object):
    # a galaxy object containing it's information and methods of retrieving and
    # returning it

#     def __init__(self,name,votable,html,fullHtml):
#         self.name = name
#         self.votable = votable
#         self.html = html
#         self.fullHtml = fullHtml

    def __init__(self,name):
        self.name = name


    def queryNED(self):
        '''
        query NED server for each object name, returning the file object representing
        the xml VOTable output from NED as well as the redshift-independent distance html file

        Note: htmlFileObject is the html page for redshift independent distance, while
        fullHtmlFileObject is the whole html page for the object

        '''
#         queue = Queue.Queue()
#         outQueue = Queue.Queue()
#         voList = []


        # queue = Queue.Queue()
        # outQueue = Queue.Queue()
        # numberQueue = Queue.Queue()
        queue = Queue()
        outQueue = Queue()
        numberQueue = Queue()

        t1 = ThreadUrl(queue,outQueue,numberQueue)
        t2 = ThreadUrl(queue,outQueue,numberQueue)
        t3 = ThreadUrl(queue,outQueue,numberQueue)
        t1.setDaemon(True)
        t2.setDaemon(True)
        t3.setDaemon(True)
        t1.start()
        t2.start()
        t3.start()

        name = self.name

        # maybe use this instead? (03/11/18)
#        https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=xml_all&zv_breaker=30000.0&list_limit=5&img_stamp=YES

        hosts = [["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname=%s"\
        %name,1],["http://ned.ipac.caltech.edu/cgi-bin/nDistance?name=%s"\
        %name,2],["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES"\
        %name,3]]

        # populate queue with data
        for host in hosts:
#             time.sleep(0.1)
            queue.put(host)

        # wait on the queue until everything has been processed
        queue.join()

        # decide which object is which in the queue
        first = outQueue.get()
        second = outQueue.get()
        third = outQueue.get()

        number1 = numberQueue.get()
        number2 = numberQueue.get()
        number3 = numberQueue.get()

        if number1 == 1:
            votable = first
        elif number1 == 2:
            html = first
        elif number1 == 3:
            fullHtml = first

        if number2 == 1:
            votable = second
        elif number2 == 2:
            html = second
        elif number2 == 3:
            fullHtml = second

        if number3 == 1:
            votable = third
        elif number3 == 2:
            html = third
        elif number3 == 3:
            fullHtml = third

        self.votable = votable
        self.html = html
        self.fullHtml = fullHtml


    def returnRedshift(self):
        # find and return redshift for each object

        redshift = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            redshift = table.array['Redshift'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return redshift. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if isNumber(redshift):
            if math.isnan(float(redshift)):
                redshift = 'x'
        return redshift


    def returnJ2000Position(self):
        # find and return equatorial J2000 position RA and Dec for each object

        ra = 'x'
        dec = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
#             ra = table.array['RA(deg)'][0]
#             dec = table.array['DEC(deg)'][0]
            ra = table.array['RA'][0]
            dec = table.array['DEC'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return J2000 position. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if isNumber(ra) or isNumber(dec):
            if math.isnan(float(ra)) or math.isnan(float(dec)):
                ra = 'x'
                dec = 'x'
        return ra,dec


    def returnGalactic(self):
        # find and return galactic position

        galacticLong = 'x'
        galacticLat = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_PositionDataTable")
            galacticLong = table.array['pos_lon_gal_d'][0]
            galacticLat = table.array['pos_lat_gal_d'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return galactic position. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if isNumber(galacticLong) or isNumber(galacticLat):
            if math.isnan(float(galacticLong)) or math.isnan(float(galacticLat)):
                galacticLong ='x'
                galacticLat = 'x'
        return galacticLong,galacticLat


    def returnRedIndependentDist(self):
        # find and return redshift-independent distance with mean, std. dev., min and max

        mean_metricDist = 'x'
        std_metricDist = 'x'
        min_metricDist = 'x'
        max_metricDist = 'x'

        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Redshift_IndependentDistances")
#             method = table.array['Statistical method'][0]
            mean_metricDist  = table.array['MetricDistance'][0]
            std_metricDist  = table.array['MetricDistance'][1]
            min_metricDist  = table.array['MetricDistance'][2]
            max_metricDist  = table.array['MetricDistance'][3]
#             median_metricDist  = table.array['MetricDistance'][4]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to find redshift-independent distance. Here's the error "\
                "message built into the exception: \n %s\n" %e)


        if isNumber(mean_metricDist):
            if math.isnan(float(mean_metricDist)):
                mean_metricDist = 'x'
        if isNumber(std_metricDist):
            if math.isnan(float(std_metricDist)):
                std_metricDist = 'x'
        if isNumber(min_metricDist):
            if math.isnan(float(min_metricDist)):
                min_metricDist = 'x'
        if isNumber(max_metricDist):
            if math.isnan(float(max_metricDist)):
                max_metricDist = 'x'
        return mean_metricDist, std_metricDist, min_metricDist, max_metricDist


    def returnMorphology(self):
        # find and return NED homogenized galaxy morphology

        morphology = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            morphology = table.array['morph_type'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return galaxy morphology. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if morphology == '' or morphology == ' ' or morphology == '\n' or morphology == ' \n':
            morphology = 'x'
        elif isNumber(morphology):
            if math.isnan(float(morphology)):
                morphology = 'x'
        return morphology


    # def returnDistIndicator(self):
    #     # find and return distance indicator (under Classifications)
    #
    #     distIndicator = 'x'
    #     try:
    #         for line in self.fullHtml:
    #             index = line.find('Distance Indicator')
    #             if index != -1:
    #                 index2 = line[index:].find('<TD>')
    #                 index3 = line[index2+index:].find('</TD>')
    #                 distIndicator = line[index2+index+4:index3+index+index2]
    #                 break
    #
    #     except Exception as e:
    #         sys.stderr.write("\n Unable to return distance indicator. Here is the error message "\
    #         "built into the exception:\n %s\n" %e)
    #
    #     if isNumber(distIndicator):
    #         if math.isnan(float(distIndicator)):
    #             distIndicator = 'x'
    #     return distIndicator



    def returnDistanceIndicator(self):
        # find and return distance indicator (under Classifications)

        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Classifications")
            # --- check if "Distance Indicator" is present
            ind = np.where(table.array['class_col1'].data == b'Distance Indicator')[0]
            if ind.any():
                # --- class_col3  is the NED Homogenized classification
                # --- decode from bytes type to unicode string
                distance_indicator = table.array['class_col3'].data[ind[0]].decode('UTF-8')
                distance_indicator_ref = table.array['class_col5'].data[ind[0]].decode('UTF-8')
            else:
                distance_indicator = 'x'
                distance_indicator_ref = 'x'

            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return Distance Indicator. Here is the error message built into the exception:\n %s\n" %e)
            distance_indicator = 'x'
            distance_indicator_ref = 'x'

        # --- try getting it from the HTML if the XML/VOTable fails
        if distance_indicator == 'x':
            try:
                for line in self.fullHtml:
                    index = line.find('Distance Indicator')
                    if index != -1:
                        index2 = line[index:].find('<td><strong>')
                        index3 = line[index2+index:].find('<td><strong>')
                        distance_indicator = line[index2 + index + len('<td><strong>'):index3 + index + index2]

                        index_ref = line[index:].find('"ned_dw">')
                        index_ref2 = line[index_ref+index:].find('</a></td></tr>')
                        distance_indicator_ref = line[index + index_ref + len('"ned_dw">'):index + index_ref + index_ref2]
                        break

            except Exception as e:
                sys.stderr.write("\n Unable to return distance indicator. Here is the error message built into the exception:\n %s\n" %e)

        if isNumber(distance_indicator):
            if math.isnan(float(distance_indicator)):
                distance_indicator = 'x'
        return distance_indicator, distance_indicator_ref



    def returnLuminosityClass(self):
        # find and return luminosity class (under Classifications)
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Classifications")
            # --- check if "Distance Indicator" is present
            ind = np.where(table.array['class_col1'].data == b'Luminosity Class')[0]
            if ind.any():
                # --- class_col3  is the NED Homogenized classification
                # --- decode from bytes type to unicode string
                luminosity_class = table.array['class_col3'].data[ind[0]].decode('UTF-8')
                luminosity_class_ref = table.array['class_col5'].data[ind[0]].decode('UTF-8')
            else:
                luminosity_class = 'x'
                luminosity_class_ref = 'x'

            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return Luminosity Class. Here is the error message built into the exception:\n %s\n" %e)
            luminosity_class = 'x'
            luminosity_class_ref = 'x'

        # --- try getting it from the HTML if the XML/VOTable fails
        if luminosity_class == 'x':
            try:
                for line in self.fullHtml:
                    index = line.find('Luminosity Class')
                    if index != -1:
                        index2 = line[index:].find('<td><strong>')
                        index3 = line[index2+index:].find('<td><strong>')
                        luminosity_class = line[index2 + index + len('<td><strong>'):index3 + index + index2]

                        index_ref = line[index:].find('"ned_dw">')
                        index_ref2 = line[index_ref+index:].find('</a></td></tr>')
                        luminosity_class_ref = line[index + index_ref + len('"ned_dw">'):index + index_ref + index_ref2]
                        break

            except Exception as e:
                sys.stderr.write("\n Unable to return Luminosity Class. Here is the error message built into the exception:\n %s\n" %e)

        return luminosity_class, luminosity_class_ref

    #
    # def returnLuminosityClass(self):
    #     # find and return luminosity class (under Classifications)
    #
    #     luminosityClass = 'x'
    #     try:
    #         for line in self.fullHtml:
    #             index = line.find('Luminosity Class')
    #             if index != -1:
    #                 index2 = line[index:].find('<TD>')
    #                 index3 = line[index2+index:].find('</TD>')
    #                 index4 = line[index3+index2+index:].find('<TD>')
    #                 index5 = line[index4+index3+index2+index:].find('</TD>')
    #                 luminosityClass = line[index4+index3+index2+index+4:index5+index4+index3+index+index2]
    #                 break
    #
    #     except Exception as e:
    #         sys.stderr.write("\n Unable to return luminosity class. Here is the error message "\
    #         "built into the exception:\n %s\n" %e)
    #
    #     if isNumber(luminosityClass):
    #         if math.isnan(float(luminosityClass)):
    #             luminosityClass = 'x'
    #     return luminosityClass


    def returnEBminusV(self):
        # find and return foreground galactic extinction E(B-V)

        EBminusV = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            EBminusV = table.array['gal_extinc_E(B-V)'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return EBminusV. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if EBminusV ==' ' or EBminusV == '' or not isNumber(EBminusV):
            EBminusV = 'x'
        elif math.isnan(float(EBminusV)):
            EBminusV = 'x'
        return EBminusV


    def returnRadialVelocity(self):
        # find and return radial velocity

        radialVelocity = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            radialVelocity = table.array['main_col6'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return Radial Velocity. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if radialVelocity == ' ' or radialVelocity == '' or not isNumber(radialVelocity):
            radialVelocity = 'x'
        elif math.isnan(float(radialVelocity)):
            radialVelocity = 'x'
        return radialVelocity


    def returnDiameters(self):
        # find and return major and minor diameters

        majorDiameter = 'x'
        minorDiameter = 'x'
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            majorDiameter = table.array['diam_maj'][0]
            minorDiameter = table.array['diam_min'][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return major or minor diameter. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if majorDiameter == ' ' or majorDiameter == '' or not isNumber(majorDiameter):
            majorDiameter = 'x'
        if minorDiameter == ' ' or minorDiameter == '' or not isNumber(minorDiameter):
            minorDiameter == 'x'
        if isNumber(majorDiameter):
            if math.isnan(float(majorDiameter)):
                majorDiameter = 'x'
        if isNumber(minorDiameter):
            if math.isnan(float(minorDiameter)):
                minorDiameter = 'x'
        return majorDiameter,minorDiameter


    def returnPhotometry(self):
        # find and return quick-look photometry information

        photometry = []
        try:
            k = 0
            lineNumber = 0
            found = False
            for line in self.fullHtml:
                k +=1
                index = line.find('Apparent Mag or Flux')
                if index != -1 and found == False:
                    lineNumber = k
                    found = True
                if found and k >= lineNumber +1:
                    if line.find('<b>NOTE: </b>The above') != -1:
                        break
                    else:
                        indexBand = line.find('<td>')+4
                        band = line[indexBand:indexBand+2]
                        for i in ['U ','B ','V ','R ','H ','NU','U_','B_','V_','R_','H_','b_']:
                            if band == i:
                                indexBandEnd = line[indexBand:].find('</td>')+indexBand
                                fullBand = line[indexBand:indexBandEnd].replace(' ','')

                                indexApparentMag = line[indexBandEnd:].find('<td>')+4 + indexBandEnd
                                indexApparentMagEnd = line[indexApparentMag:].find('</td>')+indexApparentMag
                                fullApparentMag = line[indexApparentMag:indexApparentMagEnd].replace(' ','')

#                               indexRef = line[indexApparentMagEnd:].find('<td>')+4 + indexApparentMagEnd
#                               indexRefEnd = line[indexRef:].find('</td>') + indexRef
                                indexRef = line.find('TARGET="ned_dw">') +16
                                indexRefEnd = line[indexRef:].find('</a>') + indexRef
#                               fullRef = line[indexRef:indexRefEnd].replace(' ','')
                                fullRef = line[indexRef:indexRefEnd]

                                indexAbsoluteMag = line[indexRefEnd:].find('<td>')+4 + indexRefEnd
                                indexAbsoluteMagEnd = line[indexAbsoluteMag:].find('</td>')+indexAbsoluteMag
                                fullAbsoluteMag = line[indexAbsoluteMag:indexAbsoluteMagEnd].replace(' ','')

                                indexBolometric = line[indexAbsoluteMagEnd:].find('<td>')+4 + indexAbsoluteMagEnd
                                indexBolometricEnd = line[indexBolometric:].find('</td>')+indexBolometric
                                fullBolometric = line[indexBolometric:indexBolometricEnd].replace(' ','')

                                photometry.append((fullBand,fullApparentMag,fullAbsoluteMag,\
                                    fullBolometric,fullRef))

        except Exception as e:
            sys.stderr.write("Unable to return photometry data. Here is the error message "\
            "built into the exception:\n %s\n" %e)

        if len(photometry) == 0:
            photometry.append('x')
        return photometry


    def returnNames(self):
        # find and return a list of all the names for a given object

        formattedNames = []
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_NamesTable")
            names = table.array['name_col1']
            parsedNames = parseGalaxyNames(names)
            for name in parsedNames:
                # formattedName = urllib.unquote_plus(name).replace('\n','').strip()
                formattedName = unquote_plus(name).replace('\n','').strip()
                formattedNames.append(formattedName)
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write("\n Unable to return alternative names. Here is the error message "\
            "built into the exception:\n %s\n" %e)
            return 'x'
        if len(formattedNames) ==0:
            formattedNames.append('x')
        return formattedNames


def parseGalaxyNames(nameList):
    # format galaxy names for the url

    newNameList = []
    for name in nameList:
        # print("type(name): ", type(name))
        if not isinstance(name, str):
            name = name.decode('UTF-8')

        nname = name.strip()
        nname = nname.replace('*','')
        nname = quote_plus(nname)
        nname = nname.replace('\n','')
        newNameList.append(nname)
    return newNameList


def createCSVTable(outFile,fieldnames):
    # creates and returns a DictReader object populated with header names

    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n,n) for n in fieldnames)
    writer.writerow(headers)
    return writer



def pickPreferredName(altNames, oldName):
    '''
    Chooses the best name for a galaxy

    Inputs
    -----
        altNames  : array-like
            list of all known names for a galaxy

        oldName   : string
            the original name supplied in the input file or on run. May not be in altNames

    Returns
    ------
        finalName  : string
            the chosen best name  for the galaxy
    '''

    order = ['NGC','IC','MRK','UGC','UGCA','PHL','3C','SBS','MCG','ESO','TON','TONS',\
    'PGC','PG','PB','FGC','HS','HE','KUG','IRAS','RX','CGCG','KAZ','FCC','FAIRALL',\
    'HOLM','IZw','IIZw','IIIZw','IVZw','VZw','VIZw','VIIZw','VIIIZw','IRAS','IRASF',\
    'KISS','KISSR','FBQS','LBQS','PKS','SDSS','VCC','2MASS','2DF','6DF','HIPASS','2MASX']

    # --- add oldName to the list of alternate names if it is not already there
    if len(altNames) >=1:
        if not bfind(str(altNames),str(oldName)):
            try:
                altNames.append(oldName)
            except Exception as e:
                sys.stdout.write("Issue with picking preferred name: {0}".format(e))
                sys.stdout.write("altNames list was: {0}".format(altNames))
    else:
        altNames = [oldName]

    found = False
    for i in order:
        for n in altNames:
            # n is input string, i is what i'm looking for
            if bfind(n,i) and not found and not bfind(n,':') and not bfind(n,'['):
                finalName = n
                found = True
                break

    if not found:
        finalName = oldName

    return finalName


def returnGalaxyName(ra, dec, radius):
    '''
    Inputs
    ------
        ra : string
            right ascension in sexagesimal for a galaxy

        dec : string
            declination in sexagesimal for a galaxy

        radius : float
            radius of search around that point

    Returns
    -------
        name  : string
            returns the NED preferred name for the first (closest) object in the
            search radius results
    '''

    # format RA and Dec
    newra = ra.replace('+','%2B')
    newra = newra.replace('-','%2D')
    newdec = dec.replace('+','%2B')
    newdec = newdec.replace('-','%2D')

    mainhost = 'http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon={0}&lat={1}&radius=0.5&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_main&zv_breaker=30000.0&list_limit=5&img_stamp=YES'.format(newra,newdec)
    try:
        mainurl = urlopen(mainhost)
        print('opened mainurl',mainurl)
    except Exception as e:
        sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n %s\n" %e)

    warnings.simplefilter("ignore")
    mainvotable = 'x'
    try:
        mainvotable = parse(mainurl,pedantic=False)
        print('parsed mainvotable',mainvotable)
    except Exception as e:
        sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n %s \n" %e)

    warnings.resetwarnings()
    mainurl.close()

    name = False
    try:
        warnings.simplefilter("ignore")
        maintable = mainvotable.get_table_by_id("NED_MainTable")
        mainName = maintable.array['main_col2'][0]
        mainName = str(mainName).strip()
        name = quote_plus(mainName)
        name = name.replace('\n','')
        print('found name: ',name)
        warnings.resetwarnings()
    except Exception as e:
        sys.stderr.write("\n Unable to return alternative names. Here is the error message "\
        "built into the exception:\n %s\n" %e)
        name = False

    return name


def is_file(filename):
    # decides if opts.filename is a file or the name of a galaxy

    try:
        openFile = open(filename,'rU')
    except Exception as e:
        openFile = False

    return openFile



def print_output(fieldnames, data):
    '''
        Prints out data associated with each fieldname
    '''
    print('---------------')
    for f,d in zip(fieldnames,data):
        print('{0}: {1}'.format(f,d))

    print()




##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################



def main(opts):
    # assuming 'theFile' contains one name per line, read the file

#     hubbleConstant = 71.
    hubbleConstant = params.hubble_constant


    # set a limit to how many queries to make?
#     max_retrieve = 10000
    max_retrieve = params.max_retrieve

    start = time.time()

    # open the file if it's a file, otherwise opts.filename is the name of a galaxy
    nameList = []
    coordList = []
    if not is_file(opts.filename):
        name = opts.filename
        nameList.append(name)
        coordList.append(('x','x'))

    else:
    #     fileLines = csv.DictReader(theFile,delimiter='|')
        fileLines = csv.DictReader(opts.filename)
        nameList = []
        coordList = []

        for i in fileLines:
            name = i['Object Name']
            ra = i['RA(deg)']
            dec = i['DEC(deg)']

            nameList.append(name)
            coordList.append([ra,dec])

        opts.filename.close()

    # format the names in the list and return a new list
    newNameList = parseGalaxyNames(nameList)

    # ---
    fieldnames = params.fieldnames

    # check to see if file already exists
    full_path = opts.output+'/'+opts.outname+'.csv'
    print('full_path: ',full_path)
    print()
    if os.path.exists(full_path):
        # It exists. Ask what to do.
        answer = input("%s already exists. Append results to this file? [y,n]" %opts.outname)
        while answer != 'y' and answer != 'n':
            answer = input("Please answer with a 'y' or 'n': ")
        if answer == 'y':
            # open that file and read it as a dictionary csv, allowing it to be written to
            readerFile = open(full_path,'rU')
            print('opening: ',full_path)
            reader = csv.DictReader(readerFile)
            writerOutFile = open(opts.output+opts.outname+'1'+'.csv','wt')
            opts.outname = opts.outname+'1'
            writer = createCSVTable(writerOutFile,fieldnames)

            # also determine if the last name written to this file shows up in the
            # names file. If so, assume that we should start from there and ignore
            # all names coming before. Else, start from the top
            oldName = 'x'
            for line in reader:
                oldName = line['oldName']
                lineList = [line['preferredName'],\
                line['oldName'],\
                line['redshift'],\
                line['degreesJ2000RA_Dec'],\
                line['J2000RA_Dec'],\
                line['galacticLong_Lat'],\
                line['rIndependentDistMean_sd_min_max (Mpc)'],\
                line['morphology'],\
                line['distanceIndicator'],\
                line['luminosityClass'],\
                line['EBminusV'],\
                line['radialVelocity (km/s)'],\
                line['vcorr (km/s)'],\
                line['angDiameters (arcsec)'],\
                line['linDiameters (kpc)'],\
                line['distvcorr (Mpc)'],\
                line['inclination (deg)'],\
                line['photometry'],\
                line['alternativeNames']]
                lineRow = dict((f,o) for f,o in zip(fieldnames,lineList))
                writer.writerow(lineRow)

            readerFile.close()
            print('Last entry in {0}: {1}'.format(opts.outname, oldName))
            try:
                nameListIndex = 0
                index = 0
                for n in newNameList:
                    n = unquote_plus(n).replace('\n','').replace(' ','').strip()
                    n = n.replace('*','')
                    if n == oldName:
                        nameListIndex = index
                        break
                    else:
                        index +=1
#                 nameListIndex = newNameList.index(oldName)
            except Exception as e:
                # name not found so an exception is raised
                nameListIndex = False
                print('e: ',e)

            print('nameListIndex: ',nameListIndex)

            if nameListIndex:
                print("len - index: ",len(newNameList),', ',nameListIndex)
                if len(newNameList) - nameListIndex > 2:
                    print('Starting from {0} in object name list'.format(newNameList[nameListIndex+1:nameListIndex+2][0]))
                    coordList = coordList[nameListIndex+1:]
                    newNameList = newNameList  [nameListIndex+1:]
                else:
                    print('Object {0} is the last object in the name list. Please rerun with a'.format(oldName))
                    print('new list of objects to search for. Exiting...')
                    sys.exit()

            else:
                print('Could not find {0} in object name list. Starting from the beginning.'.format(oldName))

        if answer == 'n':
            # ask for a new filename, or give the option of quitting
            aTwo = input("Please enter new filename, or 'q' to quit: ")
            if aTwo == 'q':
                # quit the program
                sys.exit()
            else:
                # create a new file based on the newly entered file name
                opts.outname = aTwo
                writerOutFile = open(opts.output+opts.outname+'.csv','wt')
                writer = createCSVTable(writerOutFile,fieldnames)

    else:
        # It does not already exist. Create it as a dictionary csv file
        writerOutFile = open(full_path,'wt')
        writer = createCSVTable(writerOutFile,fieldnames)

##########################################################################################
##########################################################################################

    total = len(newNameList)
    counter = max_retrieve
    for coord,name in zip(coordList,newNameList):
        # return basic information for each object

        totalStart = time.time()
        percentComplete = round(float((max_retrieve-counter))/total *100,1)
        sys.stdout.write("\n")
        sys.stdout.write("Percent complete: %s %% \n" %percentComplete)
        sys.stdout.write("Galaxies left: %s \n" %counter)
        sys.stdout.write("\n")
        sys.stdout.write("Starting: %s \n" %name)

        start2 = time.time()
        queryTime = time.time() - start2

#         galaxy = galaxyClass(name,votable,html,fullHtml)
        galaxy = galaxyClass(name)
        galaxy.queryNED()
        sys.stdout.write("1...")
        sys.stdout.flush()

        start3 = time.time()
        redshift = galaxy.returnRedshift()
        sys.stdout.write("2...")
        sys.stdout.flush()

        degreesJ2000RA,degreesJ2000Dec = galaxy.returnJ2000Position()
        if isNumber(degreesJ2000RA) and isNumber(degreesJ2000Dec):
            J2000RA_Dec = convertRAandDec(degreesJ2000RA,degreesJ2000Dec,'sexagesimal')
        else:
            J2000RA_Dec = ('x','x')
            degreesJ2000RA,degreesJ2000Dec = coord

        sys.stdout.write("3...")
        sys.stdout.flush()

        galacticLong_Lat = galaxy.returnGalactic()
        sys.stdout.write("4...")
        sys.stdout.flush()

        rIndependentDistMean_sd_min_max = galaxy.returnRedIndependentDist()
        sys.stdout.write("5...")
        sys.stdout.flush()

        morphology = galaxy.returnMorphology()
        sys.stdout.write("6...")
        sys.stdout.flush()

        distance_indicator,  distance_indicator_ref = galaxy.returnDistanceIndicator()
        sys.stdout.write("7...")
        sys.stdout.flush()

        luminosity_class, luminosity_class_ref = galaxy.returnLuminosityClass()

        EBminusV = galaxy.returnEBminusV()

        radialVelocity = galaxy.returnRadialVelocity()
        if isNumber(radialVelocity):
            if float(radialVelocity) >=0:
                vcorr = calculatevcorr(degreesJ2000RA,degreesJ2000Dec,radialVelocity)
                if vcorr >0:
                    distvcorr = vcorr / hubbleConstant
                else:
                    distvcorr ='x'
            else:
                vcorr = 'x'
                distvcorr = 'x'
        else:
            vcorr = 'x'
            distvcorr = 'x'


        # I am assuming that diameters are found in arcmin, and I'm converting to arcsec
        diameters = galaxy.returnDiameters()
        angMaj = diameters[0]
        if isNumber(angMaj):
            angMaj = float(angMaj)*60
        angMin = diameters[1]
        if isNumber(angMin):
            angMin = float(angMin)*60
            # switch them if the major diameter is larger than the minor
            if isNumber(angMaj):
                if angMaj < angMin:
                    angMaj,angMin = angMin, angMaj
        else:
            angMin = 'x'

        if isNumber(angMaj) and isNumber(angMin) and isNumber(distvcorr):
            linDiameters = calculateLinearDiameters(angMaj,angMin,distvcorr)
        else:
            linDiameters = ('x','x')

        sys.stdout.write("8...")
        sys.stdout.flush()

        if isNumber(angMaj) and isNumber(angMin):
            inclination = calculateInclination(angMaj,angMin)
        else:
            inclination = 'x'

#       photometry = galaxy.returnPhotometry()
        photometry = ['x']
        sys.stdout.write("9...")
        sys.stdout.flush()

        alternativeNames = galaxy.returnNames()
        sys.stdout.write("10")
        sys.stdout.flush()

        # remove weird URL coding from name, pick preferred name
        oldName = unquote_plus(name).replace('\n','').strip()
        preferredName = pickPreferredName(alternativeNames,oldName)


        strippedAlternativeNames = []
        for alternate in alternativeNames:
            stripped = alternate.replace(' ','')
            strippedAlternativeNames.append(stripped)

        print('Query time: ',queryTime)
        print('Time for other stuff: ',time.time() - start3)


        # --- create DataFrame and write it out

        df = pd.DataFrame()



        objectInfoList = [preferredName.replace(' ',''),\
        oldName.replace(' ',''),\
        redshift,\
        (float(degreesJ2000RA),float(degreesJ2000Dec)),\
        J2000RA_Dec,galacticLong_Lat,\
        rIndependentDistMean_sd_min_max,\
        morphology,\
        distance_indicator,\
        luminosity_class,\
        EBminusV,\
        radialVelocity,\
        vcorr,\
        (angMaj,angMin),\
        linDiameters,\
        distvcorr,\
        inclination,\
        photometry,\
        strippedAlternativeNames]

        # write info to file
        row = dict((f,o) for f,o in zip(fieldnames,objectInfoList))
        if opts.verbose:
            print('\n','row: ',row)

        writer.writerow(row)
        sys.stdout.write("\n")
        print('Total retrieval time: ',time.time() - totalStart)
        counter -=1
        if counter == 0:
            print('Reached max_retrieve = {0}. Exiting...'.format(max_retrieve))
            break

    writerOutFile.close()
    if opts.verbose:
        print("Done.")
        print()
        print("Results can be found in: {0}".format(full_path))
        print()
        print("Elapsed Time: {0}".format(time.time() - start))
        print()

###############################################################################

if __name__=="__main__":
    __package__ = "expected.package.name"

    # parse commandline
    commandlineOptions = parse_commandline()
    # do the work
    main(commandlineOptions)
