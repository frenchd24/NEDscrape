"""
NEDscrape Navigator

A program to grab galaxy data from the NED server by means of xml-based VOTable
distributions and html parsing, and recompile specified data into a new csv file(s)

"""
__author__ = "David M. French - frenchd24@gmail.com"
__version__ = "3.1"
__date__ = "04/06/2021"

import getpass
import csv
import math
import optparse
import os
from queue import Queue
# import string
import sys
import time

from timeit import default_timer as timer

import threading
import warnings
import re   # regular expression

# import urllib2
from urllib.request import urlopen, urlretrieve
from urllib.parse import quote, quote_plus, unquote_plus
# from urllib.parse import unquote_plus

from astropy.io.votable import parse, tree
import numpy as np
import pandas as pd


################################################################################
# --- Now program specific packages

# --- these are the parameters for this code, should be in the same directory
# from . import params
import params

# --- some useful tools for this program
from utilities import *

# --- functions for dealing with photometry choices
import photometry

################################################################################
# --- Parse the options given on the command-line.

def parse_commandline():
    parser = optparse.OptionParser(
        usage=__doc__, version="$Id: NEDscrape, v {} {}".format(__version__, __date__))
    parser.add_option(
        "-f",
        "--filename",
        dest="filename",
        action="store",
        type="string",
        help="The name a galaxy or of a file containing a list of "
             "galaxy names to be queried")
    parser.add_option(
        "-o",
        "--outdir",
        dest="outdir",
        action="store",
        help="Specify the directory in which to save newly generated csv table")
    parser.add_option(
        "-n",
        "--outfile",
        dest="outfile",
        action="store",
        help="Specify the output filename for the newly generated csv table")
    parser.add_option(
        "-v", "--verbose", action="store_true", default=False, help="Run verbosely")
    parser.add_option(
        '-m', '--multithreading', action='store_true', default=False,
        help='Turn on multithreading for faster data retrieval')

    opts, args = parser.parse_args()

    # --- test for input filename/targetname
    # --- exit and print usage message if not
    if opts.filename is None:
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
        print(" -m tells it to use multithreading to speed up data retrieval. Not for large queries!")
        print(" -v tells it to run verbosely")
        print()
        sys.exit(1)

    # --- if no output filename, create one based on input filename
    if opts.outfile is None:
        if bfind(opts.outfile, '.'):
            opts.outfile = opts.filename[ : opts.filename.find('.')] + "_results"
        else:
            opts.outfile = opts.filename + '_results'

    # --- set output to current directory if not specified
    if opts.outdir is None:
        opts.outdir = os.getcwd() + "/"


    if opts.verbose:
        print()
        print("********************** PARAMETERS ******************************")
        print("Filename: ", opts.filename)
        print("Output directory: ", opts.outdir)
        print("Output filename: ", opts.outfile)
        print()
        print("************************ RESULTS *******************************")

    return opts

################################################################################
# --- class for multi-threading URL retrieval

class ThreadUrl(threading.Thread):
    """Threaded Url Grab"""

    def __init__(self, queue, outQueue, numberQueue):
        threading.Thread.__init__(self)
        self.queue = queue
        self.outQueue = outQueue
        self.numberQueue = numberQueue

    def run(self):
        while True:
            # grabs host from queue
            host = self.queue.get()

            # grabs urls of hosts and reads or parses them according to their type
            try:
                # url = urllib2.urlopen(host[0])
                # url = urlopen(host[0])
                url_obj, headers = urlretrieve(host[0])

            except Exception as e:
                sys.stderr.write("\n Unable to return url file object. Here is the error message built into the exception: \n {}\n".format(e))

            if host[1] == 1 or host[1] == 4:
                warnings.simplefilter("ignore")
                votable1 = "x"
                try:
                    # votable1 = parse(url, pedantic=False)
                    votable1 = parse(url_obj, pedantic=False)
                except Exception as e:
                    sys.stderr.write("\n Unable to parse voTable. Here is the message built into the exception: \n {} \n".format(e))

                self.outQueue.put(votable1)
                warnings.resetwarnings()

            else:
                html = "x"
                try:
                    # html = url.readlines()
                    with open(url_obj) as f:
                        html = f.readlines()
                except Exception as e:
                    sys.stderr.write("\n Unable to read html. Here is the message built into the exception: \n {} \n".format(e))

                self.outQueue.put(html)

            self.numberQueue.put(host[1])
            self.queue.task_done()
            # url.close()

################################################################################
# --- class for retrieving the specific galaxy info I want

class galaxyClass(object):
    # a galaxy object containing it's information and methods of retrieving and
    # returning it

    #     def __init__(self,name,votable,html,fullHtml):
    #         self.name = name
    #         self.votable = votable
    #         self.html = html
    #         self.fullHtml = fullHtml

    def __init__(self, name):
        self.name = name


    def queryNED(self, multithreading=False):
        """
        query NED server for each object name, returning the file object representing
        the xml VOTable output from NED as well as the redshift-independent distance html file

        Note: htmlFileObject is the html page for redshift independent distance, while
        fullHtmlFileObject is the whole html page for the object

        """
        #         queue = Queue.Queue()
        #         outQueue = Queue.Queue()
        #         voList = []

        # queue = Queue.Queue()
        # outQueue = Queue.Queue()
        # numberQueue = Queue.Queue()
        queue = Queue()
        outQueue = Queue()
        numberQueue = Queue()

        t1 = ThreadUrl(queue, outQueue, numberQueue)
        t2 = ThreadUrl(queue, outQueue, numberQueue)
        t3 = ThreadUrl(queue, outQueue, numberQueue)
        t4 = ThreadUrl(queue, outQueue, numberQueue)
        t1.setDaemon(True)
        t2.setDaemon(True)
        t3.setDaemon(True)
        t4.setDaemon(True)
        t1.start()
        t2.start()
        t3.start()
        t4.start()

        name = self.name

        # --- maybe use this instead? (03/11/18)
        #        https://ned.ipac.caltech.edu/cgi-bin/objsearch?objname=%s&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=xml_all&zv_breaker=30000.0&list_limit=5&img_stamp=YES

        # --- host #1 = main VOTable
        # --- host #2 = redshift independent distances
        # --- host #3 = full HTML result
        # --- host #4 = photometry
        hosts = [
            ["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?extend=no&of=xml_all&objname={}".format(
                name), 1,],
            ["http://ned.ipac.caltech.edu/cgi-bin/nDistance?name={}".format(
                name), 2],
            ["http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname={}&extend=no&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES".format(
                name), 3,],
            ["http://ned.ipac.caltech.edu/cgi-bin/nph-datasearch?objname={}&meas_type=bot&ebars_spec=ebars&label_spec=no&x_spec=freq&y_spec=Fnu_jy&xr=-1&of=xml_main&search_type=Photometry".format(
                name), 4,],
                ]

        # --- if the multithreading options is on, add all the hosts
        # --- to the queue all at once
        if multithreading:
            # --- populate queue with data
            for host in hosts:
                #             time.sleep(0.1)
                queue.put(host)

            # --- wait on the queue until everything has been processed
            queue.join()

            # --- decide which object is which in the queue
            first = outQueue.get()
            second = outQueue.get()
            third = outQueue.get()
            fourth = outQueue.get()

            number1 = numberQueue.get()
            number2 = numberQueue.get()
            number3 = numberQueue.get()
            number4 = numberQueue.get()

            if number1 == 1:
                votable = first
            elif number1 == 2:
                html = first
            elif number1 == 3:
                fullHtml = first
            elif number1 == 4:
                photometry_votable = first

            if number2 == 1:
                votable = second
            elif number2 == 2:
                html = second
            elif number2 == 3:
                fullHtml = second
            elif number2 == 4:
                photometry_votable = second

            if number3 == 1:
                votable = third
            elif number3 == 2:
                html = third
            elif number3 == 3:
                fullHtml = third
            elif number3 == 4:
                photometry_votable = third

            if number4 == 1:
                votable = fourth
            elif number4 == 2:
                html = fourth
            elif number4 == 3:
                fullHtml = fourth
            elif number4 == 4:
                photometry_votable = fourth

        # --- if no multithreading, then add and remove each host one at a time
        else:
            # --- populate queue with data
            for host in hosts:
                # request_time = time.start()
                request_time = timer()
                queue.put(host)

                # --- wait on the queue until everything has been processed
                queue.join()

                # --- decide which object is which in the queue
                if host[1] == 1:
                    first = outQueue.get()
                    number1 = numberQueue.get()

                    if number1 == 1:
                        votable = first
                    else:
                        print('Problem: number1 != 1')

                elif host[1] == 2:
                    second = outQueue.get()
                    number2 = numberQueue.get()

                    if number2 == 2:
                        html = second
                    else:
                        print('Problem: number2 != 2')

                elif host[1] == 3:
                    third = outQueue.get()
                    number3 = numberQueue.get()

                    if number3 == 3:
                        fullHtml = third
                    else:
                        print('Problem: number3 != 3')

                elif host[1] == 4:
                    fourth = outQueue.get()
                    number4 = numberQueue.get()

                    if number4 == 4:
                        photometry_votable = fourth
                    else:
                        print('Problem: number4 != 4')

                else:
                    print('SOMETHING WENT TERRIBLY WRONG')
                    sys.exit()

                while timer() - request_time < 1.0:
                    time.sleep(0.1)


        self.votable = votable
        self.html = html
        self.fullHtml = fullHtml
        self.photometry_votable = photometry_votable


    def returnRedshift(self):
        # --- find and return redshift for each object

        redshift = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            redshift = table.array["Redshift"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return redshift. Here is the error message "
                "built into the exception:\n {} \n".format(e))

        if isNumber(redshift):
            if math.isnan(float(redshift)):
                redshift = "x"
        return redshift


    def returnJ2000Position(self):
        # --- find and return equatorial J2000 position RA and Dec for each object

        ra = "x"
        dec = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            #             ra = table.array['RA(deg)'][0]
            #             dec = table.array['DEC(deg)'][0]
            ra = table.array["RA"][0]
            dec = table.array["DEC"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return J2000 position. Here is the error message "
                "built into the exception:\n {} \n".format(e))

        if isNumber(ra) or isNumber(dec):
            if math.isnan(float(ra)) or math.isnan(float(dec)):
                ra = "x"
                dec = "x"
        return ra, dec


    def returnGalactic(self):
        # --- find and return galactic position

        galacticLong = "x"
        galacticLat = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_PositionDataTable")
            galacticLong = table.array["pos_lon_gal_d"][0]
            galacticLat = table.array["pos_lat_gal_d"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return galactic position. Here is the error message "
                "built into the exception:\n {} \n".format(e))

        if isNumber(galacticLong) or isNumber(galacticLat):
            if math.isnan(float(galacticLong)) or math.isnan(float(galacticLat)):
                galacticLong = "x"
                galacticLat = "x"
        return galacticLong, galacticLat


    def returnRedIndependentDist(self):
        # --- find and return redshift-independent distance with mean, std. dev., min and max

        mean_metricDist = "x"
        std_metricDist = "x"
        min_metricDist = "x"
        max_metricDist = "x"

        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Redshift_IndependentDistances")
            #             method = table.array['Statistical method'][0]
            mean_metricDist = table.array["MetricDistance"][0]
            std_metricDist = table.array["MetricDistance"][1]
            min_metricDist = table.array["MetricDistance"][2]
            max_metricDist = table.array["MetricDistance"][3]
            #             median_metricDist  = table.array['MetricDistance'][4]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to find redshift-independent distance. Here's the error "
                "message built into the exception: \n {} \n".format(e))

        if isNumber(mean_metricDist):
            if math.isnan(float(mean_metricDist)):
                mean_metricDist = "x"
        if isNumber(std_metricDist):
            if math.isnan(float(std_metricDist)):
                std_metricDist = "x"
        if isNumber(min_metricDist):
            if math.isnan(float(min_metricDist)):
                min_metricDist = "x"
        if isNumber(max_metricDist):
            if math.isnan(float(max_metricDist)):
                max_metricDist = "x"
        return mean_metricDist, std_metricDist, min_metricDist, max_metricDist


    def returnMorphology(self):
        # --- find and return NED homogenized galaxy morphology

        morphology = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            morphology = table.array["morph_type"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return galaxy morphology. Here is the error message "
                "built into the exception:\n {} \n".format(e))

        if (morphology == "" or morphology == " " or morphology == "\n" or morphology == " \n"):
            morphology = "x"
        elif isNumber(morphology):
            if math.isnan(float(morphology)):
                morphology = "x"
        return morphology


    def returnDistanceIndicator(self):
        '''
            Find and return Distance Indicator (under Classifications)

            Returns:
            -------
                distance_indicator        :   str
                    the method used for redshift independent
                    distance measurement (e.g., Tully-Fisher)

                distance_indicator_ref    :   str
                    reference for where the classification came from
        '''
        # --- find and return distance indicator (under Classifications)

        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Classifications")
            # --- check if "Distance Indicator" is present
            ind = np.where(table.array["class_col1"].data == "Distance Indicator")[0]
            if ind.size == 0:
                distance_indicator = "x"
                distance_indicator_ref = "x"
            else:
                # --- class_col3  is the NED Homogenized classification
                try:
                    # --- decode from bytes type to unicode string
                    distance_indicator = (
                        table.array["class_col3"].data[ind[0]].decode("UTF-8"))
                    distance_indicator_ref = (
                        table.array["class_col5"].data[ind[0]].decode("UTF-8"))
                except (UnicodeDecodeError, AttributeError):
                    distance_indicator = (
                        table.array["class_col3"].data[ind[0]])
                    distance_indicator_ref = (
                        table.array["class_col5"].data[ind[0]])

            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return Distance Indicator. Here is the error "
                "message built into the exception:\n {} \n".format(e))
            distance_indicator = "x"
            distance_indicator_ref = "x"

        # --- try getting it from the HTML if the XML/VOTable fails
        if distance_indicator == "x":
            try:
                for line in self.fullHtml:
                    index = line.find("Distance Indicator")
                    if index != -1:
                        index2 = line[index:].find("<td><strong>")
                        index3 = line[index2 + index :].find("<td><strong>")
                        distance_indicator = line[
                            index2
                            + index
                            + len("<td><strong>") : index3
                            + index
                            + index2]

                        index_ref = line[index:].find('"ned_dw">')
                        index_ref2 = line[index_ref + index :].find("</a></td></tr>")
                        distance_indicator_ref = line[
                            index
                            + index_ref
                            + len('"ned_dw">') : index
                            + index_ref
                            + index_ref2]
                        break

            except Exception as e:
                sys.stderr.write("\n Unable to return distance indicator. "
                "Here is the error message built into the exception:"
                "\n {} \n".format(e))

        if distance_indicator == '' or distance_indicator == '\n':
            distance_indicator = 'x'
        if isNumber(distance_indicator):
            if math.isnan(float(distance_indicator)):
                distance_indicator = "x"

        return distance_indicator, distance_indicator_ref


    def returnLuminosityClass(self):
        '''
            Find and return luminosity class (under Classifications)

            Returns:
            -------
                luminosity_class        :   str
                    the luminosity classe (usually Roman numeral)

                luminosity_class_ref    :   str
                    reference for where the classification came from
        '''

        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("Classifications")
            # --- check if "Luminosity Class" is present
            # --- This is saved as bytes type instead of unicode, hence the 'b'
            ind = np.where(table.array["class_col1"].data == "Luminosity Class")[0]
            # print(table.array['class_col1'].data)
            # print('ind: ',ind)
            # print('ind.size == 0?', ind.size == 0)
            # print('ind.type(): ',ind.type())
            if ind.size == 0:
                luminosity_class = "x"
                luminosity_class_ref = "x"
            else:
                # --- class_col3  is the NED Homogenized classification
                try:
                    # --- decode from bytes type to unicode string
                    luminosity_class = (
                        table.array["class_col3"].data[ind[0]].decode("UTF-8"))
                    luminosity_class_ref = (
                        table.array["class_col5"].data[ind[0]].decode("UTF-8"))
                except (UnicodeDecodeError, AttributeError):
                    luminosity_class = (
                        table.array["class_col3"].data[ind[0]])
                    luminosity_class_ref = (
                        table.array["class_col5"].data[ind[0]])

            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return Luminosity Class. Here is the error "
                "message built into the exception:\n {} \n".format(e))
            luminosity_class = "x"
            luminosity_class_ref = "x"

        # --- try getting it from the HTML if the XML/VOTable fails
        if luminosity_class == "x":
            try:
                for line in self.fullHtml:
                    index = line.find("Luminosity Class")
                    if index != -1:
                        index2 = line[index:].find("<td><strong>")
                        index3 = line[index2 + index :].find("<td><strong>")
                        luminosity_class = line[
                            index2
                            + index
                            + len("<td><strong>") : index3
                            + index
                            + index2]

                        index_ref = line[index:].find('"ned_dw">')
                        index_ref2 = line[index_ref + index :].find("</a></td></tr>")
                        luminosity_class_ref = line[
                            index
                            + index_ref
                            + len('"ned_dw">') : index
                            + index_ref
                            + index_ref2]
                        break

            except Exception as e:
                sys.stderr.write(
                    "\n Unable to return Luminosity Class. Here is the error "
                    "message built into the exception:\n {} \n".format(e))
                luminosity_class = 'x'
                luminosity_class_ref = 'x'

        # --- check for empty or null results
        if luminosity_class.strip() == '' or luminosity_class.strip() == '\n':
            luminosity_class = 'x'
            luminosity_class_ref = 'x'

        return luminosity_class, luminosity_class_ref


    def returnEBminusV(self):
        # --- find and return foreground galactic extinction E(B-V)

        EBminusV = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            EBminusV = table.array["gal_extinc_E(B-V)"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return EBminusV. Here is the error message "
                "built into the exception:\n {} \n".format(e))

        if EBminusV == " " or EBminusV == "" or not isNumber(EBminusV):
            EBminusV = "x"
        elif math.isnan(float(EBminusV)):
            EBminusV = "x"
        return EBminusV


    def returnRadialVelocity(self):
        # --- find and return radial velocity

        radialVelocity = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_MainTable")
            radialVelocity = table.array["main_col6"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return Radial Velocity. Here is the error message"
                " built into the exception:\n {} \n".format(e))

        if (radialVelocity == " " or radialVelocity == "" or not isNumber(radialVelocity)):
            radialVelocity = "x"
        elif math.isnan(float(radialVelocity)):
            radialVelocity = "x"
        return radialVelocity


    def returnDiameters(self):
        # --- find and return major and minor diameters

        majorDiameter = "x"
        minorDiameter = "x"
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_BasicDataTable")
            majorDiameter = table.array["diam_maj"][0]
            minorDiameter = table.array["diam_min"][0]
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return major or minor diameter. Here is the error"
                " message built into the exception:\n {} \n".format(e))

        if majorDiameter == " " or majorDiameter == "" or not isNumber(majorDiameter):
            majorDiameter = "x"
        if minorDiameter == " " or minorDiameter == "" or not isNumber(minorDiameter):
            minorDiameter == "x"
        if isNumber(majorDiameter):
            if math.isnan(float(majorDiameter)):
                majorDiameter = "x"
        if isNumber(minorDiameter):
            if math.isnan(float(minorDiameter)):
                minorDiameter = "x"
        return majorDiameter, minorDiameter



    def returnPhotometry(self):
        """
        Find and return photometry for each object.

        Returns a dictionary where the keys are:
            {}'b', 'u', 'g', 'r', 'i', 'z', 'j', 'h', 'k'}
            and the values are the corresponding photometry measurements
        """

        photometry_df = pd.DataFrame()
        warnings.simplefilter("ignore")
        try:
            table = self.photometry_votable.get_table_by_id("NED_PhotometricData")
            passbands = table.array["photo_col2"]
            phot_measurement = table.array["photo_col3"]
            uncertainty = table.array["photo_col4"]
            units = table.array["photo_col5"]
            frequency = table.array["photo_col6"]
            ned_PM = table.array["photo_col7"]
            ned_uncertainty = table.array["photo_col8"]
            ned_units = table.array["photo_col9"]
            refcode = table.array["photo_col10"]
            significance = table.array["photo_col11"]
            published_freq = table.array["photo_col12"]
            freq_mode = table.array["photo_col13"]
            spatial_mode = table.array["photo_col15"]
            qualifiers = table.array["photo_col16"]
            comments = table.array["photo_col17"]

            photometry_df = pd.DataFrame({
                            'passband':passbands,
                            'phot_measurement':phot_measurement,
                            'uncertainty':uncertainty,
                            'units':units,
                            'frequency':frequency,
                            'ned_PM':ned_PM,
                            'ned_uncertainty':ned_uncertainty,
                            'ned_units':ned_units,
                            'refcode':refcode,
                            'significance':significance,
                            'published_freq':published_freq,
                            'freq_mode':freq_mode,
                            'spatial_mode':spatial_mode,
                            'qualifiers':qualifiers,
                            'comments':comments})

        except Exception as e:
            sys.stderr.write(
            "\n Unable to return passband. Here is the error message built "
            "into the exception: \n {} \n".format(e))

        # --- sort through photometry information
        phot_dict = {
            "b": [],
            "u": [],
            "g": [],
            "r": [],
            "i": [],
            "z": [],
            "j": [],
            "h": [],
            "k": []}

        if self.photometry_votable != "x" and len(photometry_df) >0:
            # photometry_df = photometry_df.columns.str.strip()
            photometry_df.passband = photometry_df.passband.apply(str.strip).apply(str.lower)
            # photometry_df.passband = photometry_df.passband.agg('str', 'strip', 'lower')

            # print('photometry_df.head(15): ', photometry_df.head(150))

            b_df = photometry_df[photometry_df.passband.str[0] == 'b']
            j_df = photometry_df[(photometry_df.passband.str[0] == 'j') &
                                 (photometry_df.passband.str.contains('2mass'))]
            h_df = photometry_df[(photometry_df.passband.str[0] == 'h') &
                                 (photometry_df.passband.str.contains('2mass'))]
            k_df = photometry_df[(photometry_df.passband.str[0] == 'k') &
                                 (photometry_df.passband.str.contains('2mass'))]

            u_df = photometry_df[(photometry_df.passband.str[0] == 'u') &
                                 (photometry_df.passband.str.contains('sdss'))]
            g_df = photometry_df[(photometry_df.passband.str[0] == 'g') &
                                 (photometry_df.passband.str.contains('sdss'))]
            r_df = photometry_df[(photometry_df.passband.str[0] == 'r') &
                                 (photometry_df.passband.str.contains('sdss'))]
            i_df = photometry_df[(photometry_df.passband.str[0] == 'i') &
                                 (photometry_df.passband.str.contains('sdss'))]
            z_df = photometry_df[(photometry_df.passband.str[0] == 'z') &
                                 (photometry_df.passband.str.contains('sdss'))]

            phot_dict = {
                "b": b_df,
                "u": u_df,
                "g": g_df,
                "r": r_df,
                "i": i_df,
                "z": z_df,
                "j": j_df,
                "h": h_df,
                "k": k_df}

        return phot_dict



    def returnPhotometry_old(self):
        # --- DEPRECIATED; USE THE ONE ABOVE (USING PANDAS)
        """
        Find and return photometry for each object.

        Returns a dictionary where the keys are:
            {}'b', 'u', 'g', 'r', 'i', 'z', 'j', 'h', 'k'}
            and the values are the corresponding photometry measurements
        """

        totalArray = []
        photometry_df = pd.DataFrame()
        warnings.simplefilter("ignore")
        try:
            table = self.photometry_votable.get_table_by_id("NED_PhotometricData")
            passbands = table.array["photo_col2"]
            phot_measurement = table.array["photo_col3"]
            uncertainty = table.array["photo_col4"]
            units = table.array["photo_col5"]
            frequency = table.array["photo_col6"]
            ned_PM = table.array["photo_col7"]
            ned_uncertainty = table.array["photo_col8"]
            ned_units = table.array["photo_col9"]
            refcode = table.array["photo_col10"]
            significance = table.array["photo_col11"]
            published_freq = table.array["photo_col12"]
            freq_mode = table.array["photo_col13"]
            spatial_mode = table.array["photo_col15"]
            qualifiers = table.array["photo_col16"]
            comments = table.array["photo_col17"]

            totalArray = zip(
                passbands,
                phot_measurement,
                uncertainty,
                units,
                frequency,
                ned_PM,
                ned_uncertainty,
                ned_units,
                refcode,
                significance,
                published_freq,
                freq_mode,
                spatial_mode,
                qualifiers,
                comments)

        except Exception as e:
            sys.stderr.write(
            "\n Unable to return passband. Here is the error message built "
            "into the exception: \n {} \n".format(e))

        if not totalArray:
            totalArray.append("x")

        # --- sort through photometry information
        phot_dict = {
            "b": [],
            "u": [],
            "g": [],
            "r": [],
            "i": [],
            "z": [],
            "j": [],
            "h": [],
            "k": []}

        for measurement in totalArray:
            # --- this is just the name of the measurement,
            # --- e.g., 'B (m_B)' or 'g (SDSS PSF) AB'
            passband = str(measurement[0]).lower().strip()

            # --- now check if the first letter matches any of the specific
            # --- bands I'm looking for
            if passband[0] in phot_dict:
                go = True

                # --- if the first letter is h, j, k but is not '2MASS',
                # --- don't put it in the corresponding bin
                if passband[0] == "h" and not bfind(passband, "2mass"):
                    go = False
                if passband[0] == "j" and not bfind(passband, "2mass"):
                    go = False
                if passband[0] == "k" and not bfind(passband, "2mass"):
                    go = False

                # --- if the first letter is u,g,r,i,z but is not 'sdss',
                # --- don't put it in the corresponding bin
                if passband[0] == "u" and not bfind(passband, "sdss"):
                    go = False
                if passband[0] == "g" and not bfind(passband, "sdss"):
                    go = False
                if passband[0] == "r" and not bfind(passband, "sdss"):
                    go = False
                if passband[0] == "i" and not bfind(passband, "sdss"):
                    go = False
                if passband[0] == "z" and not bfind(passband, "sdss"):
                    go = False

                if go:
                    phot_dict[passband[0]].append(measurement)

        return phot_dict


    def returnNames(self):
        # find and return a list of all the names for a given object

        formattedNames = []
        try:
            warnings.simplefilter("ignore")
            table = self.votable.get_table_by_id("NED_NamesTable")
            names = table.array["name_col1"]
            parsedNames = parseGalaxyNames(names)
            for name in parsedNames:
                # formattedName = urllib.unquote_plus(name).replace('\n','').strip()
                formattedName = unquote_plus(name).replace("\n", "").strip()
                formattedNames.append(formattedName)
            warnings.resetwarnings()

        except Exception as e:
            sys.stderr.write(
                "\n Unable to return alternative names. Here is the error "
                "message built into the exception:\n %s\n".format(e))

            return "x"

        if len(formattedNames) == 0:
            formattedNames.append("x")
        return formattedNames

################################################################################
# --- other helper functions

def parseGalaxyNames(nameList):
    # format galaxy names for the url

    newNameList = []
    for name in nameList:
        # print("type(name): ", type(name))
        if not isinstance(name, str):
            name = name.decode("UTF-8")

        nname = name.strip()
        nname = nname.replace("*", "")
        nname = quote_plus(nname)
        nname = nname.replace("\n", "")
        newNameList.append(nname)
    return newNameList


def createCSVTable(outFile, fieldnames):
    # creates and returns a DictReader object populated with header names

    writer = csv.DictWriter(outFile, fieldnames=fieldnames)
    headers = dict((n, n) for n in fieldnames)
    writer.writerow(headers)
    return writer


def pickPreferredName(altNames, oldName):
    """
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
    """

    order = [
        "NGC",
        "IC",
        "MRK",
        "UGC",
        "UGCA",
        "PHL",
        "3C",
        "SBS",
        "MCG",
        "ESO",
        "TON",
        "TONS",
        "PGC",
        "PG",
        "PB",
        "FGC",
        "HS",
        "HE",
        "KUG",
        "IRAS",
        "RX",
        "CGCG",
        "KAZ",
        "FCC",
        "FAIRALL",
        "HOLM",
        "IZw",
        "IIZw",
        "IIIZw",
        "IVZw",
        "VZw",
        "VIZw",
        "VIIZw",
        "VIIIZw",
        "IRAS",
        "IRASF",
        "KISS",
        "KISSR",
        "FBQS",
        "LBQS",
        "PKS",
        "SDSS",
        "VCC",
        "2MASS",
        "2DF",
        "6DF",
        "HIPASS",
        "2MASX"
        ]

    # --- add oldName to the list of alternate names if it is not already there
    if len(altNames) >= 1:
        if not bfind(str(altNames), str(oldName)):
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
            if bfind(n, i) and not found and not bfind(n, ":") and not bfind(n, "["):
                finalName = n
                found = True
                break

    if not found:
        finalName = oldName

    return finalName


def returnGalaxyName(ra, dec, radius):
    """
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
    """

    # format RA and Dec
    newra = ra.replace("+", "%2B")
    newra = newra.replace("-", "%2D")
    newdec = dec.replace("+", "%2B")
    newdec = newdec.replace("-", "%2D")

    mainhost = "http://ned.ipac.caltech.edu/cgi-bin/objsearch?search_type=Near+Position+Search&in_csys=Equatorial&in_equinox=J2000.0&lon={0}&lat={1}&radius=0.5&hconst=73&omegam=0.27&omegav=0.73&corr_z=1&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=xml_main&zv_breaker=30000.0&list_limit=5&img_stamp=YES".format(
        newra, newdec)
    try:
        mainurl = urlopen(mainhost)
        print("opened mainurl", mainurl)
    except Exception as e:
        sys.stderr.write(
            "\n Unable to return url file object. Here is the error message "
            "built into the exception: \n {}\n".format(e))

    warnings.simplefilter("ignore")
    mainvotable = "x"
    try:
        mainvotable = parse(mainurl, pedantic=False)
        print("parsed mainvotable", mainvotable)
    except Exception as e:
        sys.stderr.write(
            "\n Unable to parse voTable. Here is the message built into "
            "the exception: \n %s \n".format(e))

    warnings.resetwarnings()
    mainurl.close()

    name = False
    try:
        warnings.simplefilter("ignore")
        maintable = mainvotable.get_table_by_id("NED_MainTable")
        mainName = maintable.array["main_col2"][0]
        mainName = str(mainName).strip()
        name = quote_plus(mainName)
        name = name.replace("\n", "")
        print("found name: ", name)
        warnings.resetwarnings()
    except Exception as e:
        sys.stderr.write(
            "\n Unable to return alternative names. Here is the error message "
            "built into the exception:\n {} \n".format(e))
        name = False

    return name


def is_file(filename):
    # decides if opts.filename is a file or the name of a galaxy

    try:
        openFile = open(filename, "rU")
    except Exception as e:
        openFile = False

    return openFile


def print_output(fieldnames, data):
    """
    Prints out data associated with each fieldname
    """
    print("---------------")
    for f, d in zip(fieldnames, data):
        print("{0}: {1}".format(f, d))

    print()


################################################################################
################################################################################
################################################################################
################################################################################


def main(opts):
    # --- define some parameters
    hubble_constant = params.hubble_constant
    B_star = params.B_star
    fieldnames = params.fieldnames

    # --- set a limit to how many queries to make
    max_retrieve = params.max_retrieve

    # --- begin the retrieval timer
    total_start = time.time()

    # --- open the file if it's a file, otherwise opts.filename is the name of a galaxy
    nameList = []
    coordList = []
    if not is_file(opts.filename):
        name = opts.filename
        nameList.append(name)
        coordList.append(("x", "x"))

    else:
        #     fileLines = csv.DictReader(theFile,delimiter='|')
        try:
            fileLines = csv.DictReader(opts.filename)
            nameList = []
            coordList = []

            for i in fileLines:
                name = i["Object Name"]
                ra = i["RA(deg)"]
                dec = i["DEC(deg)"]

                nameList.append(name.replace(' ', ''))
                coordList.append([ra, dec])

            opts.filename.close()
        except KeyError:
            # fileLines = np.loadtxt(opts.filename, unpack=False)
            nameList = []
            coordList = []
            with open(opts.filename) as f:
                for fileLine in f:
                    nameList.append(fileLine.replace(' ',''))
                    coordList.append(['x', 'x'])

        else:
            sys.stdout("Cannot read input file")
            sys.exit()

    # --- format the names in the list and return a new list
    newNameList = parseGalaxyNames(nameList)

    # --- check to see if file already exists
    full_path = opts.outdir + "/" + opts.outfile
    if "csv" not in opts.outfile[-3:]:
        full_path += ".csv"

    print("Output file path: ", full_path)
    print()
    if os.path.exists(full_path):
        # --- File already exists. Ask what to do.
        answer = input(
            "{} already exists. Append results to this file? [y,n]".format(opts.outfile)
        )
        while answer != "y" and answer != "n":
            answer = input("Please answer with 'y' or 'n': ")
        if answer == "y":
            # open that file and read it as a dictionary csv, allowing it to be written to
            readerFile = open(full_path, "r")
            print("opening: ", full_path)
            reader = csv.DictReader(readerFile)
            writerOutFile = open(opts.outdir + opts.outfile + "1" + ".csv", "wt")
            opts.outfile = opts.outfile + "1"
            writer = createCSVTable(writerOutFile, fieldnames)

            # also determine if the last name written to this file shows up in the
            # names file. If so, assume that we should start from there and ignore
            # all names coming before. Else, start from the top
            oldName = "x"
            for line in reader:
                oldName = line["oldName"]
                lineList = [
                    line["preferredName"],
                    line["oldName"],
                    line["redshift"],
                    line["degreesJ2000RA_Dec"],
                    line["J2000RA_Dec"],
                    line["galacticLong_Lat"],
                    line["rIndependentDistMean_sd_min_max (Mpc)"],
                    line["morphology"],
                    line["distanceIndicator"],
                    line["luminosityClass"],
                    line["EBminusV"],
                    line["radialVelocity (km/s)"],
                    line["vcorr (km/s)"],
                    line["angDiameters (arcsec)"],
                    line["linDiameters (kpc)"],
                    line["distvcorr (Mpc)"],
                    line["inclination (deg)"],
                    line["b_phot"],
                    line["u_phot"],
                    line["g_phot"],
                    line["r_phot"],
                    line["i_phot"],
                    line["z_phot"],
                    line["j_phot"],
                    line["h_phot"],
                    line["k_phot"],
                    line["alternativeNames"],
                ]
                lineRow = dict((f, o) for f, o in zip(fieldnames, lineList))
                writer.writerow(lineRow)

            readerFile.close()
            print("Last entry in {0}: {1}".format(opts.outfile, oldName))
            try:
                nameListIndex = 0
                index = 0
                for n in newNameList:
                    n = unquote_plus(n).replace("\n", "").replace(" ", "").strip()
                    n = n.replace("*", "")
                    if n == oldName:
                        nameListIndex = index
                        break
                    else:
                        index += 1
            #                 nameListIndex = newNameList.index(oldName)
            except Exception as e:
                # name not found so an exception is raised
                nameListIndex = False
                print("e: ", e)

            print("nameListIndex: ", nameListIndex)

            if nameListIndex:
                print("len - index: ", len(newNameList), ", ", nameListIndex)
                if len(newNameList) - nameListIndex > 2:
                    print("Starting from {0} in object name list".format(
                            newNameList[nameListIndex + 1 : nameListIndex + 2][0]))

                    coordList = coordList[nameListIndex + 1 :]
                    newNameList = newNameList[nameListIndex + 1 :]
                else:
                    print("Object {0} is the last object in the name list. "
                            "Please rerun with a".format(oldName))

                    print("new list of objects to search for. Exiting...")
                    sys.exit()

            else:
                print("Could not find {0} in object name list. Starting "
                        "from the beginning.".format(oldName))

        if answer == "n":
            # ask for a new filename, or give the option of quitting
            aTwo = input("Please enter new filename, or 'q' to quit: ")
            if aTwo == "q":
                # quit the program
                sys.exit()
            else:
                # create a new file based on the newly entered file name
                opts.outfile = aTwo
                writerOutFile = open(opts.outdir + opts.outfile + ".csv", "wt")
                writer = createCSVTable(writerOutFile, fieldnames)

    else:
        # It does not already exist. Create it as a dictionary csv file
        writerOutFile = open(full_path, "wt")
        writer = createCSVTable(writerOutFile, fieldnames)

    ##########################################################################################
    ##########################################################################################

    total = len(newNameList)
    counter = max_retrieve
    for coord, name in zip(coordList, newNameList):

        # --- Record start time and update counter output
        target_start = time.time()
        percentComplete = round(float((max_retrieve - counter)) / total * 100, 1)
        sys.stdout.write("\n")
        sys.stdout.write("Percent complete: %s %% \n" % percentComplete)
        sys.stdout.write("Galaxies left: %s \n" % counter)
        sys.stdout.write("\n")
        sys.stdout.write("Starting: %s \n" % name)

        # --- Query NED and record time
        query_start = time.time()
        sys.stdout.write("1...")
        sys.stdout.flush()
        galaxy = galaxyClass(name)

        if opts.multithreading:
            galaxy.queryNED(multithreading=True)
        else:
            galaxy.queryNED()

        query_time = time.time() - query_start

        other_start = time.time()
        sys.stdout.write("2...")
        sys.stdout.flush()
        redshift = galaxy.returnRedshift()

        sys.stdout.write("3...")
        sys.stdout.flush()
        degreesJ2000RA, degreesJ2000Dec = galaxy.returnJ2000Position()
        if isNumber(degreesJ2000RA) and isNumber(degreesJ2000Dec):
            J2000RA_Dec = convertRAandDec(
                degreesJ2000RA, degreesJ2000Dec, "sexagesimal"
            )
        else:
            J2000RA_Dec = ("x", "x")
            degreesJ2000RA, degreesJ2000Dec = coord

        sys.stdout.write("4...")
        sys.stdout.flush()
        galacticLong_Lat = galaxy.returnGalactic()

        sys.stdout.write("5...")
        sys.stdout.flush()
        rIndependentDistMean_sd_min_max = galaxy.returnRedIndependentDist()

        sys.stdout.write("6...")
        sys.stdout.flush()
        morphology = galaxy.returnMorphology()

        sys.stdout.write("7...")
        sys.stdout.flush()
        distance_indicator, distance_indicator_ref = galaxy.returnDistanceIndicator()
        print('Luminosity class')
        luminosity_class, luminosity_class_ref = galaxy.returnLuminosityClass()
        EBminusV = galaxy.returnEBminusV()

        sys.stdout.write("8...")
        sys.stdout.flush()
        radialVelocity = galaxy.returnRadialVelocity()
        if isNumber(radialVelocity):
            if float(radialVelocity) >= 0:
                vcorr = calculatevcorr(degreesJ2000RA, degreesJ2000Dec, radialVelocity)
                if vcorr > 0:
                    distvcorr = vcorr / hubble_constant
                else:
                    distvcorr = "x"
            else:
                vcorr = "x"
                distvcorr = "x"
        else:
            vcorr = "x"
            distvcorr = "x"

        # I am assuming that diameters are found in arcmin, and I'm converting to arcsec
        diameters = galaxy.returnDiameters()
        angMaj = diameters[0]
        if isNumber(angMaj):
            angMaj = float(angMaj) * 60
        angMin = diameters[1]
        if isNumber(angMin):
            angMin = float(angMin) * 60
            # switch them if the major diameter is larger than the minor
            if isNumber(angMaj):
                if angMaj < angMin:
                    angMaj, angMin = angMin, angMaj
        else:
            angMin = "x"

        if isNumber(angMaj) and isNumber(angMin) and isNumber(distvcorr):
            linDiameters = calculateLinearDiameters(angMaj, angMin, distvcorr)
        else:
            linDiameters = ("x", "x")

        if isNumber(angMaj) and isNumber(angMin):
            inclination = calculateInclination(angMaj, angMin)
        else:
            inclination = "x"

        # --- photometry returns a dictionary with keys:
        # --- photDict = {'b','u','g','r','i','z','j','h','k','all'}
        sys.stdout.write("9...")
        sys.stdout.flush()
        all_photometry = galaxy.returnPhotometry_old()

        # --- pick out the 'best' photometry values, keep all of it though
        selected_photometry = photometry.select_photometry(all_photometry)
        # selected_photometry = all_photometry['b'].to_numpy()

        # --- get names, remove weird URL coding from name, pick preferred name
        alternativeNames = galaxy.returnNames()
        oldName = unquote_plus(name).replace("\n", "").strip()
        preferredName = pickPreferredName(alternativeNames, oldName)

        strippedAlternativeNames = []
        for alternate in alternativeNames:
            stripped = alternate.replace(" ", "")
            strippedAlternativeNames.append(stripped)

        # --- create DataFrame and write it out
        # df = pd.DataFrame()

        # --- write it all out
        objectInfoList = [
            preferredName.replace(" ", ""),
            oldName.replace(" ", ""),
            redshift,
            (degreesJ2000RA, degreesJ2000Dec),
            J2000RA_Dec,
            galacticLong_Lat,
            rIndependentDistMean_sd_min_max,
            morphology,
            distance_indicator,
            luminosity_class,
            EBminusV,
            radialVelocity,
            vcorr,
            (angMaj, angMin),
            linDiameters,
            distvcorr,
            inclination,
            selected_photometry['b'],
            selected_photometry['u'],
            selected_photometry['g'],
            selected_photometry['r'],
            selected_photometry['i'],
            selected_photometry['z'],
            selected_photometry['j'],
            selected_photometry['h'],
            selected_photometry['k'],
            strippedAlternativeNames
            ]

        # write info to file
        row = dict((f, o) for f, o in zip(fieldnames, objectInfoList))
        if opts.verbose:
            print("\n", "row: ", row)

        writer.writerow(row)

        sys.stdout.write("10")
        sys.stdout.flush()
        sys.stdout.write("\n")

        print("Query time: ", query_time)
        print("Time for other stuff: ", time.time() - other_start)
        print("Total retrieval time: ", time.time() - target_start)
        counter -= 1
        if counter == 0:
            print("Reached max_retrieve = {0}. Exiting...".format(max_retrieve))
            break

    writerOutFile.close()
    if opts.verbose:
        print("Done.")
        print()
        print("Results can be found in: {0}".format(full_path))
        print()
        print("Elapsed Time: {0}".format(time.time() - total_start))
        print()

###############################################################################

if __name__ == "__main__":
    __package__ = "expected.package.name"

    # parse commandline
    commandlineOptions = parse_commandline()
    # do the work
    main(commandlineOptions)
#
