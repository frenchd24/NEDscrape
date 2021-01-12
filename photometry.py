
'''
#!/usr/bin/env python

By David M. French (frenchd@astro.wisc.edu)

$Id: GTupdate2_photometry.py, v 2.3 11/06/17

Calculate best B band magnitudes and Lstar estimates from: return_phot_full.csv

combined with distance data from: return_basic_full.csv

Reminder: return_basic_full.csv, return_phot_full.csv, return_pa_full.csv come from:

return111111111111111.csv
returnPhot11111111111111.csv
returnPA11111111111111.csv

and are the result of taking out galaxies with vel or redshift = 'x' from the original
results.

Created processedPhot2.csv on (3/28/17), but have not included E(B-V) yet.

v2.1: fixed to include E(B-V) - made processedPhot2_extinc.csv (4/25/17)

reran to make processedPhot2_extinc3.csv - to insure everything is consistent
with GTupdate2_diameters.py (06/20/17)

v2.2: updated to use the new RID distances:
    reran AGAIN to make processedPhot3.csv - again because I've changed GTupdate2_diameters.py
    into GTupdate2_diameters_rid2.py (9/25/17)

    Input was:
    return_basic_full_extinc_rc3_rid2.csv
    return_phot_full.csv


v2.3: updated to use the new RID distances:
    reran AGAIN to make processedPhot4.csv - again because I updated the distances.
    pick_sdss_photometry() also had a few issues - 'mag' -> 'asinh mag' for cutting off
    the '+/-' from errors, and the order of cModel vs petro vs model was screwed up.
    Fixed now. Also, I fixed a small typo in defining the 'bestDistance' value below

    Input was:
    return_basic_full_extinc_rc3_rid5.csv
    return_phot_full.csv

    Output: processedPhot4.csv, rejected_processedPhot4.csv
    (11/06/17)

'''

import sys
# import os
# import csv
import warnings
# from pylab import *
import numpy as np
# import getpass
import math

import pandas as pd

# --- import program specific packages
from utilities import *
import params

###########################################################################
# --- define some helper functions

def calculate_absoluteMag_noExtinc(m, dm, d, dd):
    # m is apparent magnitude, d is distance in Mpc
    #
    # dm is the error in apparent magnitude, dd is the error in distance

    M = float(m) - 5*math.log10((float(d)*10**6)/10)

    # now do the error
    dM = math.sqrt(float(dm)**2 + ((5*float(dd))**2 / (float(d) * math.log10(10))**2))

    return M,dM


def calculate_absoluteMag(m, dm, d, dd, e):
    # m is apparent magnitude, d is distance in Mpc, e is extinction E(B-V)
    #
    # dm is the error in apparent magnitude, dd is the error in distance

    M = float(m) - 5*math.log10((float(d)*10**6)/10) - 3.1*float(e)

    # now do the error
    # error is dominated by dd, so don't even bother with an extinction error

    dM = math.sqrt(float(dm)**2 + ((5*float(dd))**2 / (float(d) * math.log10(10))**2))

    return M, dM


def calculate_lstar(mstar, m, dm):
    # calculate and return L/Lstar
    # mstar is the average M_star value, m is the absolute magnitude in question
    #
    # dm is the error in the absolute magnitude

    lratio = 10**(-0.4*(m-mstar))

    # now the error
    dLratio = 0.921034 * math.sqrt(10**(-0.8*(m-mstar)) * dm**2)

    return lratio,dLratio


def convert_sdss_to_B_jester(g, r):
    # convert SDSS photometry to B-band Johnson following Jester et al. 2005

    g = float(g)
    r = float(r)
    b = g + 0.39*(g - r) + 0.21
    return b


def convert_sdss_to_B_lupton(g, e_g, r, e_r):
    '''
    Convert SDSS photometry to B-band Johnson following Lupton 2005:
     - http://www.sdss.org/dr12/algorithms/sdssUBVRITransform/#Lupton2005

    sigma = 0.0107 for this conversion, but the measurement errors are assumed
    to dominate

    Inputs:
    ------
        g   :   float
            g-band photometry value
        e_g :   float
            error in g
        r   :   float
            r-band photometry value
        e_r :   float
            error in r

    Returns:
    -------
        b   :   float
            converted b-band value
        e_b :   float
            error in b

    '''

    g = float(g)
    r = float(r)
    e_g = float(e_g)
    e_r = float(e_r)

    # --- convert to B-Johnson
    b = g + 0.3130*(g - r) + 0.2271

    # --- calculate the error
    # --- the derivatives come out as:
    dbdr = -0.313
    dbdg = 1.313

    # --- now the error formula
    e_b = math.sqrt((dbdr * e_r)**2 + (dbdg * e_g)**2)

    return b, e_b



def pick_photometry_new(photList):
    '''
    Picks out the minimum, median, and maximum photometry values from an input
    list consisting of a single band.

    Inputs:
    ------
        photList    : array-like
            a nested list of photometry measurements. all the same band, in the
            format [[key, value, error, unit], [etc]]

    Returns:
    ------
        results     :  dictionary
            dictionary like {'median':{'value':val, 'error':error, 'key':key},
                            'min':{'value':val, 'error':error, 'key':key},
                            'max':{'value':val, 'error':error, 'key':key}}

    '''
    # output is either:
    # [min, median, and max] band values, along with the [keys for each]
    #
    # or:
    # output is [min_val, min_err, min_key], [med_val, med_err, med_key], [max_val, max_err, max_key]

    phot_df = pd.DataFrame({'passband'})

    vals = []
    errs = []
    keys = []
    for i in photList:
        key = i[0]
        val = i[1]
        err = i[2]
        unit = i[3]

        if unit == 'mag' and isNumber(val):
            # --- remove the '+/-' from the error value
            try:
                err = str(err).replace('+/-','')
                if isNumber(err):
                    err = float(err)
                else:
                    err = 'x'
            except Exception as e:
                sys.stdout.write("\n Error: {0}".format(e))
                err = 'x'

            # --- add them to the appropriate lists
            vals.append(val)
            errs.append(err)
            keys.append(key)

        elif unit == 'microJy':
            # --- convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val))

                    # --- then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err))
                        else:
                            mag_err = 'x'
                    except Exception as e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = 'x'

                    # --- add them to the appropriate lists
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)

        elif unit == 'Jy':
            # --- convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val)*10**6)

                    # --- then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err)*10**6)
                        else:
                            mag_err = float(mag_val)*0.001
                    except Exception as e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = float(mag_val)*0.001

                    # --- add them to the appropriate lists
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)
        else:
            sys.stdout.write('\n Error: Unknown unit = {}'.format(unit))

    if len(vals) >0:
        # --- now sort them, largest first
        all = zip(vals,errs,keys)
        all_sorted = sorted(all,reverse=True)

        # --- unzip them now that they are all sorted
        vals_sorted, errs_sorted, keys_sorted = zip(*all_sorted)

        # grab the median, min and max values
        # MIN MEANS BRIGHTER!
        val_med = median_low(vals_sorted)
        val_min = min(vals_sorted)
        val_max = max(vals_sorted)

        # --- find the index of the median within the list
        val_med_index = vals_sorted.index(val_med)

        # --- now use that index to grab the corresponding key and error
        err_med = errs_sorted[val_med_index]
        key_med = keys_sorted[val_med_index]

        # --- the min and max are easier:
        err_min = errs_sorted[-1]
        err_max = errs_sorted[0]

        key_min = keys_sorted[-1]
        key_max = keys_sorted[0]

        # results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min]

    else:
        # results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min] = ['x','x','x'],['x','x','x'],['x','x','x']

        val_med = 'x'
        err_med = 'x'
        key_med = 'x'

        val_min = 'x'
        err_min = 'x'
        key_min = 'x'

        val_max = 'x'
        err_max = 'x'
        key_max = 'x'

    # --- define and return results
    results = {'median':{'value':val_med, 'error':err_med, 'key':key_med},
                'min':{'value':val_min, 'error':err_min, 'key':key_min},
                'max':{'value':val_max, 'error':err_max, 'key':key_max}}

    return results





def pick_photometry(photList):
    '''
    Picks out the minimum, median, and maximum photometry values from an input
    list consisting of a single band.

    Inputs:
    ------
        photList    : array-like
            a nested list of photometry measurements. all the same band, in the
            format [[key, value, error, unit], [etc]]

    Returns:
    ------
        results     :  dictionary
            dictionary like {'median':{'value':val, 'error':error, 'key':key},
                            'min':{'value':val, 'error':error, 'key':key},
                            'max':{'value':val, 'error':error, 'key':key}}

    '''
    # output is either:
    # [min, median, and max] band values, along with the [keys for each]
    #
    # or:
    # output is [min_val, min_err, min_key], [med_val, med_err, med_key], [max_val, max_err, max_key]

    vals = []
    errs = []
    keys = []
    for i in photList:
        key = i[0]
        val = i[1]
        err = i[2]
        unit = i[3]

        if unit == 'mag' and isNumber(val):
            # --- remove the '+/-' from the error value
            try:
                err = str(err).replace('+/-','')
                if isNumber(err):
                    err = float(err)
                else:
                    err = 'x'
            except Exception as e:
                sys.stdout.write("\n Error: {0}".format(e))
                err = 'x'

            # --- add them to the appropriate lists
            vals.append(val)
            errs.append(err)
            keys.append(key)

        elif unit == 'microJy':
            # --- convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val))

                    # --- then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err))
                        else:
                            mag_err = 'x'
                    except Exception as e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = 'x'

                    # --- add them to the appropriate lists
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)

        elif unit == 'Jy':
            # --- convert to mags
            if isNumber(val):
                if float(val) >0:
                    mag_val = 23.9 - 2.5*math.log10(float(val)*10**6)

                    # --- then remove the '+/-' from the error value
                    try:
                        err = str(err).replace('+/-','')
                        if isNumber(err):
                            mag_err = 23.9 - 2.5*math.log10(float(err)*10**6)
                        else:
                            mag_err = float(mag_val)*0.001
                    except Exception as e:
                        sys.stdout.write("\n Error: {0}".format(e))
                        mag_err = float(mag_val)*0.001

                    # --- add them to the appropriate lists
                    vals.append(mag_val)
                    errs.append(mag_err)
                    keys.append(key)
        else:
            sys.stdout.write('\n Error: Unknown unit = {}'.format(unit))

    if len(vals) >0:
        # --- now sort them, largest first
        all = zip(vals,errs,keys)
        all_sorted = sorted(all,reverse=True)

        # --- unzip them now that they are all sorted
        vals_sorted, errs_sorted, keys_sorted = zip(*all_sorted)

        # grab the median, min and max values
        # MIN MEANS BRIGHTER!
        val_med = median_low(vals_sorted)
        val_min = min(vals_sorted)
        val_max = max(vals_sorted)

        # --- find the index of the median within the list
        val_med_index = vals_sorted.index(val_med)

        # --- now use that index to grab the corresponding key and error
        err_med = errs_sorted[val_med_index]
        key_med = keys_sorted[val_med_index]

        # --- the min and max are easier:
        err_min = errs_sorted[-1]
        err_max = errs_sorted[0]

        key_min = keys_sorted[-1]
        key_max = keys_sorted[0]

        # results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min]

    else:
        # results = [val_max,err_max,key_max],[val_med,err_med,key_med],[val_min,err_min,key_min] = ['x','x','x'],['x','x','x'],['x','x','x']

        val_med = 'x'
        err_med = 'x'
        key_med = 'x'

        val_min = 'x'
        err_min = 'x'
        key_min = 'x'

        val_max = 'x'
        err_max = 'x'
        key_max = 'x'

    # --- define and return results
    results = {'median':{'value':val_med, 'error':err_med, 'key':key_med},
                'min':{'value':val_min, 'error':err_min, 'key':key_min},
                'max':{'value':val_max, 'error':err_max, 'key':key_max}}

    return results


def pick_sdss_photometry(photList):
    '''
    Picks out best SDSS photometry value from an input list consisting
    of a single band. 1st choice: cmodel, 2nd: petrosian, 3rd: model, 4th: rest

    Inputs:
    ------
        photList    : array-like
            a nested list of photometry measurements. all the same band, in the
            format [[key, value, error, unit], [etc]]

    Returns:
    ------
        result      :  dictionary
            dictionary like {'value':val, 'error':error, 'key':key}

    '''

    vals = []
    errs = []
    keys = []

    d = {}
    for i in photList:
        key = i[0]
        val = i[1]
        err = i[2]
        unit = i[3]

        if unit == 'asinh mag':
            # remove the '+/-' from the error value
            try:
                err = str(err).replace('+/-','')
                if isNumber(err):
                    err = float(err)
                else:
                    err = float(val)*0.001
            except Exception as e:
                sys.stdout.write("\n Error: {0}".format(e))
                err = float(val)*0.001


        # pick out 'cModel' values
        if bfind(key.lower(),'cmodel'):
            if 'cmodel' in d:
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old, err_old, key_old] = d['cmodel']
                if val < val_old:
                    d['cmodel']=[val, err, key]
            else:
                d['cmodel']=[val, err, key]

        # pick out 'Petrosian' values
        elif bfind(key.lower(),'petrosian'):
            if 'petrosian' in d:
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old, err_old, key_old] = d['petrosian']
                if val < val_old:
                    d['petrosian']=[val, err, key]
            else:
                d['petrosian']=[val, err, key]

        # pick out 'model' values
        elif bfind(key.lower(),'model'):
            if 'model' in d:
                # check if this new measurement is brighter than the last
                # take the brightest (i.e., lowest)
                [val_old, err_old, key_old] = d['model']
                if val < val_old:
                    d['model']=[val, err, key]
            else:
                d['model']=[val, err, key]

        else:
            if not 'else' in d:
                d['else']=[val, err, key]

    # best choice
    if 'cmodel' in d:
        best = d['cmodel']

    # second best choice
    elif 'petrosian' in d:
        best = d['petrosian']

    # third best choice
    elif 'model' in d:
        best = d['model']

    # anything SDSS at all
    elif 'else' in d:
        best = d['else']

    # if something went wrong
    else:
        best = ['x', 'x', 'x']

    # --- turn it into a dictionary
    result = {'value':best[0], 'error':best[1], 'key':best[2]}

    # now return the best result
    return result


def derived_photometry(phot_dict):

    '''
    Derive absolute magnitude, Lstar, other

    NOT FINISHED

    '''

    # what velocity to assign to galaxies with negative or 0 redshifts?
    # initially 110 -> leads to a distance of ~1.6 Mpc, the radius of the local group
    # reset to 71.0 -> leads to a distance of 1.0 Mpc
    zeroVelocity = 71.

################################################################################
    # --- define the best distance
    # --- use the redshift independent distance if available

    if isNumber(RID2[1]):
        bestDistance = float(RID2[1])

        # use the given error, or assume 10% error
        if isNumber(RID2[2]):
            if RID2[2] == 0:
                distErr = bestDistance*0.1
            else:
                distErr = float(RID2[2])
        else:
            distErr = bestDistance*0.1
    else:
        # otherwise use Hubble's law to define a distance
        if float(vcorr) >0:
            bestDistance = float(vcorr)/hubbleConstant
            distErr = bestDistance*0.1
        else:
            bestDistance = float(zeroVelocity)/hubbleConstant
            distErr = bestDistance*0.5

################################################################################

    # --- calculate absolute magnitude
    if isNumber(b_mins[0]):
#                 M_b_min, M_b_min_err = calculate_absoluteMag_noExtinc(b_mins[0],b_mins[1],bestDistance,distErr)
#                 M_b_med, M_b_med_err = calculate_absoluteMag_noExtinc(b_meds[0],b_meds[1],bestDistance,distErr)
#                 M_b_max, M_b_max_err = calculate_absoluteMag_noExtinc(b_maxes[0],b_maxes[1],bestDistance,distErr)
        M_b_min, M_b_min_err = calculate_absoluteMag(b_mins[0],b_mins[1],bestDistance,distErr,extinc)
        M_b_med, M_b_med_err = calculate_absoluteMag(b_meds[0],b_meds[1],bestDistance,distErr,extinc)
        M_b_max, M_b_max_err = calculate_absoluteMag(b_maxes[0],b_maxes[1],bestDistance,distErr,extinc)
    else:
        M_b_min, M_b_min_err = 'x','x'
        M_b_med, M_b_med_err = 'x','x'
        M_b_max, M_b_max_err = 'x','x'

    # absolute magnitude from SDSS data
    if isNumber(sdss_b):
#                 M_sdss_b, M_sdss_b_err = calculate_absoluteMag_noExtinc(sdss_b,sdss_b_err,bestDistance,distErr)
        M_sdss_b, M_sdss_b_err = calculate_absoluteMag(sdss_b,sdss_b_err,bestDistance,distErr,extinc)
    else:
        M_sdss_b, M_sdss_b_err = 'x','x'

    # now calculate Lstar values
    if isNumber(M_b_min):
        Lstar_min, Lstar_min_err = calculate_lstar(averageBstar,M_b_min,M_b_min_err)
        Lstar_med, Lstar_med_err = calculate_lstar(averageBstar,M_b_med,M_b_med_err)
        Lstar_max, Lstar_max_err = calculate_lstar(averageBstar,M_b_max,M_b_max_err)
    else:
        Lstar_min, Lstar_min_err = 'x','x'
        Lstar_med, Lstar_med_err = 'x','x'
        Lstar_max, Lstar_max_err = 'x','x'

    if isNumber(M_sdss_b):
        Lstar_sdss, Lstar_sdss_err = calculate_lstar(averageBstar,M_sdss_b,M_sdss_b_err)
    else:
        Lstar_sdss, Lstar_sdss_err = 'x', 'x'


    if isNumber(Lstar_med):
        Lstar_med_val = round_to_sig(Lstar_med,sig=3)
        Lstar_med_err_val = round_to_sig(Lstar_med_err,sig=3)
    else:
        Lstar_med_val = 'x'
        Lstar_med_err_val = 'x'

    if isNumber(Lstar_max):
        Lstar_max_val = round_to_sig(Lstar_max,sig=3)
        Lstar_max_err_val = round_to_sig(Lstar_max_err,sig=3)
    else:
        Lstar_max_val = 'x'
        Lstar_max_err_val = 'x'

    if isNumber(Lstar_min):
        Lstar_min_val = round_to_sig(Lstar_min,sig=3)
        Lstar_min_err_val = round_to_sig(Lstar_min_err,sig=3)
    else:
        Lstar_min_val = 'x'
        Lstar_min_err_val = 'x'

    if isNumber(Lstar_sdss):
        Lstar_sdss_val = round_to_sig(Lstar_sdss,sig=3)
        Lstar_sdss_err_val = round_to_sig(Lstar_sdss_err,sig=3)
    else:
        Lstar_sdss_val = 'x'
        Lstar_sdss_err_val = 'x'


    pass



def select_photometry(phot_dict):

    '''
    Takes in a dictionary with keys: {'b','u','g','r','i','z','j','h','k'}
    Picks out the min, max, median measurements for b, j, h, k
    Picks out the best values u, g, r, i, z

    Inputs:
    ------
        phot_dict   :   dictionary
            keys are {'b','u','g','r','i','z','j','h','k'}, values are lists
            of all the corresponding measurements in that key's passband

    Returns:
    ------
        selected_phot   :   dictionary
            keys are {'b','u','g','r','i','z','j','h','k'}, values are
            nested dictionaries like:
                {'b': {'median':{'value':val, 'error':err, 'key':key},
                        'min': {'value':val, 'error':err, 'key':key},
                        'max': {'value':val, 'error':err, 'key':key}
                        },
                'u': {'value':val, 'error':err, 'key':key},
                'g': same as 'u',
                'r': same as 'u',
                'i': same as 'u',
                'z': same as 'u',
                'j': same as 'b',
                'h': same as 'b',
                'k': same as 'b'}

    '''

    # --- Mstar to use throughout
    # averageBstar = -19.57
    averageBstar = params.B_star

    # --- Hubble constant used throughout
    # hubbleConstant = 71.0
    hubbleConstant = params.hubble_constant

    # --- First deal with b, j, h, k photometry
    b_phot = pick_photometry(phot_dict['b'])
    j_phot = pick_photometry(phot_dict['j'])
    h_phot = pick_photometry(phot_dict['h'])
    k_phot = pick_photometry(phot_dict['k'])

    # --- best SDSS values now
    u_phot = pick_sdss_photometry(phot_dict['u'])
    g_phot = pick_sdss_photometry(phot_dict['g'])
    r_phot = pick_sdss_photometry(phot_dict['r'])
    i_phot = pick_sdss_photometry(phot_dict['i'])
    z_phot = pick_sdss_photometry(phot_dict['z'])

################################################################################
    # --- now convert to B from SDSS
    if isNumber(g_phot['value']) and isNumber(r_phot['value']):
        # --- if there are no errors, set them to 5% (originally 1%)
        if not isNumber(g_phot['error']):
            e_g = float(g_phot['value']) * 0.05
        else:
            e_g = float(g_phot['error'])
        if not isNumber(r_phot['error']):
            e_r = float(r_phot['value']) * 0.05
        else:
            e_r = float(r_phot['error'])

        # --- now do it - using Lupton's 2005 tranformation
        b_sdss, e_b_sdss = convert_sdss_to_B_lupton(g_phot['value'],
                                                    e_g,
                                                    r_phot['value'],
                                                    e_r)
    else:
        b_sdss, e_b_sdss = 'x', 'x'

    # --- add it to b_phot
    b_phot['sdss'] = {'value':b_sdss, 'error':e_b_sdss}

################################################################################
    selected_phot = {'b': b_phot,
                     'u': u_phot,
                     'g': g_phot,
                     'r': r_phot,
                     'i': i_phot,
                     'z': z_phot,
                     'j': j_phot,
                     'h': h_phot,
                     'k': k_phot}

    return selected_phot
