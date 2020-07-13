#!/usr/bin/env python

'''

Parameters for NEDretriever.py

'''

#   fieldnames = ('Name',
#                 'NEDname',
#                 'redshift',
#                 'RAdeg',
#                 'DEdeg',
#                 'RA',
#                 'DEC',
#                 'GLON',
#                 'GLAT',
#                 'RID_mean',
#                 'RID_median',
#                 'RID_std',
#                 'RID_min',
#                 'RID_max',
#                 'MType',
#                 'distanceIndicator',
#                 'luminosityClass',
#                 'E(B-V)',
#                 'Vhel',
#                 'vcorr',
#                 'MajDiam_ang',
#                 'MinDiam_ang',
#                 'MajDiam',
#                 'MinDiam',
#                 'distvcorr',
#                 'inclination',
#                 'photometry',
#                 'altNames')


fieldnames = ('preferredName',\
'oldName',\
'redshift',\
'degreesJ2000RA_Dec',\
'J2000RA_Dec',\
'galacticLong_Lat',\
'rIndependentDistMean_sd_min_max (Mpc)',\
'morphology',\
'distanceIndicator',\
'luminosityClass',\
'EBminusV',\
'radialVelocity (km/s)',\
'vcorr (km/s)',\
'angDiameters (arcsec)',\
'linDiameters (kpc)',\
'distvcorr (Mpc)',\
'inclination (deg)',\
'photometry',\
'alternativeNames')

                
hubble_constant = 71.
max_retrieve = 10000.