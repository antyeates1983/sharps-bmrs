"""
    Main script to generate a list of BMRs (bipolar magnetic regions) from HMI/SHARPs data.
    
    Interpolates original map to equally-spaced grid in sin(latitude) and longitude for computation of the BMR properties. The resolution and degree of smoothing in this interpolation can be modified.
    
    Output:
    outputpath/allsharps.txt - list of all SHARPs whether selected or not
    
    A.R Yeates, Durham University, 29/5/20
"""
import os
import datetime
from datetime import datetime as dt_obj
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from sunpy.coordinates.sun import carrington_rotation_number
from scipy.io import netcdf
import pickle
import sharptools

#===========================================================
# SET INPUT PARAMETERS:

restart = True  # set to False if you want to regenerate all output

# Grid sizes in s=cos(theta) and phi:
ns = 180
nph = 360

# Required time range:
t_start = dt_obj(2010, 7, 20, 00)
t_end = dt_obj(2010, 8, 20, 23)

starttime = t_start.strftime('%Y%m%d.%H')
endtime = t_end.strftime('%Y%m%d.%H')

# Path to output directory [with trailing slash]:
outputpath = './'

# Parameters for emerging regions:
sharps_smoothing = 4
imbalance_threshold = 0.5
bmr_a = 0.56

# Parameters for removing repeat regions:
repeat_threshold = 1

#===========================================================

# Create output directory:
os.system('mkdir -p '+outputpath.replace(' ', '\ '))

# Query database and get pickle file outputpat/sharps.ps listing all SHARPs in interval:
sharptools.get_sharps(outputpath, t_start, t_end, restart)

# Read magnetograms for individual SHARPs and fit BMRs where possible:
sharptools.get_bmrs(outputpath, 'allsharps.txt', t_start, t_end, ns, nph, sharps_smoothing, imbalance_threshold, bmr_a, restart, plots=True)

# Determine repeated regions:
sharptools.get_pair_overlaps(outputpath, 'allsharps.txt', repeat_threshold=1, outfile='repeatpairs.txt', plots=True)

# Evolve SHARPs with BMRs forward in time and generate final database file:
sharptools.evolve_forward(outputpath, 'allsharps.txt', 'repeatpairs.txt', eta=466.8, v0=15.5e-3, p=2.33, outfile='bmrsharps_evol.txt', plots=True)
