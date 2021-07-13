"""
    Run this script to extend an existing dataset up to the most recent available observation
    [or change t_newend to a different desired end date]

    Note that we use the **definitive** SHARP data series, so we never need to modify a SHARP that is already in our database.

    A.R Yeates, Durham University, updated Jul-2021
"""
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import sharptools
import pickle

# Database to be extended:
outputpath = '/Users/bmjg46/Documents/SHARPs3/'

# New end date:
#t_newend = datetime.datetime.combine(datetime.date.today(), datetime.time(23,59,59))
t_newend = datetime.datetime(2010, 9, 1, 00)

#--------------------------------------------------------
# Read parameters from allsharps.txt and repeatpairs.txt:
restart = True
bmr_a = 0.56  # this isn't written in the file
with open(outputpath+'allsharps.txt', 'r') as file1:
    head1 = file1.readline().split(' ')
    head2 = file1.readline().split(' ')
    head3 = file1.readline().split(' ')
ns, nph = int(head2[2]), int(head2[4].split(',')[0])
t_start = datetime.datetime.strptime(head1[2]+head1[3], '%Y-%m-%d%H:%M:%S')
t_end = datetime.datetime.strptime(head1[5]+head1[6], '%Y-%m-%d%H:%M:%S\n')
print(head2)
sharps_smoothing = int(head2[7].split(',')[0])
magtype = head2[10].split(',')[0]
method = head2[13].split(',')[0]
maxlon = int(head2[-1])
imbalance_threshold = float(head3[-1])
with open(outputpath+'repeatpairs.txt', 'r') as file1:
    head1 = file1.readline().split('=')
repeat_threshold = float(head1[-1])

# Load lists from existing pickle file:
picklefile = open(outputpath+'sharps.p','rb')
sharps = pickle.load(picklefile)
t_rec = pickle.load(picklefile)
t_em = pickle.load(picklefile)
picklefile.close()

# Identify new SHARPs and append to lists:
sharps1, t_rec1, t_em1 = sharptools.getSHARPs(t_end, t_newend, discardfile=outputpath+'limb_discards.txt')
for k,sharp1 in enumerate(sharps1):
    if ((sharp1 in sharps)==False):
        print(sharp1, t_rec1[k], t_em1[k])
        sharps.append(sharp1)
        t_rec.append(t_rec1[k])
        t_em.append(t_em1[k])

# Update pickle file:
picklefile1 = open(outputpath+'sharps.p','wb')
pickle.dump(sharps, picklefile1)
pickle.dump(t_rec, picklefile1)
pickle.dump(t_em, picklefile1)
picklefile1.close()
    
#===========================================================

# Read magnetograms for individual SHARPs and fit BMRs where possible:
sharptools.get_bmrs(outputpath, 'allsharps.txt', t_start, t_end, ns, nph, sharps_smoothing, imbalance_threshold, bmr_a, restart, plots=True, magtype=magtype, method=method, maxlon=maxlon)

# Determine repeated regions:
# [this considers also pairs where the first is in the existing data]
sharptools.get_pair_overlaps(outputpath, 'allsharps.txt', repeat_threshold=repeat_threshold, outfile='repeatpairs.txt', plots=True, t_restart=t_end)

# Evolve SHARPs with BMRs forward in time and generate final database file:
# [at present, this regenerates all of the existing data too]
sharptools.evolve_forward(outputpath, 'allsharps.txt', 'repeatpairs.txt', eta=350, v0=15e-3, p=2.33, outfile='bmrsharps_evol.txt', plots=True)
