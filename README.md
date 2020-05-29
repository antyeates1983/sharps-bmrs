# sharps-bmrs
Python code for automated extraction of bipolar magnetic regions from the HMI/SHARPs database.

## Overview

This Python code generates a database of bipolar magnetic regions by automated analysis of the <a href="http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm">Spaceweather HMI Active Region Patch (SHARP) data</a>. The original data are pulled in automatically using <a href="https://sunpy.org">sunpy</a> and are provided by the <a href="http://hmi.stanford.edu">Helioseismic and Magnetic Imager</a> on Solar Dynamics Observatory.

Full details of the purpose and design of the code are given in the paper <a href="">A.R. Yeates, How good is the bipolar approximation of active regions for surface flux transport? (in preparation)</a>.

If you just want to use the BMR data without downloading and running the Python code, a version of `bmrharps_evol.txt` (see below) is maintained at the <a href="">Solar Dynamo Dataverse</a>.

## Dependencies

The code is tested with Python 3.7.3. The only nonstandard libraries required are astropy, sunpy and drms.

## Usage

Run the script `main.py`. You should first set the parameters in this script, especially the time ranges and the output path.

Various output files are produced at each stage of the script:
- a file `allsharps.txt` with a single entry for each detected SHARP, with parameters of the fitted BMR for all SHARPs deemed suitable for BMR fitting (i.e. with `good=1` - see the paper for details).
- .png images in the directory `outputpath` showing all of the SHARPs, and the fitted BMRs where appropriate (classified as good or bad depending on whether a BMR is fitted or not). 
- netcdf files containing the magnetic field of each good SHARP, for use in the next stage below (or to drive other numerical simulations).
- a file `repeatpairs.txt` listing pairs of BMRs identified as repeats.
- .png images showing the time evolution of the axial dipole moment predicted by the SFT model for both the original SHARP and the fitted BMR.
- a final database file `bmrsharps_evol.txt` listing all good SHARPs from `allsharps.txt` with columns for the predicted asymptotic dipole moment added (both for the SHARP and the BMR). Here is an example:
```
SHARPs from 2010-07-20 00:00:00 to 2010-08-20 23:00:00
-- Produced by anthony.yeates[at]durham.ac.uk --
7
Grid resolution: 180 x 360, smoothing_param = 4
Selection criteria: (i) sep > 1 deg,  (ii) |imbalance| <  0.5
Last two columns use 10-year 1D SFT simulation with eta=466.8 km^2/s, v0=0.0155 km/s, p=2.33, no decay term.
------------------------------------------------------
SHARP	NOAA	CM time		Latitude	Carr-Longitude	Unsgnd flux	Imbalance	Dipole		Bip-Separation	Bip-Tilt	Bip-Dipole	Pred-Dip-Real	Pred-Dip-Bip
92	11089	2010-07-25	-22.90078	203.23528	2.05991e+22	-0.02122	9.70499e-03	5.68659e+00	1.67930e+02	9.61091e-03	9.51234e-04	9.79128e-04
97	0	2010-07-24	13.34402	220.46158	3.58470e+20	 0.37365	-2.40001e-05	1.31483e+00	-6.98935e+00	-3.34932e-05	-2.36123e-05	-3.94304e-05
98	11090	2010-07-29	22.49662	149.96277	5.08242e+20	-0.32211	5.19955e-04	5.30587e+00	2.90281e+01	5.15122e-04	5.44939e-05	5.93080e-05
104	11092	2010-08-04	15.89691	76.61377	1.30975e+22	-0.14607	3.23345e-02	9.35438e+00	3.96531e+01	3.19762e-02	2.26765e-02	2.61885e-02
```

Be aware that the code will take some time to run if you are trying to cover long time periods. However, if it fails at any point during data download, then you should be able to run it again and it will continue from where it left off. You can force it to start from scratch by setting `restart=False`.



## Author

*A.R. Yeates - Durham University, UK*
