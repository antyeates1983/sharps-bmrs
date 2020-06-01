# sharps-bmrs
Python code for automated extraction of bipolar magnetic regions from the HMI/SHARPs database.

## Overview

This Python code generates a database of bipolar magnetic regions by automated analysis of the <a href="http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htm">Spaceweather HMI Active Region Patch (SHARP) data</a>. The original data are pulled in automatically using <a href="https://sunpy.org">sunpy</a> and are provided by the <a href="http://hmi.stanford.edu">Helioseismic and Magnetic Imager</a> on Solar Dynamics Observatory.

Full details of the purpose and design of the code are given in the paper <a href="">A.R. Yeates, How good is the bipolar approximation of active regions for surface flux transport? (in preparation)</a>.

If you just want to use the BMR data without downloading and running the Python code, a version of `bmrharps_evol.txt` (see below) is maintained at the <a href="">Solar Dynamo Dataverse</a>.

## Credit

If you use the code or the derived data in your published work, please cite the accompanying paper [under review].

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
SHARPs from 2010-07-01 00:00:00 to 2010-09-01 00:00:00
-- Produced by anthony.yeates[at]durham.ac.uk --
11
Grid resolution: 180 x 360, smoothing_param = 4
Selection criteria: (i) sep > 1 deg,  (ii) |imbalance| <  0.5
Last two columns use 10-year 1D SFT simulation with eta=350 km^2/s, v0=0.015 km/s, p=2.33, no decay term.
------------------------------------------------------
SHARP	NOAA	CM time		Latitude	Carr-Longitude	Unsgnd flux	Imbalance	Dipole		Bip-Separation	Bip-Tilt	Bip-Dipole	Pred-Dip-Real	Pred-Dip-Bip
86	11087	2010-07-15	19.49569	335.55594	1.83895e+22	-0.03850	1.92913e-02	8.97116e+00	1.67619e+01	1.90763e-02	4.51918e-03	3.39163e-03
87	0	2010-07-09	19.52667	100.34961	3.89123e+20	 0.31597	8.57507e-05	1.09127e+00	2.98810e+01	1.07181e-04	1.02349e-05	1.26450e-05
89	11088	2010-07-15	-20.59370	337.44025	7.32340e+20	-0.40323	-4.75846e-05	2.18557e+00	-1.75765e+02	-5.06708e-05	-3.01310e-06	-4.10831e-06
92	11089	2010-07-25	-22.90078	203.23528	2.05991e+22	-0.02122	9.70499e-03	5.68659e+00	1.67930e+02	9.61091e-03	3.75159e-04	3.86887e-04
97	0	2010-07-24	13.34402	220.46158	3.58470e+20	 0.37365	-2.40001e-05	1.31483e+00	-6.98935e+00	-3.34932e-05	-1.63547e-05	-2.89587e-05
98	11090	2010-07-29	22.49662	149.96277	5.08242e+20	-0.32211	5.19955e-04	5.30587e+00	2.90281e+01	5.15122e-04	2.10756e-05	2.42236e-05
104	11092	2010-08-04	15.89691	76.61377	1.30975e+22	-0.14607	3.23345e-02	9.35438e+00	3.96531e+01	3.19762e-02	1.51923e-02	1.88848e-02
114	11095	2010-08-09	-16.74810	 4.55401	4.01993e+21	 0.22351	3.77086e-03	7.36436e+00	1.61984e+02	3.73241e-03	1.23802e-03	1.42898e-03
115	11093	2010-08-10	14.64567	348.66903	1.22733e+22	-0.49672	-6.87333e-03	3.36876e+00	-2.35482e+01	-6.81773e-03	-3.63902e-03	-4.21641e-03
131	11098	2010-08-14	14.04513	302.18062	2.14102e+21	 0.12191	4.55119e-04	3.55713e+00	8.23495e+00	4.47399e-04	3.66343e-04	3.29530e-04
146	11102	2010-08-29	26.81754	103.17549	2.51881e+21	-0.33635	5.07201e-04	3.23274e+00	9.33950e+00	4.99081e-04	5.70436e-06	3.65075e-06
```

Be aware that the code will take some time to run if you are trying to cover long time periods. However, if it fails at any point during data download, then you should be able to run it again and it will continue from where it left off. You can force it to start from scratch by setting `restart=False`.

By default the code is set up to run for a two-month period, and the corresponding `allsharps.txt`, `repeatpairs.txt` and `bmrsharps_evol.txt` files are provided here for verification.

## Author

*A.R. Yeates - Durham University, UK*
