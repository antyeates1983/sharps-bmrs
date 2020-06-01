"""
    Python module with routines associated to sharps-bmrs.
    
    A.R. Yeates, Durham University, 1/6/20
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from astropy.io import fits
import sunpy.map
from sunpy.coordinates.sun import carrington_rotation_number
from sunpy.coordinates import get_sun_L0
from scipy.interpolate import interp2d
from scipy.io import netcdf
from scipy.ndimage.filters import gaussian_filter as gauss
import drms
import datetime
import pickle

#--------------------------------------------------------------------------------
def bmr2d(sc, pc, lat0, lon0, sep0, tilt0, B0, xithresh=9):
    """
        Return magnetic field for 2d tilted magnetic bipole with peak strength B0.
        **Note that all input angles should be in radians.
        
        In an untilted frame where the bipole is at the equator it has the form
        Br = -B0*(lon/sep0)*exp[ 0.5*( 1 - xi )]
        where xi = (lon**2 + 2*lat**2)/sep0**2
        where the strength B1 is chosen so that the maximum of Br is B0.
        
        Note: cutoff based on threshold of xi.
    """
    
    ns = np.size(sc)
    nph = np.size(pc)
    ds = sc[1] - sc[0]
    dp = pc[1] - pc[0]
    sc2, pc2 = np.meshgrid(sc, pc, indexing='ij')
    sinc2 = np.sqrt(1 - sc2**2)
    
    # Cartesian coordinates:
    x = np.cos(pc2)*sinc2
    y = np.sin(pc2)*sinc2
    z = sc2
    
    # Rotate to frame where BMR is on equator and untilted:
    xb = x*np.cos(lat0)*np.cos(lon0) + y*np.cos(lat0)*np.sin(lon0) + z*np.sin(lat0)
    yb = x*(-np.cos(tilt0)*np.sin(lon0) + np.sin(tilt0)*np.sin(lat0)*np.cos(lon0)) + y*(np.cos(tilt0)*np.cos(lon0) + np.sin(tilt0)*np.sin(lat0)*np.sin(lon0)) - z*np.sin(tilt0)*np.cos(lat0)
    zb = x*(-np.sin(tilt0)*np.sin(lon0) - np.cos(tilt0)*np.sin(lat0)*np.cos(lon0)) + y*(np.sin(tilt0)*np.cos(lon0)-np.cos(tilt0)*np.sin(lat0)*np.sin(lon0)) + z*np.cos(tilt0)*np.cos(lat0)
    zb[zb > 1] = 1
    zb[zb < -1] = -1
    
    # Magnetic field of BMR in this frame:
    thb = np.arccos(zb)
    phb = np.arctan2(yb, xb)
    xi = (phb**2 + 2*(0.5*np.pi-thb)**2)/sep0**2
    brb = -B0*phb/sep0*np.exp(-xi)
    
    # Cutoff at threshold:
    brb[xi > xithresh] = 0
    msk = (xi <= xithresh).astype('int')
    
    # Remove any flux imbalance:
    brb = correct_flux_multiplicative(brb)
    
    return brb, msk
    
#--------------------------------------------------------------------------------
def correct_flux_multiplicative(f):
    """
        Corrects the flux balance in the map f (assumes that cells have equal area).
    """
    
    # Compute positive and negative fluxes:
    ipos = f > 0
    ineg = f < 0
    fluxp = np.abs(np.sum(f[ipos]))
    fluxn = np.abs(np.sum(f[ineg]))
    
    # Rescale both polarities to mean:
    fluxmn = 0.5*(fluxn + fluxp)
    f1 = f.copy()
    f1[ineg] *= fluxmn/fluxn
    f1[ipos] *= fluxmn/fluxp
    
    return f1
    
#--------------------------------------------------------------------------------
def derotate(map, days):
    """
    Derotates 2D br map on computational grid according to differential rotation over period
    given by "days" argument.
    """
    # Coordinates:
    ns, nph = np.shape(map)
    ds = 2.0/ns
    dp = 2*np.pi/nph
    sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
    # Interpolator for original map (add ghost cells for periodicity):
    pcm = np.linspace(-0.5*dp, 2*np.pi + 0.5*dp, nph+2)
    mapm = np.zeros((ns, nph+2))
    mapm[:,1:-1] = map
    mapm[:,0] = mapm[:,-2]
    mapm[:,-1] = mapm[:,1]
    mapi = interp2d(pcm, sc, mapm, kind='linear', copy=True, bounds_error=False, fill_value=0)
    # Differential rotation profile:
    omA = np.deg2rad(0.18)
    omB = -np.deg2rad(2.396)
    omC = -np.deg2rad(1.787)
    om = omA + omB*sc**2 + omC*sc**4
    # Rotate back in time by "days":
    map1 = map*0
    for j in range(ns):
        pc1 = (pcm[1:-1] + om[j]*days.days + 2*np.pi) % (2*np.pi)
        for k in range(nph):
            map1[j,k] = mapi(pc1[k], sc[j]).flatten()
        # THE FOR LOOP SLOWS IT DOWN BUT HAVEN'T SUCCEEDED IN VECTORISING IT
    return map1
    
#--------------------------------------------------------------------------------
def evolve_forward(outputpath, sharpsfile, repeatfile, eta=466.8, v0=15.5e-3, p=2.33, outfile='bmrsharps_evol.txt', plots=True):
    """
    Reads list of regions from get_bmrs and list of repeats from get_pair_overlaps, selects only
    "good" regions that are not repeats, and evolves them forward in time according to the SFT model.
    """
    # GET LIST OF GOOD SHARPS:
    # - read full data:
    dat0 = np.loadtxt(outputpath+sharpsfile, skiprows=5, delimiter='\t', dtype={'names':('SHARP', 'NOAA', 'CM time', 'Latitude', 'Carr-Longitude', 'Unsgnd flux', 'Imbalance', 'Good', 'Dipole', 'Bip-Separation', 'Bip-Tilt', 'Bip-Dipole'), 'formats':('i4', 'i4', 'S10', 'f8', 'f8', 'f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'f8')})

    t0 = np.array([datetime.datetime.strptime(dat0['CM time'][i].decode('UTF-8'), '%Y-%m-%d') for i in range(len(dat0))])
    years0 = np.array([tt.year for tt in t0])
    # - remove bad regions:
    t1 = t0[dat0['Good']==1]
    years1 = years0[dat0['Good']==1]
    dat1 = dat0[dat0['Good']==1]
    # - remove repeats:
    repdat = np.loadtxt(outputpath+repeatfile, skiprows=2, delimiter='\t', dtype={'names':('reg1', 'reg2'), 'formats':('i4', 'i4')})
    uniq = np.full(np.size(t1), True)
    for k, sharpnum in enumerate(dat1['SHARP']):
        if np.isin(sharpnum, repdat['reg2']):   # repeated region
            uniq[k] = False
    t = t1[uniq]
    years = years1[uniq]
    dat = dat1[uniq]

    sharps = ['sharp%5.5i.nc' % k for k in dat['SHARP']]
    sharps.sort()
    print('Evolving %i regions forward by SFT.' % len(sharps))
     
    # Read in first region file to get resolution and set up flux transport object:
    fh = netcdf.netcdf_file(outputpath+sharps[0], 'r', mmap=False)
    br0 = fh.variables['br'][:]
    fh.close()
    ns, nph = np.shape(br0)

    sft = SFT1(ns)
    sftbip = SFT1(ns)

    # Initialize meridional flow, diffusion and decay:
    tau = 0
    def vs(s, v0=v0, p=p):
        Du = v0*(1+p)**(0.5*(p+1))/p**(0.5*p)
        return Du*s*(np.sqrt(1 - s**2))**p
    sft.prep_sft(vs, eta, tau)
    sftbip.prep_sft(vs, eta, tau)

    # Prepare new ASCII file:
    # - Load previous ASCII file:
    dat = np.loadtxt(outputpath+sharpsfile, skiprows=5, usecols=(0,1,2,3,4,5,6,7,8,9,10,11), delimiter='\t', dtype={'names':('SHARP', 'NOAA', 'CM time', 'Latitude', 'Carr-Longitude', 'Unsgnd flux', 'Imbalance', 'Good', 'Dipole', 'Bip-Separation', 'Bip-Tilt', 'Bip-Dipole'), 'formats':('i4', 'i4', 'S10', 'f8', 'f8', 'f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'f8')})
    # - Get header rows:
    with open(outputpath+sharpsfile) as myfile:
        heads = [next(myfile) for x in range(5)]
    # - Create new ASCII file incorporating evolution information:
    newfile = open(outputpath+outfile, 'w')
    newfile.write(heads[0])
    newfile.write('-- Produced by anthony.yeates[at]durham.ac.uk --\n')
    newfile.write('%i\n' % (len(sharps)))
    newfile.write(heads[1])
    newfile.write(heads[2])
    newfile.write('Last two columns use 10-year 1D SFT simulation with eta=%g km^2/s, v0=%g km/s, p=%g, no decay term.\n' % (eta, v0, p))
    newfile.write('------------------------------------------------------\n')
    heads[4] = 'SHARP\tNOAA\tCM time\t\tLatitude\tCarr-Longitude\tUnsgnd flux\tImbalance\tDipole\t\tBip-Separation\tBip-Tilt\tBip-Dipole\n'  # remove Good column
    newfile.write(heads[4][:-1]+'\tPred-Dip-Real\tPred-Dip-Bip\n')

    # Run simulations for each good SHARP:
    for sharp in sharps:
        fh = netcdf.netcdf_file(outputpath+sharp, 'r', mmap=False)
        try:
            brb = fh.variables['br_bipole'][:]
            print(sharp)
            good = True   # a bipole array was only created if it was a "good" SHARP
        except:
            good = False
        if (good):
            # Simulate evolution of bipole approximation:
            sftbip.setbr(np.mean(brb, axis=1))
                   
            Dmbip = [1.5*np.sum(sftbip.br*sftbip.sc*sftbip.ds)]
            for k in range(10*365):
                sftbip.evolve(sftbip.ndt*1)
                Dmbip.append(1.5*np.sum(sftbip.br*sftbip.sc*sftbip.ds))
                
            # Simulate evolution of original region:
            br0 = fh.variables['br'][:]
            sft.setbr(np.mean(br0, axis=1))
            Dm = [1.5*np.sum(sft.br*sft.sc*sft.ds)]
            t = [0]
            for k in range(10*365):
                sft.evolve(sft.ndt*1)
                Dm.append(1.5*np.sum(sft.br*sft.sc*sft.ds))
                t.append(k/365.)
                
            if (plots):
                plt.figure(figsize=(6,4))
                plt.plot(t, Dm, 'tab:blue', label='Original Br')
                plt.plot(t, Dmbip, 'tab:red', label='Bipole')
                dmax = max([np.max(np.abs(Dm)), np.max(np.abs(Dmbip))])
                plt.ylim(-dmax, dmax)
                plt.xlim(t[0], t[-1])
                plt.plot([t[0], t[-1]], [0,0], 'k-', linewidth=0.5)
                plt.xlabel('Years')
                plt.ylabel('Dipole')
                plt.title('SHARP '+sharp[5:10])
                plt.legend()
                plt.savefig(outputpath+'Devol_'+sharp[5:10]+'.png', bbox_inches='tight')
            
            Df = Dm[-1]
            Dfbip = Dmbip[-1]
        else:
            Df = 0
            Dfbip = 0
        
        k = np.where(dat['SHARP'] == int(sharp[5:10]))[0]
        newfile.write('%i\t%i\t%s\t%8.5f\t%8.5f\t%8.5e\t%8.5f\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\t%8.5e\n' % (dat['SHARP'][k], dat['NOAA'][k], dat['CM time'][k][0].decode('UTF-8'), dat['Latitude'][k], dat['Carr-Longitude'][k], dat['Unsgnd flux'][k], dat['Imbalance'][k], dat['Dipole'][k], dat['Bip-Separation'][k], dat['Bip-Tilt'][k], dat['Bip-Dipole'][k], Df, Dfbip))

        fh.close()
        
    newfile.close()

    
#--------------------------------------------------------------------------------
def get_bmrs(outputpath, outputfile, t_start, t_end, ns=180, nph=360, sharps_smoothing=4, imbalance_threshold=0.5, bmr_a=0.56, restart=True, plots=False):
    """
    Read magnetograms for list of SHARPs in pre-existing pickle file, and fit BMRs to those
    that are "good" (see paper). Saves output to text file.
    """
    
    print('Reading SHARP list and fitting BMRs...')

    picklefile = open(outputpath+'sharps.p','rb')
    sharps = pickle.load(picklefile)
    t_rec = pickle.load(picklefile)
    t_em = pickle.load(picklefile)
    picklefile.close()
    
    # Coordinate arrays etc:
    ds = 2.0/ns
    dph = 2*np.pi/nph
    sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
    pc = np.linspace(0.5*dph, 2*np.pi - 0.5*dph, nph)
    sc2, pc2 = np.meshgrid(sc, pc, indexing='ij')
    R = 6.96e10   # radius of Sun in cm

    # Initialise the file unless it already exists and restart=True:
    try:
        # Determine where to restart from by counting entries in file:
        with open(outputpath+outputfile) as f:
            count = 0
            for line in f:
                count += 1
        if (restart == True):
            allfile = open(outputpath+outputfile, 'a')
            print('Continuing existing file...')
            krestart = count - 5
    except:
        restart = False
    if (restart == False):
        allfile = open(outputpath+outputfile, 'w')
        allfile.write('SHARPs from %s to %s\n' % (t_start, t_end))
        allfile.write('Grid resolution: %i x %i, smoothing_param = %i\n' % (ns, nph, sharps_smoothing))
        allfile.write('Selection criteria: (i) sep > %g deg,  (ii) |imbalance| <  %g\n' % (np.rad2deg(dph), imbalance_threshold))
        allfile.write('\n')
        allfile.write('SHARP\tNOAA\tCM time\t\tLatitude\tCarr-Longitude\tUnsgnd flux\tImbalance\tGood\tDipole\t\tBip-Separation\tBip-Tilt\tBip-Dipole\n')
        krestart = 0

    if (plots):
        plt.figure(figsize=(10,6))
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        bmax = 50

    for k in range(krestart, len(sharps)):
        print(t_em[k], sharps[k])
        # - generate Br map for SHARP, and record imbalance (before correction):
        br1, pcen1, imbalance, keys = readSHARP(sharps[k], t_rec[k], ns, nph, sm=sharps_smoothing)
        flux1 = np.sum(np.abs(br1))
        # - get overall centroid position:
        ph1 = np.sum(np.abs(br1)*pc2)/np.sum(np.abs(br1)) + pcen1 - np.pi
        s1 = np.sum(np.abs(br1)*sc2)/np.sum(np.abs(br1))
        lat1 = np.arcsin(s1)
        # - compute axial dipole moment:
        ax_sharp = 0.75/np.pi*np.sum(br1*sc2)*ds*dph
        # - generate approximating BMR:
        phpos = np.sum(br1[br1 > 0]*pc2[br1 > 0])/np.sum(br1[br1 > 0])
        spos = np.sum(br1[br1 > 0]*sc2[br1 > 0])/np.sum(br1[br1 > 0])
        phneg = np.sum(br1[br1 < 0]*pc2[br1 < 0])/np.sum(br1[br1 < 0])
        sneg = np.sum(br1[br1 < 0]*sc2[br1 < 0])/np.sum(br1[br1 < 0])
        phcen = 0.5*(phpos + phneg)
        scen = 0.5*(spos + sneg)
        lat0 = np.arcsin(scen)
        lon0 = phcen - np.pi + pcen1
        sep0 = np.arccos(spos*sneg + np.sqrt(1-spos**2)*np.sqrt(1-sneg**2)*np.cos(phpos - phneg))
        tilt0 = np.arctan2(np.arcsin(spos)-np.arcsin(sneg), np.sqrt(1-scen**2)*(phneg - phpos))
        # - select "good" regions:
        if (~np.isnan(sep0) & (sep0 >= dph) & (np.abs(imbalance) < imbalance_threshold)):
            brb, msk = bmr2d(sc, pc, lat0, lon0, bmr_a*sep0, tilt0, 1, xithresh=9)
            brb *= flux1/np.sum(np.abs(brb))
            ax_bmr = 0.75/np.pi*np.sum(brb*sc2)*ds*dph
            good = True
        else:
            tilt0, ax_bmr = 0, 0
            good = False
        # - output data for this SHARP:
        allfile.write('%i\t%i\t%s\t%8.5f\t%8.5f\t%8.5e\t%8.5f\t%1.1i\t%8.5e\t%8.5e\t%8.5e\t%8.5e\n' % (sharps[k], keys.NOAA_AR, t_em[k].strftime('%Y-%m-%d'), np.rad2deg(lat1), np.rad2deg(ph1), flux1*ds*dph*R**2, imbalance, int(good), ax_sharp, np.rad2deg(sep0), np.rad2deg(tilt0), ax_bmr))
        # Rotate br to correct longitude:
        br1 = np.roll(br1, int((pcen1 - np.pi)/dph), axis=1)
        # Plot region:
        if (plots):
            plt.clf()
            ax = plt.subplot(211)
            lat = 0.5*np.pi - np.arccos(sc)
            pm = ax.pcolormesh(np.rad2deg(pc), np.rad2deg(lat), br1, cmap='bwr')
            pm.set_clim(vmin=-bmax, vmax=bmax)
            plt.title('Time %s -- SHARP %i, NOAA %i, FLUX = %g' % (t_em[k], sharps[k], keys.NOAA_AR, flux1*ds*dph*(6.96e10)**2))
            cb1 = plt.colorbar(pm)
            cb1.set_clim(vmin=-bmax, vmax=bmax)
            if (good):
                ax = plt.subplot(212)
                lat = 0.5*np.pi - np.arccos(sc)
                pm = ax.pcolormesh(np.rad2deg(pc), np.rad2deg(lat), brb, cmap='bwr')
                pm.set_clim(vmin=-bmax, vmax=bmax)
                cb1 = plt.colorbar(pm)
                cb1.set_clim(vmin=-bmax, vmax=bmax)
                ax.set_title('D(orig) = %g, D(BMR) = %g, sep = %1.2f, imb = %1.3f' % (ax_sharp, ax_bmr, np.rad2deg(sep0), imbalance))
            plt.tight_layout()
            if (good):
                plt.savefig(outputpath+'b_good_sharp%5.5i.png' % sharps[k], bbox_inches='tight')
            else:
                plt.savefig(outputpath+'b_bad_sharp%5.5i.png' % sharps[k], bbox_inches='tight')
            plt.clf()
            
        # Output to file:
        fid = netcdf.netcdf_file(outputpath+'sharp%5.5i.nc' % sharps[k], 'w')
        fid.createDimension('sdim', ns)
        fid.createDimension('pdim', nph)
        brid = fid.createVariable('br', 'd', ('sdim','pdim'))
        brid[:] = br1
        if (good):
            brbid = fid.createVariable('br_bipole', 'd', ('sdim','pdim'))
            brbid[:] = brb
        fid.close()
           
    allfile.close()
    
#--------------------------------------------------------------------------------
def get_pair_overlaps(outputpath, sharpsfile, repeat_threshold=1, outfile='repeatpairs.txt', plots=True, bmin=1e-12, t_restart=''):
    """
    Reads data file produced by get_bmrs, identifies pairs of repeat regions, and saves list of their numbers to text file.
    """
    
    print('Identifying repeat regions...')
    
    # Read "good" data (i.e. ones with BMRs):
    dat = np.loadtxt(outputpath+sharpsfile, skiprows=5, delimiter='\t', dtype={'names':('SHARP', 'NOAA', 'CM time', 'Latitude', 'Carr-Longitude', 'Unsgnd flux', 'Imbalance', 'Good', 'Dipole', 'Bip-Separation', 'Bip-Tilt', 'Bip-Dipole'), 'formats':('i4', 'i4', 'S10', 'f8', 'f8', 'f8', 'f8', 'i4', 'f8', 'f8', 'f8', 'f8')})
    dat = dat[:][np.where(dat['Good'] == 1)]
    times = [datetime.datetime.strptime(d.decode("utf-8") , '%Y-%m-%d') for d in dat['CM time']]

    # Prepare new ASCII file:

    if (plots):
        plt.figure(figsize=(12,4))
        gs = gridspec.GridSpec(2, 3, width_ratios=[1, 1, 1])
        bmax = 200
        # - Choose plotting boundaries:
        pmin = 0
        pmax = 360
        smin = -1
        smax = 1
        # - get coordinates:
        f = netcdf.netcdf_file(outputpath+'sharp%5.5i.nc' % dat['SHARP'][0], 'r', mmap=False)
        br = f.variables['br'][:]
        f.close()
        ns, nph = np.shape(br)
        del(br)
        ds = 2.0/ns
        dp = 2*np.pi/nph
        sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)
        pc = np.linspace(0.5*dp, 2*np.pi - 0.5*dp, nph)
    
    cnt = 0
    
    # If not restarting, set initial time and initialise ASCII file:
    if (t_restart==''):
        t_restart = np.min(times)
        repeatfile = open(outputpath+outfile, 'w')
        repeatfile.write('List of pairs of repeat SHARPs with repeat_threshold=%g\n' % repeat_threshold)
        repeatfile.write('-- Produced by anthony.yeates[at]durham.ac.uk --\n')
    else:
        repeatfile = open(outputpath+outfile, 'a')

    # Loop through each region:
    for k1, sharp1 in enumerate(dat['SHARP']):
        print(sharp1)
        # Load footprint of region on computational grid:
        f = netcdf.netcdf_file(outputpath+'sharp%5.5i.nc' % sharp1, 'r', mmap=False)
        br1 = f.variables['br'][:]
        br1b = f.variables['br_bipole'][:]
        f.close()
        reg1 = (np.abs(br1) > bmin).astype('float')
        reg1b = (np.abs(br1b) > bmin).astype('float')
        flux1 = np.sum(np.abs(br1))

        for k2 in range(k1+1, len(dat)):
            # Select regions later in list with CM passage from 20 to 34 days later:
            if ((times[k2] > t_restart) & (times[k2] - times[k1] >= datetime.timedelta(days=20)) & (times[k2] - times[k1] <= datetime.timedelta(days=34))):
                
                sharp2 = dat['SHARP'][k2]
            
                # - Load footprint of region on computational grid:
                f = netcdf.netcdf_file(outputpath+'sharp%5.5i.nc' % sharp2, 'r', mmap=False)
                br2 = f.variables['br'][:]
                br2b = f.variables['br_bipole'][:]
                f.close()
                reg2 = (np.abs(br2) > bmin).astype('float')
                flux2 = np.sum(np.abs(br2))

                if (flux2 > flux1):
                    continue

                # - Rotate back differentially:
                dt = (times[k2] - times[k1]).days
                br2d = derotate(br2, times[k2] - times[k1])
                reg2d = (np.abs(br2d) > bmin).astype('float')

                # - Flux of BMR 1 in derotated BMR 2:
                ofluxb1 = np.sum(np.abs(br1)*reg2d)

                olap = ofluxb1/flux2
                
                # - Get number of pixels in BMR representations
                npix1b = np.sum(reg1)
                npix2bd = np.sum(reg2d)
                
                # - If this is greater than overlap_threshold, record and plot:
                if ((olap > repeat_threshold) & (npix1b*npix2bd > 0)):
                
                    if (plots):
                        # Load corresponding HMI synoptic maps for context:
                        crot1 = carrington_rotation_number(times[k1])
                        brm1, scm1, pcm1 = readmap(crot1)
                        crot2 = carrington_rotation_number(times[k2])
                        brm2, scm2, pcm2 = readmap(crot2)
                        
                        # Compare two regions on grid:
                        ax = plt.subplot(gs[0])
                        pm = ax.pcolormesh(np.rad2deg(pcm1), scm1, brm1, cmap='bwr')
                        plt.contour(np.rad2deg(pc), sc, reg1, [0.5], colors='g', linewidths=0.75)
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(10, 0.8, '(a) CR%4.4i' % crot1)
                        ax.set_ylabel('Sine Latitude')

                        ax = plt.subplot(gs[1])
                        pm = ax.pcolormesh(np.rad2deg(pc), sc, br1, cmap='bwr')
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(pmin+1./36*(pmax-pmin), smax-0.1*(smax-smin), '(b) SHARP %i' % sharp1)
                        ax.set_xlim(pmin, pmax)
                        ax.set_ylim(smin, smax)

                        ax = plt.subplot(gs[2])
                        pm = ax.pcolormesh(np.rad2deg(pc), sc, br1b, cmap='bwr')
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(pmin+1./36*(pmax-pmin), smax-0.1*(smax-smin), r'(c) $\Phi$ = %6.2e' % (flux1*ds*dp*6.96e10**2))
                        # plt.contour(np.rad2deg(pc), sc, reg2bd, [0.5], colors='k', linewidths=0.75, linestyles='--')
                        ax.text(pmin+1./36*(pmax-pmin), smax-0.95*(smax-smin), r'$R_{12}$ = %4.2f' % olap)
                        ax.set_xlim(pmin, pmax)
                        ax.set_ylim(smin, smax)
                        
                        ax = plt.subplot(gs[3])
                        pm = ax.pcolormesh(np.rad2deg(pcm2), scm2, brm2, cmap='bwr')
                        plt.contour(np.rad2deg(pc), sc, reg2, [0.5], colors='k', linewidths=0.75)
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(10, 0.8, '(d) CR%4.4i' % crot2)
                        ax.set_xlabel('Carrington Longitude')
                        ax.set_ylabel('Sine Latitude')

                        ax = plt.subplot(gs[4])
                        pm = ax.pcolormesh(np.rad2deg(pc), sc, br2d, cmap='bwr')
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(pmin+1./36*(pmax-pmin), smax-0.1*(smax-smin), '(e) SHARP %i' % sharp2)
                        ax.set_xlabel('Carrington Longitude')
                        ax.set_xlim(pmin, pmax)
                        ax.set_ylim(smin, smax)
                        
                        ax = plt.subplot(gs[5])
                        br2bd = derotate(br2b, times[k2] - times[k1])
                        pm = ax.pcolormesh(np.rad2deg(pc), sc, br2bd, cmap='bwr')
                        pm.set_clim(vmin=-bmax, vmax=bmax)
                        ax.text(pmin+1./36*(pmax-pmin), smax-0.1*(smax-smin), r'(f) $\Phi$ = %6.2e' % (flux2*ds*dp*6.96e10**2))
                        ax.set_xlabel('Carrington Longitude')
                        ax.set_xlim(pmin, pmax)
                        ax.set_ylim(smin, smax)
                    
                        plt.savefig(outputpath+'repeat1_%5.5i.png' % cnt, bbox_inches='tight')
                    
                        plt.clf()
                    
                        cnt += 1
                    
                    # Output to list of repeat regions:
                    repeatfile.write('%i\t%i\n' % (sharp1, sharp2))
    repeatfile.close()
    
#--------------------------------------------------------------------------------
def get_sharps(outputpath, t_start, t_end, restart=True):
    """
    Generate pickle file listing SHARPs within a timeframe.
    """

    # Start time:
    cr_start = int(carrington_rotation_number(t_start))

    # Identify all SHARP regions present within simulation timeframe:
    # [cache in pickle file for later use]
    try:
        if (restart):
            picklefile = open(outputpath+'sharps.p','rb')
            picklefile.close()
            print('Reading pre-computed SHARPs from sharps.p')
    except:
        restart = False
    if (restart == False):
        sharps, t_rec, t_em = getSHARPs(t_start, t_end)
        picklefile1 = open(outputpath+'sharps.p','wb')
        pickle.dump(sharps, picklefile1)
        pickle.dump(t_rec, picklefile1)
        pickle.dump(t_em, picklefile1)
        picklefile1.close()
        
#--------------------------------------------------------------------------------
def getSHARPs(t_start, t_end):
    """
    Identify all SHARPs within given time range and identify timestamp of their closest
    frame to central meridian.
    Returns (1) sorted list of SHARP numbers, (2) timestamp of frame, and (3) corresponding
    emergence time (next noon) as a datetime object.
    """

    # Identify SHARP numbers in time range:
    t1 = t_start.strftime('%Y.%m.%d_%H:%M:%S')
    t2 = t_end.strftime('%Y.%m.%d_%H:%M:%S')
    print('Identifying SHARPs to be considered between '+t1+' and '+t2+'...')
    c = drms.Client()
    k = c.query('hmi.sharp_cea_720s[]['+t1+'-'+t2+'@1d]', key='HARPNUM, T_REC')
    try:
        sharps1 = list(set(k.HARPNUM))
        sharps1.sort()
    except:
        print('ERROR: no SHARPs found in time range.')
        sys.exit()

    # Find closest time to central meridian for each HARP and corresponding "emergence time":
    t_rec = []
    t_em = []
    sharps = []
    for h in sharps1:
        print(h)
        t_rec1, t_em1 = SHARPtime(h)
        if (t_em1 >= t_start):
            t_rec.append(t_rec1)
            t_em.append(t_em1)
            sharps.append(h)

    # Sort into order:
    idx = sorted(range(len(t_em)), key=t_em.__getitem__)
    sharps = [sharps[i] for i in idx]
    t_rec = [t_rec[i] for i in idx]
    t_em = [t_em[i] for i in idx]

    print('Number of SHARPs to consider: %i' % len(sharps))
    return sharps, t_rec, t_em
    
#--------------------------------------------------------------------------------
def readmap(rot):
    """
        Reads the synoptic map for Carrington rotation rot.
        [no flux correction etc, so just for display]
        
        ARGUMENTS:
            rot is the number of the required Carrington rotation (e.g. 2190)
    """
    
    # (1) READ IN DATA
    # ----------------
    # The seg='Mr_polfil' downloads the polar field corrected data --> without NaN values and errors due to projection effect
    # Read "Polar Field Correction for HMI Line-of-Sight Synoptic Data" by Xudong Sun, 2018 to know more
    # Link: https://arxiv.org/pdf/1801.04265.pdf
    try:
        c = drms.Client()
        seg = c.query(('hmi.synoptic_mr_polfil_720s[%4.4i]' % rot), seg='Mr_polfil')
    except:
        print('Error downloading HMI synoptic map -- required rotation: %4.4i' % rot)
    
    # Extract data array: (data is stored in 2nd slot with No. 1 and not in the PRIMARY (No. 0), thus [1].data)
    # typical file structure:
    # No.    Name      Ver    Type        Cards   Dimensions     Format
    # 0     PRIMARY    1   PrimaryHDU       6     ()
    # 1                1   CompImageHDU     13    (3600, 1440)   int32
    brm = (fits.open('http://jsoc.stanford.edu' + seg.Mr_polfil[0]))[1].data

    # Coordinates of original map:
    nsm = np.size(brm, axis=0)
    npm = np.size(brm, axis=1)
    dsm = 2.0/nsm
    dpm = 2*np.pi/npm
    scm = np.linspace(-1 + 0.5*dsm, 1 - 0.5*dsm, nsm)
    pcm = np.linspace(0.5*dpm, 2*np.pi - 0.5*dpm, npm)
                
    return brm, scm, pcm
        
#--------------------------------------------------------------------------------
def readSHARP(sharpnum, t_cm, ns, nph, sm=4, nocomp=False):
    """
    Read magnetogram for SHARP region at specified time, and map to computational grid. Return br array and imbalance fraction (0 = flux balanced, 1 = unipolar).
    
    Currently used B-LOS.
    
    Argument sm is the standard deviation of the Gaussian filter used to smooth the magnetogram
    before interpolation.
    
    If nocomp=True, just returns original HMI magnetogram and mask arrays.
    
    Outputs:
        br -- magnetogram with only this region (shifted to 180 degrees longitude)
        pcen -- amount of longitudinal shift
        imbalance -- imbalance in polarity
        k -- metadata
    """
    
    # Get time series of longitudes of this SHARP region (0 is central meridian):
    c = drms.Client()
    k, seg = c.query(('hmi.sharp_cea_720s[%i]' % sharpnum)+'['+t_cm+']', key='USFLUX, CRPIX1, CRVAL1, CRPIX2, CRVAL2, CDELT1, CDELT2, NOAA_AR, AREA', seg='BITMAP, MAGNETOGRAM')

    # Download bounding data for region:
    im = fits.open('http://jsoc.stanford.edu' + seg.BITMAP[0]) 
    mask = im[0].data
    
    # Download l-o-s magnetogram data:
    im = fits.open('http://jsoc.stanford.edu' + seg.MAGNETOGRAM[0]) 
    blos = im[1].data
       
    # Get heliographic (Carrington) coordinates of image:
    ny, nx = blos.shape
    xmin = (1 - k.CRPIX1)*k.CDELT1 + k.CRVAL1
    xmax = (nx - k.CRPIX1)*k.CDELT1 + k.CRVAL1
    ymin = (1 - k.CRPIX2)*k.CDELT2 + k.CRVAL2
    ymax = (ny - k.CRPIX2)*k.CDELT2 + k.CRVAL2    
    lon = np.linspace(xmin+0.5*k.CDELT1, xmax-0.5*k.CDELT1, nx)
    lat = np.linspace(ymin+0.5*k.CDELT2, ymax-0.5*k.CDELT2, ny)
    
    if (nocomp):
        return mask, blos, lon, lat
    
    # Remove data outside SHARP masked region:
    blos[mask < 30] = 0.0
    
    # Smooth with Gaussian filter:    
    blos = gauss(blos, sm)
    
    # Interpolate to global map in computational coordinates:
    pcm = np.deg2rad(lon)
    scm = np.sin(np.deg2rad(lat))
    
    # Shift to Carrington longitude 180 (to avoid problems at edge of map later):
    sc2, pc2 = np.meshgrid(scm, pcm, indexing='ij')
    pcen = np.sum(np.abs(blos)*pc2)/np.sum(np.abs(blos))
    scen = np.sum(np.abs(blos)*sc2)/np.sum(np.abs(blos))
    pcm += np.pi - pcen    
    
    # Interpolate to the computational grid:
    ds = 2.0/ns
    dph = 2*np.pi/nph
    sc = np.linspace(-1 + 0.5*ds, 1 - 0.5*ds, ns)  
    pc = np.linspace(0.5*dph, 2*np.pi - 0.5*dph, nph)
    
    bri = interp2d(pcm, scm, blos, kind='cubic', copy=True, bounds_error=False, fill_value=0)
    br = np.zeros((ns, nph))
    for i in range(ns):
        br[i,:] = bri(pc, sc[i]).flatten() + bri(pc + 2*np.pi, sc[i]).flatten() + bri(pc - 2*np.pi, sc[i]).flatten()
    del(bri)

    absflux = np.sum(np.abs(br))*ds*dph*(6.96e10)**2
    netflux = np.sum(br)*ds*dph*(6.96e10)**2
    imbalance = netflux/absflux
    
    # Correct flux balance:
    br = correct_flux_multiplicative(br)
    
    return br, pcen, imbalance, k

#--------------------------------------------------------------------------------
def SHARPtime(sharpnum):
    """
    For a given SHARP, return (1) timestamp of closest frame to central meridian, and (2) corresponding emergence time (next noon) as a datetime object.
    """
    
    # Get time series of longitudes of this SHARP (0 is central meridian):
    c = drms.Client()
    k = c.query('hmi.sharp_cea_720s[%i][]' % sharpnum, key='HARPNUM, T_REC, LON_FWT')

    # Find record closest to central meridian:
    rec_cm = k.LON_FWT.abs().idxmin()
    k_cm = k.loc[rec_cm]

    t_cm = drms.to_datetime(k.T_REC[rec_cm])

    # Identify emergence (completion) time - next noon:
    twelve_hrs = datetime.timedelta(hours=12)
    t_em = t_cm + twelve_hrs
    t_em = t_em.replace(hour=12, minute=0, second=0)

    return k.T_REC[rec_cm], t_em
    
#--------------------------------------------------------------------------------
class SFT1:
    """
        1D Surface flux transport model.
    """

    def __init__(self, ns):
        self.ns = ns
        self.ds = 2.0/ns
        self.sg = np.linspace(-1 + self.ds, 1 - self.ds, ns-1)
        self.sc = np.linspace(-1 + 0.5*self.ds, 1 - 0.5*self.ds, ns)  # excluding ghost cells

    def prep_sft(self, vs, eta, tau, L=6.96e5, cfl=0.2):
        """
            Prepare for surface flux transport.
            Input should be 1D arrays of vs and vp, and a constant eta.
        """
        self.eta = eta/L**2
        self.tau = tau*86400.0*365.25
        self.vs = vs(self.sg)/L
        
        # Timestep:
        dt_eta = self.ds**2/self.eta
        dt_mf = np.min(self.ds/np.abs(self.vs))
        dt = cfl*min([dt_eta, dt_mf])
        # - modify to fit exactly into one day:
        if (dt < 86400):
            self.ndt = int(86400.0/dt)
            self.dt = 86400.0/self.ndt
        else:
            self.ndt = 1
            self.dt = 86400.0
        print('Timestep dt = %g secs' % self.dt)
        
        # Include factor sin(theta) in vs:
        self.vs *= np.sqrt(1 - self.sg**2)
        
    def setbr(self, br0):
        """
            Set br.
        """
        self.br = br0.copy()

    def evolve(self, nsteps):
        """
            Evolve br by flux transport evolution for nsteps time steps.
            Uses simple finite-volume method with upwinding for advection term.
            Boundary conditions are zero flux at the poles.
        """

        for step in range(nsteps):
            # Evaluate FV fluxes [at interior ribs]:
            f = np.zeros(self.ns+1)
            # - diffusion:
            f[1:-1] = self.eta*(1 - self.sg**2)*(self.br[1:] - self.br[:-1])/self.ds
            # - meridional flow (by upwinding):
            f[1:-1] -= 0.5*(1 + np.sign(self.vs))*self.vs*self.br[:-1]
            f[1:-1] -= 0.5*(1 - np.sign(self.vs))*self.vs*self.br[1:]
            # Update including exponential decay:
            if (self.tau > 0):
                self.br += self.dt/self.ds*(f[1:] - f[:-1]) - self.dt/self.tau*self.br
            else:
                self.br += self.dt/self.ds*(f[1:] - f[:-1])
