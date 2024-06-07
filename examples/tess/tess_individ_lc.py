## To run script after installing EXOTIC locally from the Main branch do the following in terminal
# cd EXOTIC
# cd examples
# python tess_individ_lc.py -t "WASP-77 A b"
#To run a different target simply call the code again with your new target:
# python tess_individ_lc.py -t "K2-132 b"
import os
import copy
import json
import pickle
import argparse
import requests
import numpy as np
from io import StringIO
import pandas as pd
from pandas import read_csv
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation
from scipy.signal import savgol_filter, medfilt
from astropy import constants as const
from astropy import units as u
from astropy.time import Time
from scipy.ndimage import label
import statsmodels.api as sm
from pylightcurve import exotethys
import lightkurve as lk
from wotan import flatten
from exotic.api.elca import transit, lc_fitter
from exotic.api.output_aavso import OutputFiles
from transitleastsquares import transitleastsquares

def tap_query(base_url, query, dataframe=True):
    # table access protocol query

    # build url
    uri_full = base_url
    for k in query:
        if k != "format":
            uri_full += f"{k} {query[k]} "

    uri_full = f"{uri_full[:-1]} &format={query.get('format', 'csv')}"
    uri_full = uri_full.replace(' ', '+')
    print(uri_full)

    response = requests.get(uri_full, timeout=300)
    # TODO check status_code?

    if dataframe:
        return read_csv(StringIO(response.text))
    else:
        return response.text

def nea_scrape(target=None):

    # scrape_new()
    uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
    uri_ipac_query = {
        # Table columns: https://exoplanetarchive.ipac.caltech.edu/docs/API_PS_columns.html
        "select": "pl_name,hostname,pl_radj,pl_radjerr1,ra,dec,"#pl_massj,"
                    "pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,"
                    "pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbeccen,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,"
                    "pl_orblper,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,"
                    "pl_ratror,pl_ratrorerr1,pl_ratrorerr2,"
                    "st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,"
                    "st_logg,st_loggerr1,st_loggerr2,st_mass,st_rad,st_raderr1",
        "from": "pscomppars",  # Table name
        "where": "tran_flag = 1",
       #"order by": "pl_pubdate desc",s
        "format": "csv"
    }

    if target:
        uri_ipac_query["where"] += f" and pl_name = '{target}'"
    return tap_query(uri_ipac_base, uri_ipac_query)

def sigma_clip(ogdata,dt, iterations=1):

    mask = np.ones(ogdata.shape,dtype=bool)
    for i in range(iterations):
        mdata = savgol_filter(ogdata[mask], dt, 2)
        res = ogdata[mask] - mdata
        try:
            std = np.nanmedian([np.nanstd(np.random.choice(res,100)) for i in range(250)])
        except:
            std = np.nanstd(res)
        mask[mask] = np.abs(res) < 3*std

    mdata = savgol_filter(ogdata[mask], dt, 2)
    res = ogdata[mask] - mdata

    data = copy.deepcopy(ogdata)
    data[~mask] = np.nan
    return data, np.std(res)

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Choose a target to process"
    parser.add_argument("-t", "--target", help=help_, type=str, default="WASP-18 b")

    help_ = "Directory for saving results"
    parser.add_argument("-o", "--output", help=help_, default="tess_individ_lc_output/", type=str)

    help_ = "Sectors (0 = all)"
    parser.add_argument("-s", "--sector", help=help_, default=0, type=int)

    help_ = "Choose a target to re-process. Need to do this when after changing the prior."
    parser.add_argument("-r", "--reprocess", help=help_, action='store_true', default=False)

    help_ = "Add a/Rs to the fitting process"
    parser.add_argument("--ars", help=help_, action='store_true', default=False)

    help_ = "Perform transit least squares search over residuals of global fit"
    parser.add_argument("--tls", help=help_, action='store_true', default=False)

    return parser.parse_args()

def check_std(time, flux, dt=0.5): # dt = [hr]
    tdt = np.diff(np.sort(time)).mean()
    si = np.argsort(time)
    sflux = savgol_filter(flux[si], 1+2*int(max(15,dt/24/tdt)), 2)
    return np.nanstd(flux - sflux)

# compute stellar mass from logg and radius
stellar_mass = lambda logg,rs: ((rs*u.R_sun)**2 * 10**logg*(u.cm/u.s**2) / const.G).to(u.M_sun).value

# keplerian semi-major axis (au)
sa = lambda m,P : ((const.G*m*u.M_sun*P*u.day**2/(4*np.pi**2) )**(1./3)).to(u.AU).value

if __name__ == "__main__":
    infoDict = {}
    args = parse_args()
    
    if not os.path.exists(args.output):
        os.mkdir(args.output)

    planetdir = os.path.join(args.output, args.target.replace(' ','_').replace('-','_') )
    if not os.path.exists(os.path.join(planetdir,"global_fit.png")):
        if not os.path.exists(planetdir):
            os.mkdir(planetdir)
    else:
        if not args.reprocess:
            raise(Exception("Target already processed. Skipping..."))
    
    planetname = args.target.lower().replace(' ','')
    
    search_result = lk.search_targetpixelfile(args.target[:-1], mission='TESS')
    print(search_result)

    if len(search_result) == 0:
        raise(Exception(f"no data for: {args.target}"))
        # https://exo.mast.stsci.edu/
    # load prior from disk or download
    if os.path.exists(os.path.join(planetdir,planetname+"_prior.json")):
        prior = json.load(open(os.path.join(planetdir,planetname+"_prior.json"),"r"))
    else:
        # download prior from web
        if "TOI" in args.target:
            tdf = read_csv("toi_table.csv")
            try:
                if '-' in args.target:
                    num = float(args.target.split("-")[1])
                else:
                    num = float(args.target.split(" ")[1])
                    
            except:
                num = float(args.target.split("-")[1][:-2])
                planet = args.target.split("-")[1][-1]
                if planet == 'b':
                    num += 0.01
                elif planet == 'c':
                    num += 0.02
                elif planet == 'd':
                    num += 0.03
                elif planet == 'e':
                    num += 0.04

            ti = np.argwhere(tdf.TOI.values == num)[0][0]
            row = tdf.iloc[ti]

            prior = {
                "pl_name": num,
                "hostname": int(row['TIC ID']),
                "pl_radj": float(row['Planet Radius (R_Earth)']*u.R_earth.to(u.R_jupiter)),
                "pl_radjerr1": float(row['Planet Radius (R_Earth) err']*u.R_earth.to(u.R_jupiter)),
                "ra": row['RA'],
                "dec": row['Dec'],
                #"pl_massj"
                # "pl_ratdor": -1, # TODO kepler
                # "pl_ratdorerr1": +1,
                # "pl_ratdorerr2": -1,
                "pl_orbincl": 90,
                "pl_orbinclerr1": 0,
                "pl_orbinclerr2": -10,
                "pl_orbper": float(row['Period (days)']),
                "pl_orbpererr1": float(row['Period (days) err']),
                "pl_orbpererr2": float(-row['Period (days) err']),
                "pl_orbeccen": 0.0,
                # "pl_orbsmax": -1,
                # "pl_orbsmaxerr1": 0.00022,
                # "pl_orbsmaxerr2": -0.00022999999999999998,
                "pl_orblper": 0,
                "pl_tranmid": float(row['Epoch (BJD)']),
                "pl_tranmiderr1": float(row['Epoch (BJD) err']),
                "pl_tranmiderr2": float(-row['Epoch (BJD) err']),
                #"pl_ratror": -1,
                # "pl_ratrorerr1": 0,
                # "pl_ratrorerr2": 0,
                "st_teff": float(row['Stellar Eff Temp (K)']),
                "st_tefferr1": float(row['Stellar Eff Temp (K) err']),
                "st_tefferr2": float(-row['Stellar Eff Temp (K) err']),
                "st_met": float(row['Stellar Metallicity']),
                #"st_meterr1": float(row[' Stellar Metallicity err']),
                "st_meterr1": 0,
                #"st_meterr2": float(-row[' Stellar Metallicity err']),
                "st_meterr2": 0,
                "st_logg": float(row['Stellar log(g) (cm/s^2)']),
                "st_loggerr1": float(row['Stellar log(g) (cm/s^2) err']),
                "st_loggerr2": float(-row['Stellar log(g) (cm/s^2) err']),
                "st_mass": float(stellar_mass( row['Stellar log(g) (cm/s^2)'],row['Stellar Radius (R_Sun)'])),
                "st_rad": float(row['Stellar Radius (R_Sun)']),
                "st_raderr1": float(row['Stellar Radius (R_Sun) err'])
            }

            prior['pl_orbsmax'] = sa(prior['st_mass'], row['Period (days)'])
            prior['pl_orbsmaxerr1'] = prior['pl_orbsmax']*0.1
            prior['pl_orbsmaxerr2'] = -prior['pl_orbsmax']*0.1
            
            prior['pl_ratdor'] = ((prior['pl_orbsmax']*u.AU)/(prior['st_rad']*u.R_sun)).to(u.AU/u.AU).value
            prior['pl_ratdorerr1'] = prior['pl_ratdor']*0.1
            prior['pl_ratdorerr2'] = -prior['pl_ratdor']*0.1

            prior['pl_ratror'] = float((prior['pl_radj']*const.R_jup / (prior['st_rad']*const.R_sun)))
            prior['pl_ratrorerr1'] = prior['pl_ratror']*0.1
            prior['pl_ratrorerr2'] = -prior['pl_ratror']*0.1

            if np.isnan(prior['st_met']):
                prior['st_met'] = 0
                prior['st_meterr1'] = 0.1
                prior['st_meterr2'] = -0.1

            prior['st_loggerr1'] = 0.01*prior['st_logg']
            prior['st_loggerr2'] = -0.01*prior['st_logg']

            prior['st_tefferr1'] = 0.01*prior['st_teff']
            prior['st_tefferr2'] = -0.01*prior['st_teff']

            for k in prior:
                if isinstance(prior[k],str):
                    pass
                else:
                    if np.isnan(prior[k]):
                        if "err1" in k:
                            prior[k]=1
                        elif "err2" in k:
                            prior[k]=-1
                        else:
                            prior[k]=None

        else:
            nea_df = nea_scrape(args.target)
            prior = {} # convert to dict
            for key in nea_df.keys():
                prior[key] = nea_df[key][0]
        prior["pl_name"] = planetname.replace('-','')
        if np.isnan(prior['pl_ratdor']):
            prior['pl_ratdor'] = np.round(((prior['pl_orbsmax']*u.AU)/(prior['st_rad']*u.R_sun)).to(u.AU/u.AU).value,2)
            prior['pl_ratdorerr1'] = prior['pl_ratdor']*0.1
            prior['pl_ratdorerr2'] = -prior['pl_ratdor']*0.1
        # make directories
        if not os.path.exists(args.output):
            os.mkdir(args.output)

        #if not os.path.exists(os.path.join(planetdir,"lightcurves")):
            #os.mkdir(os.path.join(planetdir,"lightcurves"))

        # save prior to disk
        with open(os.path.join(planetdir,planetname+"_prior.json"), 'w', encoding ='utf8') as json_file:
            json.dump(prior, json_file, indent=4)

    if len(prior) == 0:
        raise(Exception(f"no target named: {args.target}"))
        # https://exoplanetarchive.ipac.caltech.edu/

    # find sectors
    sectors = [int(sector.split(' ')[-1]) for sector in search_result.mission]
    if args.sector:
        if args.sector not in sectors:
            #os.rmdir(args.output)
            raise(Exception(f"no data for sector: {args.sector}"))
    
    #if np.isnan(prior['pl_orbpererr1']):
            #prior['pl_orbpererr1'] = 0.
            #prior['pl_orbpererr2'] = 0. 

    # compute limb darkening
    if np.isnan(prior['st_meterr1']):
        prior['st_meterr1'] = prior['st_met']*0.1
        prior['st_meterr2'] = -prior['st_met']*0.1

    u0,u1,u2,u3 = exotethys(prior['st_logg'], prior['st_teff'], prior['st_met'], 'TESS', method='claret', stellar_model='phoenix')

    # alloc state vector
    sv = {
        'lightcurves':[], # dicts of time,flux,fluxerr,pars,errors from light curve fit
        'sectors':{},      # exp time of each sector that gets processed
    
        # global timeseries
        'time':[],
        'flux':[],
        'flux_err':[],
        'trend':[],
        'sector':[]
    }

    #  loop over data sets
    for i, sector in enumerate(sectors):
        
        # just process single sector
        if args.sector:
            if sector != args.sector:
                continue

        # don't reprocess same sector at lower res
        if search_result.exptime[i].value > sv['sectors'].get(sector,1800):
            print(f" Skipping Sector {sector} @ {search_result.exptime[i].value} s")
            print(f"  Already processed @ {sv['sectors'][sector]} s")
            continue
        if i == 0:
            exptime = search_result.exptime[i].value
        if exptime < search_result.exptime[i].value:
            exptime = search_result.exptime[i].value

        # TODO cmd line arg for skipping cadences >= t
        #if search_result.exptime[i].value == 1800: # 30 minute cadence
        #    print("Skipping 30 minute cadence data")
        #    continue

        print(f"Downloading Sector {sector} @ {search_result.exptime[i].value} s ...")
        try:
            tpf = search_result[i].download(quality_bitmask='hard')
        except Exception as err:
            print(f" Failed to download {search_result.exptime[i].value} s data for sector {sector}")
            print(f" {err}")
            continue

        #  store exp time
        sv['sectors'][sector] = search_result.exptime[i].value #EXPOSURE TIME!!!    
        infoDict['exposure'] = search_result.exptime[i].value
        # aper selection
        lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
        nmask = np.isnan(lc.flux.value)
        lstd = check_std(lc.time.value[~nmask], lc.flux.value[~nmask])

        # test two more apertures, slightly bigger each iteration
        for it in [1,2]:
            lcd = tpf.to_lightcurve(aperture_mask=binary_dilation(tpf.pipeline_mask, iterations=it))
            nmaskd = np.isnan(lcd.flux.value)
            lstdd = check_std(lcd.time.value[~nmask], lcd.flux.value[~nmask])

            if lstd < lstdd:
                aper_final = tpf.pipeline_mask
            else:
                aper_final = binary_dilation(tpf.pipeline_mask, iterations=it)

        # final aperture
        lc = tpf.to_lightcurve(aperture_mask=aper_final)

        # aperture plot
        tpf.plot(aperture_mask=aper_final)
        plt.savefig( os.path.join(planetdir, planetname+f"_sector_{sector}_aperture.png") )
        plt.close()

        # remove first ~30 min of data after any big gaps
        tmask = np.ones(lc.time.shape).astype(bool)
        smask = np.argsort(lc.time.value)
        dts = np.diff(lc.time.value[smask])
        dmask = dts > (1./(24))
        ndt = int(15./(24*60*dts.mean()))*2+1
        tmask[0:int(2*ndt)] = False # mask first 60 minutes of data

        for idx in np.argwhere(dmask).flatten(): 
            tmask[idx:idx+int(ndt)] = False

        nmask = ~np.isnan(lc.flux.value)
        time = lc.time.value[tmask&nmask]
        flux = lc.flux.value[tmask&nmask]

        # remove outliers using a median filter
        mflux = medfilt(flux, kernel_size=15)
        rflux = flux/mflux
        newflux, std = sigma_clip(rflux, dt=15)
        mask = np.isnan(newflux)

        # clip outliers
        time = time[~mask]
        flux = flux[~mask]

        # remove stellar variability
        pdur = 2*np.arctan(1/prior['pl_ratdor']) / (2*np.pi)
        dt = np.median(np.diff(np.sort(time)))
        wl = 2# float(prior['pl_orbper'])*0.75
        #flc, tlc = flatten(time, flux, window_length=wl, return_trend=True, method='biweight', robust=True)
        dflux = np.copy(flux)
        dtrend = np.ones(len(time))
        # piece wise flattening
        diff = np.diff(time)
        dmask = np.concatenate([~(diff > 1./24),[True]])
        dlabel,dcount = label(dmask)
        for j, d in enumerate(np.unique(dlabel)):
            if d == 0:
                continue
            dmask = dlabel == d

            # flatten
            flc, tlc = flatten(time[dmask], flux[dmask], window_length=wl, return_trend=True, method='biweight')
            dflux[dmask] = flc
            dtrend[dmask] = tlc
        
        # add to master list
        sv['time'].append(time)
        sv['flux'].append(dflux)
        sv['flux_err'].append((dtrend**0.5)/np.nanmedian(dtrend))
        sv['trend'].append(dtrend)
        sv['sector'].append(np.ones(len(time))*sector)

        # create plot with trend that was removed
        fig,ax = plt.subplots(1,figsize=(10,6))
        ax.plot(time, flux,'k.')
        ax.plot(time, dtrend,'r--')
        ax.set_title(f"{args.target} - Sector {sector}")
        ax.set_ylabel("PDC Flux")
        ax.set_xlabel("Time [TBJD]")
        ax.set_ylim([np.percentile(flux, 0.1), np.percentile(flux,99.9)])
        plt.tight_layout()
        #if not os.path.exists(os.path.join(planetdir, "lightcurves")):
            #os.makedirs(os.path.join(planetdir, "lightcurves"))
        plt.savefig( os.path.join(planetdir, planetname+f"_sector_{sector}_trend.png") )
        plt.close()

    # combine all sectors
    time = np.concatenate(sv['time'])
    alls = np.concatenate(sv['sector'])
    flux = np.concatenate(sv['flux'])
    trend = np.concatenate(sv['trend'])
    flux_err = np.concatenate(sv['flux_err'])

    # remove nans
    nanmask = np.isnan(flux) | np.isnan(flux_err)
    time = time[~nanmask]
    flux = flux[~nanmask]
    flux_err = flux_err[~nanmask]
    alls = alls[~nanmask]
    trend = trend[~nanmask]

    # create dataframe for entire light curve
    df = pd.DataFrame({'time':time, 'flux':flux*trend, 'flux_err':flux_err*trend, 'sector':alls})
    df.to_csv( os.path.join(planetdir, planetname+"_lightcurve.csv"), index=False)

    # fit transit for each epoch
    period = prior['pl_orbper']
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    pdur = 2*np.arctan(1/prior['pl_ratdor']) / (2*np.pi)
    ophase = (tphase+0.25)%1-0.25 # offset phase
    # mask off all the transits
    tmask = (ophase > -1*pdur) & (ophase < 1*pdur)

    # sigma clip time series
    dt = np.median(np.diff(np.sort(time)))
    ndt = int(15./(24*60*dt))*2+1 # number of data points in 15 minutes
    try:
        flux[tmask], _ = sigma_clip( flux[tmask], max(7,ndt), iterations=2)
    except:
        pass
    std = np.std(flux[tmask])

    # remove data points from sigma clip that are nan'ed and outliers
    nanmask = np.isnan(flux) | (flux > 1.03) | (flux < 0.94)

    flux = flux[~nanmask]
    time = time[~nanmask]
    alls = alls[~nanmask]

    # recompute transit mask with new time array 
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    ophase = (tphase+0.25)%1-0.25 # offset phase
    # mask off all the transits
    tmask = (ophase > -2*pdur) & (ophase < 2*pdur)

    # estimate uncertainties based on scatter out of transit
    phot_std = np.std(flux[~tmask])
    tpars = { 
        # transit 
        'rprs': ((prior['pl_radj']*u.R_jup)/(prior['st_rad']*u.R_sun)).to(u.m/u.m).value,
        'ars': float(prior['pl_ratdor']),
        'per': float(prior['pl_orbper']),
        'inc': float(prior['pl_orbincl']),
        'tmid': float(prior['pl_tranmid']), 

        # eclipse 
        'omega': float(prior['pl_orblper']), 
        'ecc': float(prior['pl_orbeccen']),

        'a1':0,
        'a2':0,

        # limb darkening (linear, quadratic)
        'u0':u0, 'u1':u1, 'u2':u2, 'u3':u3
    }

    # check for nan in omega 
    if np.isnan(tpars['omega']):
        tpars['omega'] = 0
    if np.isnan(tpars['ecc']):
        tpars['ecc'] = 0

    # check for nans and replace with 0
    for key in prior.keys():
        try:
            if np.isnan(prior[key]):
                prior[key] = 0
                print(f"WARNING: {key} is nan, setting to 0")
        except:
            # stupid numpy, can't check if str is nan
            pass

    # create bounds for light curve fitting
    mybounds = {
        'rprs':[0,3*tpars['rprs']],
        #'per':[tpars['per']*0.999, tpars['per']*1.001],
        'tmid':[tpars['tmid']-0.1, tpars['tmid']+0.1],
        'inc':[90 - np.rad2deg(np.arctan(1./tpars['ars']))-1,90],
    }

    if args.ars:
        mybounds['ars'] = [tpars['ars']*0.01, tpars['ars']*3]
        # degenerate if fit with Rp/Rs and Inclination but could try anyways


    # prepare for individual fits
    period = prior['pl_orbper']
    tmid = prior['pl_tranmid']
    ars = prior['pl_ratdor']

    tphase = (time + 2457000.0 - tmid) / period
    pdur = 2*np.arctan(1/ars) / (2*np.pi)
    events = np.unique(np.floor(tphase))

    sv['timeseries'] = {
        'time': time,
        'flux': flux,
        'flux_err': np.copy(phot_std/flux),
        'residuals': np.copy(flux)
    }

    # for each individual transit
    for e in events:

        # mask data in phase
        tmask = (tphase > e - (2*pdur)) & (tphase < e + (2*pdur) )

        if tmask.sum() == 0:
            continue

        # incomplete transit
        #if (time[tmask].max() - time[tmask].min()) < (0.5*(4*pdur)):
        #    continue

        # update priors for fitting
        tpars['tmid'] = tmid + e*period
        tpars['per'] = period
        tpars['ars'] = ars

        # bounds for individual fits
        mybounds = {
            'rprs':[0,3*tpars['rprs']],
            'tmid':[tpars['tmid']-0.25*pdur*tpars['per'], tpars['tmid']+0.25*pdur*tpars['per']],
            'inc':[tpars['inc']-5, min(90,tpars['inc']+5)],
            #'ars':[tpars['ars']*0.9, tpars['ars']*2.1]
        }

        # skip light curves with large gaps bigger than 35 minutes
        if np.sum(np.diff(np.sort(time[tmask]))*24*60 > 35):
            continue

        # fit data
        try:
            airmass = np.zeros(len(time[tmask]))
            myfit = lc_fitter(time[tmask]+2457000.0, flux[tmask], phot_std/flux[tmask], airmass, tpars, mybounds)
        except:
            print(f"Failed to fit transit at phase: {e}")
            continue

        

        rprs2 = myfit.parameters['rprs']**2
        rprs2err = myfit.parameters['rprs']*2*myfit.errors['rprs']
        depth = myfit.transit.max() - myfit.transit.min()

        # don't add transits consistent with 0
#        if (rprs2-3*rprs2err <= 0):
        if (rprs2-1*rprs2err <= 0):
            print(f"Skipping transit at phase: {e} b.c. depth = {rprs2:.2e} +- {rprs2err:.2e}")
            continue

        # remove partial transits
        if myfit.transit[0] < 1 or myfit.transit[-1] < 1:
           continue

        # save results
        lcdata = {
            'time':myfit.time,
            'flux':myfit.data,
            'fluxerr':myfit.dataerr,
            'phase':myfit.phase,
            'residuals':myfit.residuals,
            'pars':myfit.parameters,
            'errors':myfit.errors,
            'rchi2':myfit.chi2/len(myfit.time),
            'quality':myfit.quality,
            'sector':int(np.median(alls[tmask])),
            'partial_ratio':myfit.duration_measured / myfit.duration_expected,
            'snr':rprs2/rprs2err,
            'snr2':depth/np.std(myfit.residuals)
        }

        # save results to state vector
        sv['lightcurves'].append(lcdata)

        # detrend timeseries
        sv['timeseries']['residuals'][tmask] /= myfit.transit

        # create string for plot legend
        tmidstr = str(np.round(myfit.parameters['tmid'],2)).replace('.','_')

        # save bestfit
        fig,ax = myfit.plot_bestfit(title=f"{args.target} - Sector {lcdata['sector']}", bin_dt=0.5/24.)
        plt.savefig( os.path.join(planetdir, f"{tmidstr}_"+planetname+"_lightcurve.png") )
        plt.close()

        # save posterior
        fig = myfit.plot_triangle()
        plt.savefig( os.path.join(planetdir, f"{tmidstr}_"+planetname+"_posterior.png") )
        plt.close()

        csv_data = {
            'Time (JD)':myfit.time,
            'Relative Flux':myfit.data,
            'Relative Flux Error':myfit.dataerr,
        }

        # create AAVSO formatted csv
        csv_lk = OutputFiles(myfit, prior, infoDict, planetdir)
        try:
            csv_lk.aavso_csv(airmass,u0,u1,u2,u3,tmidstr)
            csv_lk.aavso(airmass,u0,u1,u2,u3, tmidstr)
        except:
            print(f"Failed to create AAVSO csv for {args.target} - Sector {lcdata['sector']}")

    # save sv pickle
    pickle.dump(sv, open(os.path.join(planetdir, planetname+"_data.pkl"),"wb"))
        
    # O-C plot
    tmids = np.array([lc['pars']['tmid'] for lc in sv['lightcurves']])
    tmiderr = np.array([lc['errors']['tmid'] for lc in sv['lightcurves']])
    sectors = np.array([lc['sector'] for lc in sv['lightcurves']])
    ratios = np.array([lc['partial_ratio'] for lc in sv['lightcurves']])

    # mask out partial transits
    dmask = ratios > 0.75
    tmids = tmids[dmask]
    tmiderr = tmiderr[dmask]
    sectors = sectors[dmask]
    ratios = ratios[dmask]

    # TODO finish making O-C plot?
    # use example from exotic.api.nested_linear_fitter