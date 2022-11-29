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
from scipy.signal import savgol_filter
from astropy import constants as const
from astropy import units as u
from astropy.time import Time
from scipy.ndimage import label
import statsmodels.api as sm
from pylightcurve import exotethys
import lightkurve as lk
from wotan import flatten
from transitleastsquares import transitleastsquares
from ldtk import LDPSetCreator, BoxcarFilter
from ldtk.filters import create_tess
from pylightcurve.models.exoplanet_lc import planet_orbit
from exotic.api.elca import transit, lc_fitter
from exotic.exotic import LimbDarkening
from tools import lomb_scargle
from nested_linear_fitter import linear_fitter
from output_aavso import OutputFiles
# import pdb; pdb.set_trace
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
    parser.add_argument("-o", "--output", help=help_, default="output/", type=str)

    help_ = "Sectors (0 = all)"
    parser.add_argument("-s", "--sector", help=help_, default=0, type=int)

    help_ = "Choose a target to re-process"
    parser.add_argument("-r", "--reprocess", help=help_, action='store_true', default=False)

    help_ = "Add a/Rs to the fitting process"
    parser.add_argument("--ars", help=help_, action='store_true', default=False)

    return parser.parse_args()

def check_std(time, flux, dt=0.5): # dt = [hr]
    tdt = np.diff(np.sort(time)).mean()
    si = np.argsort(time)
    sflux = savgol_filter(flux[si], 1+2*int(max(15,dt/24/tdt)), 2)
    return np.nanstd(flux - sflux)


stellar_mass = lambda logg,rs: ((rs*u.R_sun)**2 * 10**logg*(u.cm/u.s**2) / const.G).to(u.M_sun).value

# # keplerian semi-major axis (au)
sa = lambda m,P : ((const.G*m*u.M_sun*P*u.day**2/(4*np.pi**2) )**(1./3)).to(u.AU).value

# # approximate radial velocity semi-amplitude 
# rv_semi = lambda M,m,a : (G*m*m/(M*a))**0.5

# # rv_semi( 1.22, 90*mearth/msun, sa(1.22,2.155)) * au/(24*60*60)

# # assume a relationship between stellar mass and radius from 
# # (http://ads.harvard.edu/books/hsaa/toc.html)
# mms = [40,18, 6.5,3.2,2.1,1.7,1.3,1.1 ,1,0.93,0.78,0.69,0.47,0.21,0.1]  # M/Msun
# rrs = [18,7.4,3.8,2.5,1.7,1.3,1.2,1.05,1,0.93,0.85,0.74,0.63,0.32,0.13] # R/Rsun
# m2r_star = interp1d(mms,rrs) # convert solar mass to solar radius


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
    
    search_result = lk.search_targetpixelfile(args.target, mission='TESS')
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
                "st_meterr1": float(row[' Stellar Metallicity err']),
                "st_meterr2": float(-row[' Stellar Metallicity err']),
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
            import pdb; pdb.set_trace()
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
        prior["pl_name"] = planetname_space
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

        # remove stellar variability
        pdur = 2*np.arctan(1/prior['pl_ratdor']) / (2*np.pi)
        dt = np.median(np.diff(np.sort(time)))
        wl = 2
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
            flc, tlc = flatten(time[dmask], flux[dmask], window_length=wl, return_trend=True, method='biweight', robust=True)
            dflux[dmask] = flc
            dtrend[dmask] = tlc
        
        # add to master list
        sv['time'].append(time)
        sv['flux'].append(dflux)
        sv['trend'].append(dtrend)
        sv['sector'].append(np.ones(len(time))*sector)


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

    time = np.concatenate(sv['time'])
    flux = np.concatenate(sv['flux'])
    alls = np.concatenate(sv['sector'])

    # fit transit for each epoch
    period = prior['pl_orbper']
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    pdur = 2*np.arctan(1/prior['pl_ratdor']) / (2*np.pi)
    ophase = (tphase+0.25)%1-0.25 # offset phase
    # mask off all the transits
    tmask = (ophase > -1*pdur) & (ophase < 1*pdur)

    # sigma clip time series
    dt = np.median(np.diff(np.sort(time)))
    ndt = int(15./(24*60*dt))*2+1
    flux[tmask], _ = sigma_clip( flux[tmask], max(7,ndt), iterations=2)
    std = np.std(flux[tmask])
    nanmask = np.isnan(flux) | (flux > 1.03) | (flux < 0.94)

    flux = flux[~nanmask]
    time = time[~nanmask]
    alls = alls[~nanmask]
    # recompute transit mask with new time array 
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    ophase = (tphase+0.25)%1-0.25 # offset phase
    # mask off all the transits
    tmask = (ophase > -3*pdur) & (ophase < 3*pdur)

    # estimate uncertainties based on scatter out of transit
    phot_std = np.std(flux[~tmask])
    tpars = { 
        # transit 
        'rprs': ((prior['pl_radj']*u.R_jup)/(prior['st_rad']*u.R_sun)).to(u.m/u.m).value,
        'ars': float(prior['pl_ratdor']),
        'per': float(prior['pl_orbper']),
        'inc': float(prior['pl_orbincl']),
        'tmid': float(prior['pl_tranmid']), 
        #'tmid': np.median(time[tmask])+2457000.0, 


        # eclipse 
        'omega': float(prior['pl_orblper']), 
        'ecc': float(prior['pl_orbeccen']),

        'a1':0,
        'a2':0,

        # limb darkening (linear, quadratic)
        'u0':u0, 'u1':u1, 'u2':u2, 'u3':u3
    }

    mybounds = {
        'rprs':[0,3*tpars['rprs']],
        #'per':[tpars['per']*0.999, tpars['per']*1.001],
        'tmid':[tpars['tmid']-0.1, tpars['tmid']+0.1],
        'inc':[90 - np.rad2deg(np.arctan(1./tpars['ars']))-1,90],
    }

    if args.ars:
        #mybounds['omega'] = [tpars['omega']-1, tpars['omega']+1]
        mybounds['ars'] = [tpars['ars']*0.01, tpars['ars']*3]
        #mybounds['ars'] = [80,100]#[tpars['ars']*0.01, tpars['ars']*3]

        #del mybounds['inc']

    # import pdb; pdb.set_trace()
    print('performing global fit...')
    #myfit = lc_fitter(time[tmask]+2457000.0, flux[tmask], phot_std/flux[tmask], tpars, mybounds, verbose=True)
    #import pdb ; pdb.set_trace()
    airmass = np.zeros(len(time[tmask]))
    myfit = lc_fitter(time[tmask]+2457000.0, flux[tmask], phot_std/flux[tmask], airmass, tpars, mybounds, verbose=True)

    myfit.plot_bestfit(title=f"{args.target} Global Fit")
    plt.savefig( os.path.join( planetdir, planetname+"_global_fit.png"))
    plt.close()

    myfit.plot_triangle()
    plt.savefig( os.path.join( planetdir, planetname+"_global_triangle.png"))
    plt.close()

    # update priors
    for k in myfit.errors.keys():
        tpars[k] = myfit.parameters[k]

    # update priors
    try:
        prior['pl_orbper'] = myfit.parameters['per']
        prior['pl_orbpererr1'] = myfit.errors['per']
        prior['pl_orbpererr2'] = -myfit.errors['per']
    except:
        pass

    prior['pl_tranmid'] = myfit.parameters['tmid']
    prior['pl_tranmiderr1'] = myfit.errors['tmid']
    prior['pl_tranmiderr2'] = -myfit.errors['tmid']

    if not args.ars:
        prior['pl_orbincl'] = myfit.parameters['inc']
        prior['pl_orbinclerr1'] = myfit.errors['inc']
        prior['pl_orbinclerr2'] = -myfit.errors['inc']

    # save prior to disk
    with open(os.path.join(planetdir,planetname+"_prior.json"), 'w', encoding ='utf8') as json_file:
        json.dump(prior, json_file, indent=4)

    # save results to state vector
    sv['global'] = {
        'pars': myfit.parameters,
        'errs': myfit.errors,
        'time': myfit.time,
        'data': myfit.data,
        'data_err': myfit.dataerr,
    }

    rprs2 = myfit.parameters['rprs']**2
    rprs2err = myfit.parameters['rprs']*2*myfit.errors['rprs']
    snr = rprs2/rprs2err

    if snr < 4:
        with open("notes.txt", 'w') as f:
            f.write(f"Skipping individual light curve fits b.c SNR = {snr:.2f}")
        # save global fit data
        pickle.dump(sv, open(os.path.join(planetdir, planetname+"_data.pkl"),"wb"))
        raise(Exception(f"Skipping individual light curve fits b.c SNR = {snr:.2f}"))

    # fit individual light curves
    #lightcurvedir = os.path.join(planetdir,'lightcurves')
    #if not os.path.exists(lightcurvedir):
        #os.mkdir(lightcurvedir)

    period = myfit.parameters['per']
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    pdur = 2*np.arctan(1/prior['pl_ratdor']) / (2*np.pi)
    events = np.unique(np.floor(tphase))
    
    for e in events:
        tmask = (tphase > e - (3*pdur)) & (tphase < e + (3*pdur) )

        if tmask.sum() == 0:
            continue

        # incomplete transit
        if (time[tmask].max() - time[tmask].min()) < (0.75*(4*pdur)):
            continue

        tpars['tmid'] = prior['pl_tranmid'] + e*period

        mybounds = {
            'rprs':[0,3*tpars['rprs']],
            'tmid':[tpars['tmid']-0.25*pdur*tpars['per'], tpars['tmid']+0.25*pdur*tpars['per']],
            'ars':[tpars['ars']*0.9, tpars['ars']*1.1]
        }

    
        # # estimate uncertainties based on scatter
        # dt = np.median(np.diff(time[tmask]))
        # sflc = savgol_filter(flux[tmask], my1+2*int(0.75/24/dt), 2)
        # res = flux[tmask] - sflc
        # phot_std = np.std(res)
        # tmask[tmask] = (np.abs(res) <3*phot_std)

        # skip light curves with large gaps
        if np.sum(np.diff(np.sort(time[tmask]))*24*60 > 35):
            continue

        # fit data
        try:
            airmass = np.zeros(len(time[tmask]))
            myfit = lc_fitter(time[tmask]+2457000.0, flux[tmask], phot_std/flux[tmask], airmass, tpars, mybounds)
        except:
            continue

        for k in myfit.bounds.keys():
            print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

        rprs2 = myfit.parameters['rprs']**2
        rprs2err = myfit.parameters['rprs']*2*myfit.errors['rprs']
        depth = myfit.transit.max() - myfit.transit.min()

        # don't add transits consistent with 0
        if (rprs2-rprs2err < 0):
            continue

        # remove partial transits
        if myfit.transit[0] < 1 or myfit.transit[-1] < 1:
            continue

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
        #NEED TO ALSO PRINT OUT AAVSO META DATA
        csv_lk = OutputFiles(myfit, prior, infoDict,planetdir)
        csv_lk.aavso_csv(airmass,u0,u1,u2,u3,tmidstr)
        csv_lk.aavso(airmass,u0,u1,u2,u3, tmidstr)

        #pd.DataFrame(csv_data).to_csv( os.path.join(planetdir, f"{tmidstr}_"+planetname+"_lightcurve.csv"), index=False )

    # O-C plot
    tmids = np.array([lc['pars']['tmid'] for lc in sv['lightcurves']])
    tmiderr = np.array([lc['errors']['tmid'] for lc in sv['lightcurves']])
    sectors = np.array([lc['sector'] for lc in sv['lightcurves']])
    ratios = np.array([lc['partial_ratio'] for lc in sv['lightcurves']])

    #duration mask
    dmask = ratios > 0.75
    tmids = tmids[dmask]
    tmiderr = tmiderr[dmask]
    sectors = sectors[dmask]
    ratios = ratios[dmask]

    pickle.dump(sv, open(os.path.join(planetdir, planetname+"_data.pkl"),"wb"))

    assert(len(tmids) > 1)

    # linear fit to transit times
    epochs = np.round((tmids - np.max(tmids)) / period)
    A = np.vstack([np.ones(len(tmids)), epochs]).T
    res = sm.WLS(tmids, A, weights=1 / tmiderr).fit()
    ephem = res.params[0] + res.params[1] * epochs
    oc = tmids - ephem

    # 3 sigma clip
    meddiff = np.median(np.abs(np.diff(oc)))
    maxdiff = np.max(np.abs(np.diff(oc)))
    mask = np.diff(oc) > (meddiff + 3* np.mean(tmiderr))
    mask = ~np.concatenate([mask,[False]])

    epochs = np.round((tmids[mask] - np.max(tmids[mask])) / period)
    A = np.vstack([np.ones(len(tmids[mask])), epochs]).T
    res = sm.WLS(tmids[mask], A, weights=1 / tmiderr[mask]).fit()
    ephem = res.params[0] + res.params[1] * epochs
    oc = tmids[mask] - ephem


    bounds = {
        # maintain the order: m,b
        'm': [res.params[1]-5.*max(tmiderr), res.params[1]+5.*max(tmiderr)],
        'b': [res.params[0]-5.*max(tmiderr), res.params[0]+5.*max(tmiderr)]
    }
    lf = linear_fitter(epochs, tmids[mask], tmiderr[mask], bounds)

    sv['ephemeris'] = {
        'tmids':list(tmids[mask]),
        'epochs':list(epochs),
        'tmids_err':list(tmiderr[mask]),
        'tmid': lf.parameters['b'],
        'tmid_err': lf.errors['b'],
        'per': lf.parameters['m'],
        'per_err': lf.errors['m'],
        'residuals': list(lf.residuals),
    }

    # save to disk
    json.dump(sv['ephemeris'], open(os.path.join(planetdir,planetname+"_ephemeris.json"),"w"), indent=4 )
    #pickle.dump(sv, open(os.path.join(planetdir, planetname+"_data.pkl"),"wb"))
    #pickle.dump([sectors,mask,tmids, lf.time, lf.residuals, lf.dataerr, lf.parameters['m'], lf.errors['m'], lf.parameters['b'], lf.errors['b'], airmass, u0,u1,u2,u3, myfit, prior, infoDict], open(os.path.join(planetdir, planetname+"_extraPdata.pkl"),"wb"))
    pickle.dump([mask,tmids,epochs,lf.residuals,lf.dataerr,], open(os.path.join(planetdir, planetname+"_data.pkl"),"wb"))
    # OC figure
    fig, ax = plt.subplots(1, ncols=len(set(sectors)), figsize=(13,5))
    plt.subplots_adjust(hspace=0,wspace=0,left=0.08,bottom=0.13,right=0.975)
    for i, sector in enumerate(np.sort(list(set(sectors)))):
        smask = sectors[mask] == sector

        # dates for header
        dmin = Time(tmids[mask][smask],format='jd').min().isot.split('T')[0]
        dmax = Time(tmids[mask][smask],format='jd').max().isot.split('T')[0]
        try:
            ax[i].errorbar(lf.time[smask], lf.residuals[smask]*24*60, yerr=lf.dataerr[smask]*24*60, ls='none', marker='o',label='Data')
            ax[i].axhline(0,ls='--',label="P: {:.5f}+-{:.5f}, Tmid: {:.5f}+-{:.5f}".format(lf.parameters['m'], lf.errors['m'], lf.parameters['b'], lf.errors['b']))
            ax[i].set_ylim([min(lf.residuals*24*60)-1, max(lf.residuals*24*60)+1])
            ax[i].grid(True, ls='--')
            ax[i].set_title(f"S{sector}")
            x_ticks = ax[i].xaxis.get_major_ticks()
            x_ticks[0].label1.set_visible(False)
            x_ticks[-1].label1.set_visible(False)
            for tick in ax[i].get_xticklabels():
                tick.set_rotation(45)
            if i > 0:
                ax[i].yaxis.set_ticklabels([])
        except:
            ax.errorbar(lf.time[smask], lf.residuals[smask]*24*60, yerr=lf.dataerr[smask]*24*60, ls='none', marker='o',label='Data')
            ax.axhline(0,ls='--',label="P: {:.5f}+-{:.5f}, Tmid: {:.5f}+-{:.5f}".format(lf.parameters['m'], lf.errors['m'], lf.parameters['b'], lf.errors['b']))
            ax.set_ylim([min(lf.residuals*24*60)-1, max(lf.residuals*24*60)+1])
            ax.set_xlabel("Epoch")
            ax.grid(True, ls='--')
            ax.set_title(f"Sector {sector} ({dmin} - {dmax})")

    try:
        ax[0].set_ylabel("Observed - Calculated [min]")
        ax[int(0.5*len(ax))].set_xlabel("Epoch")

    except:
        ax.set_ylabel("Observed - Calculated [min]")
    plt.savefig(os.path.join(planetdir, planetname+"_ephemeris.png") )
    plt.close()
