#Federico R Noguer 2024
#Created for candidate exoplanets from tess.py
#to Run
# cd /EXOTIC/examples/tess/candidates
# python toi_individ_lc.py -i <input_inits_file_path> -o <output_directory_path>

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

def sigma_clip(ogdata, dt, iterations=1):
    mask = np.ones(ogdata.shape, dtype=bool)
    for i in range(iterations):
        mdata = savgol_filter(ogdata[mask], dt, 2)
        res = ogdata[mask] - mdata
        try:
            std = np.nanmedian([np.nanstd(np.random.choice(res, 100)) for i in range(250)])
        except:
            std = np.nanstd(res)
        mask[mask] = np.abs(res) < 3 * std

    mdata = savgol_filter(ogdata[mask], dt, 2)
    res = ogdata[mask] - mdata

    data = copy.deepcopy(ogdata)
    data[~mask] = np.nan
    return data, np.std(res)

def parse_args():
    parser = argparse.ArgumentParser()

    help_ = "Path to the initialization file"
    parser.add_argument("-i", "--init", help=help_, type=str, required=True)

    help_ = "Directory for saving results"
    parser.add_argument("-o", "--output", help=help_, default="output/", type=str)

    help_ = "Sectors (0 = all)"
    parser.add_argument("-s", "--sector", help=help_, default=0, type=int)

    help_ = "Choose a target to re-process. Need to do this when after changing the prior."
    parser.add_argument("-r", "--reprocess", help=help_, action='store_true', default=False)

    help_ = "Add a/Rs to the fitting process"
    parser.add_argument("--ars", help=help_, action='store_true', default=False)

    help_ = "Perform transit least squares search over residuals of global fit"
    parser.add_argument("--tls", help=help_, action='store_true', default=False)

    args = parser.parse_args()

    # Ensure the init file path is correctly processed
    if not os.path.isfile(args.init):
        raise FileNotFoundError(f"Initialization file not found: {args.init}")

    return args

def check_std(time, flux, dt=0.5):  # dt = [hr]
    tdt = np.diff(np.sort(time)).mean()
    si = np.argsort(time)
    sflux = savgol_filter(flux[si], 1 + 2 * int(max(15, dt / 24 / tdt)), 2)
    return np.nanstd(flux - sflux)

# compute stellar mass from logg and radius
stellar_mass = lambda logg, rs: ((rs * u.R_sun) ** 2 * 10 ** logg * (u.cm / u.s ** 2) / const.G).to(u.M_sun).value

# keplerian semi-major axis (au)
sa = lambda m, P: ((const.G * m * u.M_sun * P * u.day ** 2 / (4 * np.pi ** 2)) ** (1. / 3)).to(u.AU).value

if __name__ == "__main__":
    infoDict = {}
    args = parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output, exist_ok=True)

    # Load initialization file
    with open(args.init, 'r') as file:
        prior = json.load(file)

    # Use the planet name from the initialization file
    planetname = prior['pl_name'].lower().replace(' ', '').replace('-', '_')
    planetdir = os.path.join(args.output, planetname)

    if not os.path.exists(os.path.join(planetdir, "global_fit.png")):
        os.makedirs(planetdir, exist_ok=True)
    elif not args.reprocess:
        raise Exception("Target already processed. Use -r to reprocess.")

    # Use the planet name for searching light curve data
    search_result = lk.search_targetpixelfile(prior['pl_name'], mission='TESS')
    print(search_result)

    if len(search_result) == 0:
        raise Exception(f"No data for: {prior['pl_name']}")

    ## load prior from JSON file
    #json_file_path = os.path.join(planetdir, planetname + "_prior.json")
    #if os.path.exists(json_file_path):
    #    with open(json_file_path, 'r') as file:
    #        prior = json.load(file)
    #else:
    #    raise(Exception(f"JSON file not found for: {prior['pl_name']}"))

    if len(prior) == 0:
        raise(Exception(f"no target named: {prior['pl_name']}"))

    # find sectors
    sectors = [int(sector.split(' ')[-1]) for sector in search_result.mission]
    if args.sector:
        if args.sector not in sectors:
            raise(Exception(f"no data for sector: {args.sector}"))

    # compute limb darkening
    if np.isnan(prior['st_meterr1']):
        prior['st_meterr1'] = prior['st_met'] * 0.1
        prior['st_meterr2'] = -prior['st_met'] * 0.1

    u0, u1, u2, u3 = exotethys(prior['st_logg'], prior['st_teff'], prior['st_met'], 'TESS', method='claret', stellar_model='phoenix')

    # alloc state vector
    sv = {
        'lightcurves': [],  # dicts of time, flux, fluxerr, pars, errors from light curve fit
        'sectors': {},  # exp time of each sector that gets processed

        # global timeseries
        'time': [],
        'flux': [],
        'flux_err': [],
        'trend': [],
        'sector': []
    }

    # loop over data sets
    for i, sector in enumerate(sectors):

        # just process single sector
        if args.sector:
            if sector != args.sector:
                continue

        # don't reprocess same sector at lower res
        if search_result.exptime[i].value > sv['sectors'].get(sector, 1800):
            print(f" Skipping Sector {sector} @ {search_result.exptime[i].value} s")
            print(f"  Already processed @ {sv['sectors'][sector]} s")
            continue
        if i == 0:
            exptime = search_result.exptime[i].value
        if exptime < search_result.exptime[i].value:
            exptime = search_result.exptime[i].value

        print(f"Downloading Sector {sector} @ {search_result.exptime[i].value} s ...")
        try:
            tpf = search_result[i].download(quality_bitmask='hard')
        except Exception as err:
            print(f" Failed to download {search_result.exptime[i].value} s data for sector {sector}")
            print(f" {err}")
            continue

        # store exp time
        sv['sectors'][sector] = search_result.exptime[i].value
        infoDict['exposure'] = search_result.exptime[i].value

        # aper selection
        lc = tpf.to_lightcurve(aperture_mask=tpf.pipeline_mask)
        nmask = np.isnan(lc.flux.value)
        lstd = check_std(lc.time.value[~nmask], lc.flux.value[~nmask])

        # test two more apertures, slightly bigger each iteration
        for it in [1, 2]:
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
        plt.savefig(os.path.join(planetdir, planetname + f"_sector_{sector}_aperture.png"))
        plt.close()

        # remove first ~30 min of data after any big gaps
        tmask = np.ones(lc.time.shape).astype(bool)
        smask = np.argsort(lc.time.value)
        dts = np.diff(lc.time.value[smask])
        dmask = dts > (1. / (24))
        ndt = int(15. / (24 * 60 * dts.mean())) * 2 + 1
        tmask[0:int(2 * ndt)] = False  # mask first 60 minutes of data

        for idx in np.argwhere(dmask).flatten():
            tmask[idx:idx + int(ndt)] = False

        nmask = ~np.isnan(lc.flux.value)
        time = lc.time.value[tmask & nmask]
        flux = lc.flux.value[tmask & nmask]

        # remove outliers using a median filter
        mflux = medfilt(flux, kernel_size=15)
        rflux = flux / mflux
        newflux, std = sigma_clip(rflux, dt=15)
        mask = np.isnan(newflux)

        # clip outliers
        time = time[~mask]
        flux = flux[~mask]

        # remove stellar variability
        pdur = 2 * np.arctan(1 / prior['pl_ratdor']) / (2 * np.pi)
        dt = np.median(np.diff(np.sort(time)))
        wl = 2  # float(prior['pl_orbper'])*0.75
        dflux = np.copy(flux)
        dtrend = np.ones(len(time))
        # piece wise flattening
        diff = np.diff(time)
        dmask = np.concatenate([~(diff > 1. / 24), [True]])
        dlabel, dcount = label(dmask)
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
        sv['flux_err'].append((dtrend ** 0.5) / np.nanmedian(dtrend))
        sv['trend'].append(dtrend)
        sv['sector'].append(np.ones(len(time)) * sector)

        # create plot with trend that was removed
        fig, ax = plt.subplots(1, figsize=(10, 6))
        ax.plot(time, flux, 'k.')
        ax.plot(time, dtrend, 'r--')
        ax.set_title(f"{prior['pl_name']} - Sector {sector}")
        ax.set_ylabel("PDC Flux")
        ax.set_xlabel("Time [TBJD]")
        ax.set_ylim([np.percentile(flux, 0.1), np.percentile(flux, 99.9)])
        plt.tight_layout()
        plt.savefig(os.path.join(planetdir, planetname + f"_sector_{sector}_trend.png"))
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
    df = pd.DataFrame({'time': time, 'flux': flux * trend, 'flux_err': flux_err * trend, 'sector': alls})
    df.to_csv(os.path.join(planetdir, planetname + "_lightcurve.csv"), index=False)

    # fit transit for each epoch
    period = prior['pl_orbper']
    tphase = (time + 2457000.0 - prior['pl_tranmid']) / period
    pdur = 2 * np.arctan(1 / prior['pl_ratdor']) / (2 * np.pi)
    ophase = (tphase + 0.25) % 1 - 0.25  # offset phase
    # mask off all the transits
    tmask = (ophase > -1 * pdur) & (ophase < 1 * pdur)

    # sigma clip time series
    dt = np.median(np.diff(np.sort(time)))
    ndt = int(15. / (24 * 60 * dt)) * 2 + 1  # number of data points in 15 minutes
    try:
        flux[tmask], _ = sigma_clip(flux[tmask], max(7, ndt), iterations=2)
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
    ophase = (tphase + 0.25) % 1 - 0.25  # offset phase
    # mask off all the transits
    tmask = (ophase > -2 * pdur) & (ophase < 2 * pdur)

    # estimate uncertainties based on scatter out of transit
    phot_std = np.std(flux[~tmask])
    tpars = {
        # transit
        'rprs': ((prior['pl_radj'] * u.R_jup) / (prior['st_rad'] * u.R_sun)).to(u.m / u.m).value,
        'ars': float(prior['pl_ratdor']),
        'per': float(prior['pl_orbper']),
        'inc': float(prior['pl_orbincl']),
        'tmid': float(prior['pl_tranmid']),

        # eclipse
        'omega': float(prior['pl_orblper']),
        'ecc': float(prior['pl_orbeccen']),

        'a1': 0,
        'a2': 0,

        # limb darkening (linear, quadratic)
        'u0': u0, 'u1': u1, 'u2': u2, 'u3': u3
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
            pass

    # create bounds for light curve fitting
    mybounds = {
        'rprs': [0, 3 * tpars['rprs']],
        'tmid': [tpars['tmid'] - 0.1, tpars['tmid'] + 0.1],
        'inc': [80, 90],
    }

    if args.ars:
        mybounds['ars'] = [tpars['ars'] * 0.01, tpars['ars'] * 3]

    # prepare for individual fits
    period = prior['pl_orbper']
    tmid = prior['pl_tranmid']
    ars = prior['pl_ratdor']

    tphase = (time + 2457000.0 - tmid) / period
    pdur = 2 * np.arctan(1 / ars) / (2 * np.pi)
    events = np.unique(np.floor(tphase))

    sv['timeseries'] = {
        'time': time,
        'flux': flux,
        'flux_err': np.copy(phot_std / flux),
        'residuals': np.copy(flux)
    }

    # for each individual transit
    for e in events:

        # mask data in phase
        tmask = (tphase > e - (2 * pdur)) & (tphase < e + (2 * pdur))

        if tmask.sum() == 0:
            continue

        # update priors for fitting
        tpars['tmid'] = tmid + e * period
        tpars['per'] = period
        tpars['ars'] = ars

        # bounds for individual fits
        mybounds = {
            'rprs': [0, 3 * tpars['rprs']],
            'tmid': [tpars['tmid'] - 0.25 * pdur * tpars['per'], tpars['tmid'] + 0.25 * pdur * tpars['per']],
            'inc': [tpars['inc'] - 5, min(90, tpars['inc'] + 5)],
        }

        # skip light curves with large gaps bigger than 35 minutes
        if np.sum(np.diff(np.sort(time[tmask])) * 24 * 60 > 35):
            continue

        # fit data
        try:
            airmass = np.zeros(len(time[tmask]))
            myfit = lc_fitter(time[tmask] + 2457000.0, flux[tmask], phot_std / flux[tmask], airmass, tpars, mybounds)
        except:
            print(f"Failed to fit transit at phase: {e}")
            continue

        for k in myfit.bounds.keys():
            print("{:.6f} +- {}".format(myfit.parameters[k], myfit.errors[k]))

        rprs2 = myfit.parameters['rprs'] ** 2
        rprs2err = myfit.parameters['rprs'] * 2 * myfit.errors['rprs']
        depth = myfit.transit.max() - myfit.transit.min()

#        if (rprs2 - 1 * rprs2err <= 0):
#            print(f"Skipping transit at phase: {e} b.c. depth = {rprs2:.2e} +- {rprs2err:.2e}")
#            continue

        if myfit.transit[0] < 1 or myfit.transit[-1] < 1:
            continue
        print ("depth =",depth)
        print ("duration_measured =",myfit.duration_measured)
        print ("duration_expected =",myfit.duration_expected)
        print ("Tmid =",myfit.parameters['tmid'])
        #print ("Sector =",int(np.median(alls[tmask]))
        print ("Phase =",myfit.phase)
        print ("SNR = ",rprs2 / rprs2err)
        # save results
        lcdata = {
            'time': myfit.time,
            'flux': myfit.data,
            'fluxerr': myfit.dataerr,
            'phase': myfit.phase,
            'residuals': myfit.residuals,
            'pars': myfit.parameters,
            'errors': myfit.errors,
            'rchi2': myfit.chi2 / len(myfit.time),
            'quality': myfit.quality,
            'sector': int(np.median(alls[tmask])),
            'partial_ratio': myfit.duration_measured / myfit.duration_expected,
            'snr': rprs2 / rprs2err,
            'snr2': depth / np.std(myfit.residuals)
        }

        sv['lightcurves'].append(lcdata)

        sv['timeseries']['residuals'][tmask] /= myfit.transit

        tmidstr = str(np.round(myfit.parameters['tmid'], 2)).replace('.', '_')

        fig, ax = myfit.plot_bestfit(title=f"{prior['pl_name']} - Sector {lcdata['sector']}", bin_dt=0.5 / 24.)
        plt.savefig(os.path.join(planetdir, f"{tmidstr}_" + planetname + "_lightcurve.png"))
        plt.close()

        fig = myfit.plot_triangle()
        plt.savefig(os.path.join(planetdir, f"{tmidstr}_" + planetname + "_posterior.png"))
        plt.close()

        csv_data = {
            'Time (JD)': myfit.time,
            'Relative Flux': myfit.data,
            'Relative Flux Error': myfit.dataerr,
        }

        csv_lk = OutputFiles(myfit, prior, infoDict, planetdir)
        csv_lk.aavso_csv(airmass, u0, u1, u2, u3, tmidstr)
        csv_lk.aavso(airmass, u0, u1, u2, u3, tmidstr)

    pickle.dump(sv, open(os.path.join(planetdir, planetname + "_data.pkl"), "wb"))
