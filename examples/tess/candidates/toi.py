#Federico R Noguer 2024
#Created for candidate exoplanets from tess.py
#to Run
# cd /EXOTIC/examples/tess/candidates
# python toi.py -i <input_inits_file_path> -o <output_directory_path>

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
    # Table Access Protocol (TAP) query
    uri_full = base_url
    for k in query:
        if k != "format":
            uri_full += f"{k} {query[k]} "
    uri_full = f"{uri_full[:-1]} &format={query.get('format', 'csv')}"
    uri_full = uri_full.replace(' ', '+')
    print(uri_full)
    response = requests.get(uri_full, timeout=300)
    if dataframe:
        return read_csv(StringIO(response.text))
    else:
        return response.text

def nea_scrape(target=None):
    uri_ipac_base = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="
    uri_ipac_query = {
        "select": "pl_name,hostname,pl_radj,pl_radjerr1,ra,dec,"
                  "pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,"
                  "pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbeccen,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,"
                  "pl_orblper,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,"
                  "pl_ratror,pl_ratrorerr1,pl_ratrorerr2,"
                  "st_teff,st_tefferr1,st_tefferr2,st_met,st_meterr1,st_meterr2,"
                  "st_logg,st_loggerr1,st_loggerr2,st_mass,st_rad,st_raderr1",
        "from": "pscomppars",
        "where": "tran_flag = 1",
        "format": "csv"
    }
    if target:
        uri_ipac_query["where"] += f" and pl_name = '{target}'"
    return tap_query(uri_ipac_base, uri_ipac_query)

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

    return parser.parse_args()

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
    search_result = lk.search_targetpixelfile(prior['hostname'], mission='TESS')
    print(search_result)

    if len(search_result) == 0:
        raise Exception(f"No data for: {prior['pl_name']}")

    # Load prior data if it exists, else scrape from NEA
    if len(prior) == 0:
        nea_df = nea_scrape(prior['pl_name'])
        prior = {key: nea_df[key][0] for key in nea_df.keys()}
        prior["pl_name"] = planetname.replace('-', '')

    # find sectors
    sectors = [int(sector.split(' ')[-1]) for sector in search_result.mission]
    if args.sector:
        if args.sector not in sectors:
            raise(Exception(f"No data for sector: {args.sector}"))

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
        if args.sector:
            if sector != args.sector:
                continue

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
        'per':[tpars['per'] - 0.01, tpars['per'] + 0.01]
    }

    if args.ars:
        mybounds['ars'] = [tpars['ars'] * 0.01, tpars['ars'] * 3]

    print('performing global fit...')
    airmass = np.zeros(len(time[tmask]))
    myfit = lc_fitter(time[tmask] + 2457000.0, flux[tmask], phot_std / flux[tmask], airmass, tpars, mybounds, verbose=True)

    # create plots
    fig, ax = myfit.plot_bestfit(title=f"{prior['pl_name']} Global Fit")
    ax[0].set_ylim([np.percentile(flux, 1) * 0.99, np.percentile(flux, 99) * 1.01])
    plt.savefig(os.path.join(planetdir, planetname + "_global_fit.png"))
    plt.close()

    myfit.plot_triangle()
    plt.savefig(os.path.join(planetdir, planetname + "_global_triangle.png"))
    plt.close()

    # update priors from best fit
    prior['pl_orbper'] = myfit.parameters['per']
    prior['pl_orbpererr1'] = myfit.errors['per']
    prior['pl_orbpererr2'] = -myfit.errors['per']
    prior['pl_tranmid'] = myfit.parameters['tmid']
    prior['pl_tranmiderr1'] = myfit.errors['tmid']
    prior['pl_tranmiderr2'] = -myfit.errors['tmid']

    if not args.ars:
        prior['pl_orbincl'] = myfit.parameters['inc']
        prior['pl_orbinclerr1'] = myfit.errors['inc']
        prior['pl_orbinclerr2'] = -myfit.errors['inc']

    # save prior to disk
    with open(os.path.join(planetdir, planetname + "_prior.json"), 'w', encoding='utf8') as json_file:
        json.dump(prior, json_file, indent=4)

    # save results to state vector
    sv['global'] = {
        'pars': myfit.parameters,
        'errs': myfit.errors,
        'time': myfit.time,
        'data': myfit.data,
        'data_err': myfit.dataerr,
    }

    # remove transit signal from timeseries
    residuals = np.copy(flux)
    residuals[tmask] /= myfit.model
    mask = (residuals < 1.01) & (residuals > 0.99)

    if args.tls:
        model = transitleastsquares(time[mask], residuals[mask])
        results = model.power(R_star=0.5, M_star=0.5)

        plt.figure()
        ax = plt.gca()
        ax.axvline(results.period, alpha=0.4, lw=3)
        plt.xlim(np.min(results.periods), np.max(results.periods))
        for n in range(2, 10):
            ax.axvline(n * results.period, alpha=0.4, lw=1, linestyle="dashed")
            ax.axvline(results.period / n, alpha=0.4, lw=1, linestyle="dashed")
        plt.ylabel(r'SDE')
        plt.xlabel('Period (days)')
        plt.plot(results.periods, results.power, color='black', lw=0.5)
        plt.xlim(0, max(results.periods))
        plt.savefig(os.path.join(planetdir, planetname + "_periodogram.png"))
        plt.close()

    # save global fit data if above certain SNR
    rprs2 = myfit.parameters['rprs'] ** 2
    rprs2err = myfit.parameters['rprs'] * 2 * myfit.errors['rprs']
    snr = rprs2 / rprs2err

    if snr < 1:
        with open("notes.txt", 'w') as f:
            f.write(f"Skipping individual light curve fits b.c SNR = {snr:.2f}")
        pickle.dump(sv, open(os.path.join(planetdir, planetname + "_data.pkl"), "wb"))
        raise(Exception(f"Skipping individual light curve fits b.c SNR = {snr:.2f}"))

    period = myfit.parameters['per']
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

    for e in events:
        tmask = (tphase > e - (2 * pdur)) & (tphase < e + (2 * pdur))

        if tmask.sum() == 0:
            continue

        tpars['tmid'] = tmid + e * period
        tpars['per'] = period
        tpars['ars'] = ars

        mybounds = {
            'rprs': [0, 3 * tpars['rprs']],
            'tmid': [tpars['tmid'] - 0.25 * pdur * tpars['per'], tpars['tmid'] + 0.25 * pdur * tpars['per']],
            'inc': [tpars['inc'] - 5, min(90, tpars['inc'] + 5)],
        }

        if np.sum(np.diff(np.sort(time[tmask])) * 24 * 60 > 35):
            continue

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

        if (rprs2 - 1 * rprs2err <= 0):
            print(f"Skipping transit at phase: {e} b.c. depth = {rprs2:.2e} +- {rprs2err:.2e}")
            continue

        if myfit.transit[0] < 1 or myfit.transit[-1] < 1:
            continue

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

        csv_lk = OutputFiles(myfit, prior, infoDict, planetdir)
        csv_lk.aavso_csv(airmass, u0, u1, u2, u3, tmidstr)
        csv_lk.aavso(airmass, u0, u1, u2, u3, tmidstr)

    pickle.dump(sv, open(os.path.join(planetdir, planetname + "_data.pkl"), "wb"))
