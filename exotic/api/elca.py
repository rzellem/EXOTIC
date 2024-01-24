# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #
# Exoplanet light curve analysis
#
# Fit an exoplanet transit model to time series data.
# ########################################################################### #
from astropy.time import Time
import copy
from itertools import cycle
import matplotlib.pyplot as plt
import numpy as np
from pylightcurve.models.exoplanet_lc import transit as pytransit
from scipy import spatial
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
try:
    from ultranest import ReactiveNestedSampler
except ImportError:
    import dynesty
    import dynesty.plotting
    from dynesty.utils import resample_equal
    from scipy.stats import gaussian_kde

try:
    from plotting import corner
except ImportError:
    from .plotting import corner


def weightedflux(flux, gw, nearest):
    return np.sum(flux[nearest] * gw, axis=-1)


def gaussian_weights(X, w=1, neighbors=50, feature_scale=1000):
    Xm = (X - np.median(X, 0)) * w
    kdtree = spatial.cKDTree(Xm * feature_scale)
    nearest = np.zeros((X.shape[0], neighbors))
    gw = np.zeros((X.shape[0], neighbors), dtype=float)
    for point in range(X.shape[0]):
        ind = kdtree.query(kdtree.data[point], neighbors + 1)[1][1:]
        dX = Xm[ind] - Xm[point]
        Xstd = np.std(dX, 0)
        gX = np.exp(-dX ** 2 / (2 * Xstd ** 2))
        gwX = np.product(gX, 1)
        gw[point, :] = gwX / gwX.sum()
        nearest[point, :] = ind
    gw[np.isnan(gw)] = 0.01
    return gw, nearest.astype(int)


def transit(times, values):
    model = pytransit([values['u0'], values['u1'], values['u2'], values['u3']],
                      values['rprs'], values['per'], values['ars'],
                      values['ecc'], values['inc'], values['omega'],
                      values['tmid'], times, method='claret', precision=3)
    return model


def get_phase(times, per, tmid):
    return (times - tmid + 0.25 * per) / per % 1 - 0.25


def mc_a1(m_a2, sig_a2, transit, airmass, data, n=10000):
    a2 = np.random.normal(m_a2, sig_a2, n)
    model = transit * np.exp(np.repeat(np.expand_dims(a2, 0), airmass.shape[0], 0).T * airmass)
    # change to linear model
    #model = transit * np.repeat(np.expand_dims(a2, 0), airmass.shape[0], 0).T * airmass
    detrend = data / model
    return np.mean(np.median(detrend, 0)), np.std(np.median(detrend, 0))


def round_to_2(*args):
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if np.floor(y) >= 1.:
        roundval = 2
    else:
        try:
            roundval = -int(np.floor(np.log10(abs(y)))) + 1
        except:
            roundval = 1
    return round(x, roundval)


# average data into bins of dt from start to finish
def time_bin(time, flux, dt=1. / (60 * 24)):
    bins = int(np.floor((max(time) - min(time)) / dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    bstds = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time) + i * dt)) & (time < (min(time) + (i + 1) * dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
            bstds[i] = np.nanstd(flux[mask]) / (mask.sum() ** 0.5)
    zmask = (bflux == 0) | (btime == 0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask], bstds[~zmask]


# Function that bins an array
def binner(arr, n, err=''):
    if len(err) == 0:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant',
                      constant_values=np.NaN).reshape(-1, n)
        arr = np.nanmean(ecks, axis=1)
        return arr
    else:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant',
                      constant_values=np.NaN).reshape(-1, n)
        why = np.pad(err.astype(float), (0, ((n - err.size % n) % n)), mode='constant', constant_values=np.NaN).reshape(
            -1, n)
        weights = 1. / (why ** 2.)
        # Calculate the weighted average
        arr = np.nansum(ecks * weights, axis=1) / np.nansum(weights, axis=1)
        err = np.array([np.sqrt(1. / np.nansum(1. / (np.array(i) ** 2.))) for i in why])
        return arr, err


class lc_fitter(object):

    def __init__(self, time, data, dataerr, airmass, prior, bounds, neighbors=200, mode='ns', verbose=True):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.airmass = airmass
        self.prior = prior
        self.bounds = bounds
        self.max_ncalls = 2e5
        self.verbose = verbose
        self.mode = mode
        self.neighbors = neighbors
        self.results = None
        if self.mode == "lm":
            self.fit_LM()
        elif self.mode == "ns":
            self.fit_nested()

    def fit_LM(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])

        # trim data around predicted transit/eclipse time
        if np.ndim(self.airmass) == 2:
            print(f'Computing nearest neighbors and gaussian weights for {len(self.time)} npts...')
            self.gw, self.nearest = gaussian_weights(self.airmass, neighbors=self.neighbors)

        def lc2min_nneighbor(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            lightcurve = transit(self.time, self.prior)
            detrended = self.data / lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve * wf
            return ((self.data - model) / self.dataerr) ** 2

        def lc2min_airmass(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= self.prior['a1'] * np.exp(self.prior['a2'] * self.airmass)
            return ((self.data - model) / self.dataerr) ** 2

        try:
            if np.ndim(self.airmass) == 2:
                res = least_squares(lc2min_nneighbor, x0=[self.prior[k] for k in freekeys],
                                    bounds=[boundarray[:, 0], boundarray[:, 1]], jac='3-point', loss='linear')
            else:
                res = least_squares(lc2min_airmass, x0=[self.prior[k] for k in freekeys],
                                    bounds=[boundarray[:, 0], boundarray[:, 1]], jac='3-point', loss='linear')
        except Exception as e:
            print(f"{e} \nbounded light curve fitting failed...check priors "
                  "(e.g. estimated mid-transit time + orbital period)")

            for i, k in enumerate(freekeys):
                if not boundarray[i, 0] < self.prior[k] < boundarray[i, 1]:
                    print(f"bound: [{boundarray[i, 0]}, {boundarray[i, 1]}] prior: {self.prior[k]}")

            print("removing bounds and trying again...")

            if np.ndim(self.airmass) == 2:
                res = least_squares(lc2min_nneighbor, x0=[self.prior[k] for k in freekeys],
                                    method='lm', jac='3-point', loss='linear')
            else:
                res = least_squares(lc2min_airmass, x0=[self.prior[k] for k in freekeys],
                                    method='lm', jac='3-point', loss='linear')

        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}

        for i, k in enumerate(freekeys):
            self.parameters[k] = res.x[i]
            self.errors[k] = 0

        self.create_fit_variables()

    def create_fit_variables(self):
        self.phase = get_phase(self.time, self.parameters['per'], self.parameters['tmid'])
        self.transit = transit(self.time, self.parameters)
        self.time_upsample = np.linspace(min(self.time), max(self.time), 1000)
        self.transit_upsample = transit(self.time_upsample, self.parameters)
        self.phase_upsample = get_phase(self.time_upsample, self.parameters['per'], self.parameters['tmid'])
        if self.mode == "ns":
            self.parameters['a1'], self.errors['a1'] = mc_a1(self.parameters.get('a2', 0), self.errors.get('a2', 1e-6),
                                                             self.transit, self.airmass, self.data)
        if np.ndim(self.airmass) == 2:
            detrended = self.data / self.transit
            self.wf = weightedflux(detrended, self.gw, self.nearest)
            self.model = self.transit * self.wf
            self.detrended = self.data / self.wf
            self.detrendederr = self.dataerr / self.wf
        else:
            self.airmass_model = self.parameters['a1'] * np.exp(self.parameters.get('a2', 0) * self.airmass)
            self.model = self.transit * self.airmass_model
            self.detrended = self.data / self.airmass_model
            self.detrendederr = self.dataerr / self.airmass_model

        self.residuals = self.data - self.model
        self.res_stdev = np.std(self.residuals)/np.median(self.data)
        self.chi2 = np.sum(self.residuals ** 2 / self.dataerr ** 2)
        self.bic = len(self.bounds) * np.log(len(self.time)) - 2 * np.log(self.chi2)

        # compare fit chi2 to smoothed data chi2
        dt = np.diff(np.sort(self.time)).mean()
        si = np.argsort(self.time)
        try:
            self.sdata = savgol_filter(self.data[si], 1 + 2 * int(0.5 / 24 / dt), 2)
        except:
            self.sdata = np.ones(len(self.time))

        schi2 = np.sum((self.data[si] - self.sdata) ** 2 / self.dataerr[si] ** 2)
        self.quality = schi2 / self.chi2

        # measured duration
        tdur = (self.transit < 1).sum() * np.median(np.diff(np.sort(self.time)))

        # test for partial transit
        newtime = np.linspace(self.parameters['tmid'] - 0.2, self.parameters['tmid'] + 0.2, 10000)
        newtran = transit(newtime, self.parameters)
        masktran = newtran < 1
        newdur = np.diff(newtime).mean() * masktran.sum()

        self.duration_measured = tdur
        self.duration_expected = newdur

    def fit_nested(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)

        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        def loglike(pars):
            # chi-squared
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= np.exp(self.prior['a2'] * self.airmass)
            # linear model
            #model *= self.prior['a2'] * (self.time-min(self.time))
            # check for all zeros
            #if np.all(model == 0):
            #    return -1e10
            detrend = self.data / model  # used to estimate a1
            model *= np.median(detrend)
            return -0.5 * np.sum(((self.data - model) / self.dataerr) ** 2)

        def prior_transform(upars):
            # transform unit cube to prior volume
            return boundarray[:, 0] + bounddiff * upars

        try:
            self.ns_type = 'ultranest'
            test = ReactiveNestedSampler(freekeys, loglike, prior_transform)

            noop = lambda *args, **kwargs: None
            if self.verbose is True:
                self.results = test.run(max_ncalls=int(self.max_ncalls))
            else:
                self.results = test.run(max_ncalls=int(self.max_ncalls), show_status=False, viz_callback=noop)

            for i, key in enumerate(freekeys):
                self.parameters[key] = self.results['maximum_likelihood']['point'][i]
                self.errors[key] = self.results['posterior']['stdev'][i]
                self.quantiles[key] = [
                    self.results['posterior']['errlo'][i],
                    self.results['posterior']['errup'][i]]
        except NameError:
            self.ns_type = 'dynesty'
            dsampler = dynesty.DynamicNestedSampler(loglike, prior_transform, ndim=len(freekeys),
                                                    bound='multi', sample='unif')
            dsampler.run_nested(maxcall=int(1e5), dlogz_init=0.05,
                                maxbatch=10, nlive_batch=100, print_progressbool=self.verbose)
            self.results = dsampler.results

            tests = [copy.deepcopy(self.prior) for i in range(5)]

            # Derive kernel density estimate for best fit
            weights = np.exp(self.results.logwt - self.results.logz[-1])
            samples = self.results['samples']
            logvol = self.results['logvol']
            wt_kde = gaussian_kde(resample_equal(-logvol, weights))  # KDE
            logvol_grid = np.linspace(logvol[0], logvol[-1], 1000)  # resample
            wt_grid = wt_kde.pdf(-logvol_grid)  # evaluate KDE PDF
            self.weights = np.interp(-logvol, -logvol_grid, wt_grid)  # interpolate

            # errors + final values
            mean, cov = dynesty.utils.mean_and_cov(self.results.samples, weights)
            mean2, cov2 = dynesty.utils.mean_and_cov(self.results.samples, self.weights)
            for i in range(len(freekeys)):
                self.errors[freekeys[i]] = cov[i, i] ** 0.5
                tests[0][freekeys[i]] = mean[i]
                tests[1][freekeys[i]] = mean2[i]

                counts, bins = np.histogram(samples[:, i], bins=100, weights=weights)
                mi = np.argmax(counts)
                tests[4][freekeys[i]] = bins[mi] + 0.5 * np.mean(np.diff(bins))

                # finds median and +- 2sigma, will vary from mode if non-gaussian
                self.quantiles[freekeys[i]] = dynesty.utils.quantile(self.results.samples[:, i], [0.025, 0.5, 0.975],
                                                                     weights=weights)
                tests[2][freekeys[i]] = self.quantiles[freekeys[i]][1]

            # find minimum near weighted mean
            mask = (samples[:, 0] < self.parameters[freekeys[0]] + 2 * self.errors[freekeys[0]]) & (
                    samples[:, 0] > self.parameters[freekeys[0]] - 2 * self.errors[freekeys[0]])
            bi = np.argmin(self.weights[mask])

            for i in range(len(freekeys)):
                tests[3][freekeys[i]] = samples[mask][bi, i]
                # tests[4][freekeys[i]] = np.average(samples[mask][:, i], weights=self.weights[mask], axis=0)

            # find best fit from chi2 minimization
            chis = []
            for i in range(len(tests)):
                lightcurve = transit(self.time, tests[i])
                tests[i]['a1'] = mc_a1(tests[i].get('a2', 0), self.errors.get('a2', 1e-6),
                                       lightcurve, self.airmass, self.data)[0]
                airmass = tests[i]['a1'] * np.exp(tests[i].get('a2', 0) * self.airmass)
                residuals = self.data - (lightcurve * airmass)
                chis.append(np.sum(residuals ** 2))

            mi = np.argmin(chis)
            self.parameters = copy.deepcopy(tests[mi])

        # final model
        self.create_fit_variables()

    def plot_bestfit(self, title="", bin_dt=30. / (60 * 24), zoom=False, phase=True):
        f = plt.figure(figsize=(9, 6))
        f.subplots_adjust(top=0.92, bottom=0.09, left=0.14, right=0.98, hspace=0)
        ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
        ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True, ls='--')

        rprs2 = self.parameters['rprs'] ** 2
        rprs2err = 2 * self.parameters['rprs'] * self.errors['rprs']
        lclabel1 = r"$R^{2}_{p}/R^{2}_{s}$ = %s $\pm$ %s" % (
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )

        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" % (
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid', 0))),
            str(round_to_2(self.errors.get('tmid', 0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2

        if zoom:
            axs[0].set_ylim([1 - 1.25 * self.parameters['rprs'] ** 2, 1 + 0.5 * self.parameters['rprs'] ** 2])
        else:
            if phase:
                axs[0].errorbar(self.phase, self.detrended, yerr=np.std(self.residuals) / np.median(self.data),
                                ls='none', marker='.', color='black', zorder=1, alpha=0.2)
            else:
                axs[0].errorbar(self.time, self.detrended, yerr=np.std(self.residuals) / np.median(self.data),
                                ls='none', marker='.', color='black', zorder=1, alpha=0.2)

        if phase:
            si = np.argsort(self.phase)
            bt2, br2, _ = time_bin(self.phase[si] * self.parameters['per'],
                                   self.residuals[si] / np.median(self.data) * 1e2, bin_dt)
            axs[1].plot(self.phase, self.residuals / np.median(self.data) * 1e2, 'k.', alpha=0.2,
                        label=r'$\sigma$ = {:.2f} %'.format(np.std(self.residuals / np.median(self.data) * 1e2)))
            axs[1].plot(bt2 / self.parameters['per'], br2, 'bs', alpha=1, zorder=2)
            axs[1].set_xlim([min(self.phase), max(self.phase)])
            axs[1].set_xlabel("Phase", fontsize=14)

            si = np.argsort(self.phase)
            bt2, bf2, bs = time_bin(self.phase[si] * self.parameters['per'], self.detrended[si], bin_dt)
            axs[0].errorbar(bt2 / self.parameters['per'], bf2, yerr=bs, alpha=1, zorder=2, color='blue', ls='none',
                            marker='s')
            # axs[0].plot(self.phase[si], self.transit[si], 'r-', zorder=3, label=lclabel)
            sii = np.argsort(self.phase_upsample)
            axs[0].plot(self.phase_upsample[sii], self.transit_upsample[sii], 'r-', zorder=3, label=lclabel)
            axs[0].set_xlim([min(self.phase), max(self.phase)])
            axs[0].set_xlabel("Phase ", fontsize=14)
        else:
            bt, br, _ = time_bin(self.time, self.residuals / np.median(self.data) * 1e2, bin_dt)
            axs[1].plot(self.time, self.residuals / np.median(self.data) * 1e2, 'k.', alpha=0.2,
                        label=r'$\sigma$ = {:.2f} %'.format(np.std(self.residuals / np.median(self.data) * 1e2)))
            axs[1].plot(bt, br, 'bs', alpha=1, zorder=2, label=r'$\sigma$ = {:.2f} %'.format(np.std(br)))
            axs[1].set_xlim([min(self.time), max(self.time)])
            axs[1].set_xlabel("Time [day]", fontsize=14)

            bt, bf, bs = time_bin(self.time, self.detrended, bin_dt)
            si = np.argsort(self.time)
            sii = np.argsort(self.time_upsample)
            axs[0].errorbar(bt, bf, yerr=bs, alpha=1, zorder=2, color='blue', ls='none', marker='s')
            axs[0].plot(self.time_upsample[sii], self.transit_upsample[sii], 'r-', zorder=3, label=lclabel)
            axs[0].set_xlim([min(self.time), max(self.time)])
            axs[0].set_xlabel("Time [day]", fontsize=14)

        axs[0].get_xaxis().set_visible(False)
        axs[1].legend(loc='best')
        axs[0].legend(loc='best')
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True, ls='--', axis='y')
        return f, axs

    def plot_triangle(self):
        if self.ns_type == 'ultranest':
            ranges = []
            mask1 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
            mask2 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
            mask3 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
            titles = []
            labels = []
            flabels = {
                'rprs': r'R$_{p}$/R$_{s}$',
                'per': r'Period [day]',
                'tmid': r'T$_{mid}$',
                'ars': r'a/R$_{s}$',
                'inc': r'Inc. [deg]',
                'u1': r'u$_1$',
                'fpfs': r'F$_{p}$/F$_{s}$',
                'omega': r'$\omega$ [deg]',
                'mplanet': r'M$_{p}$ [M$_{\oplus}$]',
                'mstar': r'M$_{s}$ [M$_{\odot}$]',
                'ecc': r'$e$',
                'c0': r'$c_0$',
                'c1': r'$c_1$',
                'c2': r'$c_2$',
                'c3': r'$c_3$',
                'c4': r'$c_4$',
                'a0': r'$a_0$',
                'a1': r'$a_1$',
                'a2': r'$a_2$'
            }
            for i, key in enumerate(self.quantiles):
                labels.append(flabels.get(key, key))
                titles.append(f"{self.parameters[key]:.5f} +- {self.errors[key]:.5f}")
                ranges.append([
                    self.parameters[key] - 5 * self.errors[key],
                    self.parameters[key] + 5 * self.errors[key]
                ])

                if key == 'a2' or key == 'a1':
                    continue

                mask3 = mask3 & \
                    (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - 3 * self.errors[key])) & \
                    (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 3 * self.errors[key]))

                mask1 = mask1 & \
                    (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - self.errors[key])) & \
                    (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + self.errors[key]))

                mask2 = mask2 & \
                    (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - 2 * self.errors[key])) & \
                    (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 2 * self.errors[key]))

            chi2 = self.results['weighted_samples']['logl'] * -2
            fig = corner(self.results['weighted_samples']['points'],
                         labels=labels,
                         bins=int(np.sqrt(self.results['samples'].shape[0])),
                         range=ranges,
                         # quantiles=(0.1, 0.84),
                         plot_contours=True,
                         levels=[np.percentile(chi2[mask1], 95), np.percentile(chi2[mask2], 95),
                                 np.percentile(chi2[mask3], 95)],
                         plot_density=False,
                         titles=titles,
                         data_kwargs={
                             'c': chi2,
                             'vmin': np.percentile(chi2[mask3], 1),
                             'vmax': np.percentile(chi2[mask3], 95),
                             'cmap': 'viridis'
                         },
                         label_kwargs={
                             'labelpad': 15,
                         },
                         hist_kwargs={
                             'color': 'black',
                         }
                         )
        else:
            fig, axs = dynesty.plotting.cornerplot(self.results, labels=list(self.bounds.keys()),
                                                   quantiles_2d=[0.4, 0.85],
                                                   smooth=0.015, show_titles=True, use_math_text=True, title_fmt='.2e',
                                                   hist2d_kwargs={ 'fill_contours': False})
            dynesty.plotting.cornerpoints(self.results, labels=list(self.bounds.keys()),
                                          fig=[fig, axs[1:, :-1]], plot_kwargs={'alpha': 0.1, 'zorder': 1, })
        return fig

# simultaneously fit multiple data sets with global and local parameters
class glc_fitter(lc_fitter):
    # needed for lc_fitter
    ns_type = 'ultranest'

    def __init__(self, input_data, global_bounds, local_bounds, individual_fit=False, stdev_cutoff=0.03, verbose=False):
        # keys for input_data: time, flux, ferr, airmass, priors all numpy arrays
        self.lc_data = copy.deepcopy(input_data)
        self.global_bounds = global_bounds
        self.local_bounds = local_bounds
        self.individual_fit = individual_fit
        self.stdev_cutoff = stdev_cutoff
        self.verbose = verbose
        self.results = None
        self.fit_nested()

    def fit_nested(self):

        # create bound arrays for generating samples
        nobs = len(self.lc_data)
        gfreekeys = list(self.global_bounds.keys())

        # if isinstance(self.local_bounds, dict):
        #     lfreekeys = list(self.local_bounds.keys())
        #     boundarray = np.vstack([ [self.global_bounds[k] for k in gfreekeys], [self.local_bounds[k] for k in lfreekeys]*nobs ])
        # else:
        #     # if list type
        lfreekeys = []
        boundarray = [self.global_bounds[k] for k in gfreekeys]
        for i in range(nobs):
            lfreekeys.append(list(self.local_bounds[i].keys()))
            boundarray.extend([self.local_bounds[i][k] for k in lfreekeys[-1]])
        boundarray = np.array(boundarray)

        # fit individual light curves to constrain priors
        if self.individual_fit:
            for i in range(nobs):

                print(f"Fitting individual light curve {i+1}/{nobs}")
                try:
                    mybounds = dict(**self.local_bounds[i], **self.global_bounds)
                except:
                    mybounds = {}
                    for k in self.local_bounds[i]:
                        mybounds[k] = self.local_bounds[i][k]
                    for k in self.global_bounds:
                        mybounds[k] = self.global_bounds[k]
                if 'per' in mybounds: del(mybounds['per'])
                if 'inc' in mybounds and 'rprs' in mybounds: del(mybounds['inc'])
                if 'tmid' in mybounds:
                    # find the closet mid transit time to the last observation
                    phase = (self.lc_data[i]['time'][-1] - self.lc_data[i]['priors']['tmid']) / self.lc_data[i]['priors']['per']
                    nepochs = np.round(phase)
                    newtmid = self.lc_data[i]['priors']['tmid'] + nepochs * self.lc_data[i]['priors']['per']
                    err = np.diff(mybounds['tmid'])[0]/2.
                    mybounds['tmid'] = [newtmid - err, newtmid + err ]

                # fit individual light curve
                myfit = lc_fitter(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'],
                    self.lc_data[i]['ferr'],
                    self.lc_data[i]['airmass'],
                    self.lc_data[i]['priors'],
                    mybounds
                )

                # check stdev_cutoff and residuals
                if myfit.res_stdev > self.stdev_cutoff:
                    print(f"WARNING: Stdev of residuals is large! {myfit.res_stdev:.3f} > {self.stdev_cutoff:.3f}")
                    #raise ValueError(f"Stdev of residuals is too large, please remove id: {self.lc_data[i]['name']}

                # copy data over for individual fits
                self.lc_data[i]['individual'] = myfit.parameters.copy()
                self.lc_data[i]['individual_err'] = myfit.errors.copy()
                self.lc_data[i]['res_stdev'] = myfit.res_stdev
                self.lc_data[i]['quality'] = myfit.quality

                ti = sum([len(self.local_bounds[k]) for k in range(i)])
                # update local priors
                for j, key in enumerate(self.local_bounds[i].keys()):

                    boundarray[j+ti+len(gfreekeys),0] = myfit.parameters[key] - 5*myfit.errors[key]
                    boundarray[j+ti+len(gfreekeys),1] = myfit.parameters[key] + 5*myfit.errors[key]

                    if key == 'rprs':
                        boundarray[j+ti+len(gfreekeys),0] = max(0,myfit.parameters[key] - 5*myfit.errors[key])

                # print name and stdev of residuals
                mint = np.min(self.lc_data[i]['time'])
                maxt = np.max(self.lc_data[i]['time'])
                try:
                    print(f"{self.lc_data[i]['name']} & {Time(mint,format='jd').isot} & {Time(maxt,format='jd').isot} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}")
                except:
                    print(f"{self.lc_data[i]['name']} & {mint} & {maxt} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}")

                del(myfit)

        # transform unit cube to prior volume
        bounddiff = np.diff(boundarray,1).reshape(-1)
        def prior_transform(upars):
            return (boundarray[:,0] + bounddiff*upars)

        def loglike(pars):
            chi2 = 0

            # for each light curve
            for i in range(nobs):

                # global keys
                for j, key in enumerate(gfreekeys):
                    self.lc_data[i]['priors'][key] = pars[j]

                # local keys
                ti = sum([len(self.local_bounds[k]) for k in range(i)])
                for j, key in enumerate(lfreekeys[i]):
                    self.lc_data[i]['priors'][key] = pars[j+ti+len(gfreekeys)]

                # compute model
                model = transit(self.lc_data[i]['time'], self.lc_data[i]['priors'])
                model *= np.exp(self.lc_data[i]['priors']['a2']*self.lc_data[i]['airmass'])
                detrend = self.lc_data[i]['flux']/model
                model *= np.mean(detrend)

                # add to chi2
                chi2 += np.sum( ((self.lc_data[i]['flux']-model)/self.lc_data[i]['ferr'])**2 )

            # maximization metric for nested sampling
            return -0.5*chi2

        freekeys = []+gfreekeys
        for n in range(nobs):
            for k in lfreekeys[n]:
                #clean_name = self.lc_data[n].get('name', n).replace(' ','_').replace('(','').replace(')','').replace('[','').replace(']','').replace('-','_').split('-')[0]
                freekeys.append(f"local_{k}_{n}")

        noop = lambda *args, **kwargs: None
        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e6, show_status=True)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=1e6, show_status=False, viz_callback=noop)

        self.quantiles = {}
        self.errors = {}
        self.parameters = self.lc_data[0]['priors'].copy()

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            #self.errors[key] = self.results['posterior']['median'][i]

            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # create an average Rp/Rs if it is not in global keys
        # check if 'rprs' is in lfreekeys
        rprs_in_local = False
        for i in range(nobs):
            if 'rprs' in lfreekeys[i]:
                rprs_in_local = True
                break

        if rprs_in_local:
            local_rprs = [] # used for creating an average value
            local_rprs_err = []

        # loop over observations
        for n in range(nobs):
            self.lc_data[n]['errors'] = {}
            
            # copy global parameters
            self.lc_data[n]['priors'] = copy.deepcopy(self.parameters)

            # loop over local keys and save best fit values
            for k in lfreekeys[n]:

                # create key to get results
                pkey = f"local_{k}_{n}"
                
                # overwrite priors with best fit value
                self.lc_data[n]['priors'][k] = self.parameters[pkey]
                self.lc_data[n]['errors'][k] = self.errors[pkey]

                # update key for final bestfit plot if needed
                if k == 'rprs':
                    local_rprs.append(self.lc_data[n]['priors'][k])
                    local_rprs_err.append(self.lc_data[n]['errors'][k])

            # solve for a1
            model = transit(self.lc_data[n]['time'], self.lc_data[n]['priors'])
            airmass = np.exp(self.lc_data[n]['airmass']*self.lc_data[n]['priors']['a2'])
            detrend = self.lc_data[n]['flux']/(model*airmass)
            self.lc_data[n]['priors']['a1'] = np.mean(detrend)
            self.lc_data[n]['residuals'] = self.lc_data[n]['flux'] - model*airmass*self.lc_data[n]['priors']['a1']
            self.lc_data[n]['detrend'] = self.lc_data[n]['flux']/(airmass*self.lc_data[n]['priors']['a1'])

            # phase
            self.lc_data[n]['phase'] = get_phase(self.lc_data[n]['time'], self.lc_data[n]['priors']['per'], self.lc_data[n]['priors']['tmid'])
            self.lc_data[n]['time_upsample'] = np.linspace(min(self.lc_data[n]['time']), max(self.lc_data[n]['time']), 1000)
            self.lc_data[n]['phase_upsample'] = get_phase(self.lc_data[n]['time_upsample'], self.lc_data[n]['priors']['per'], self.lc_data[n]['priors']['tmid'])
            self.lc_data[n]['transit_upsample'] = transit(self.lc_data[n]['time_upsample'], self.lc_data[n]['priors'])

        # create an average value from all the local fits
        if rprs_in_local:
            self.parameters['rprs'] = np.mean(local_rprs)
            self.errors['rprs'] = np.std(local_rprs)

        #import pdb; pdb.set_trace()

    def plot_bestfits(self):
        nrows = len(self.lc_data)//4+1
        # make sure there isn't an extra row
        if len(self.lc_data)%4 == 0:
            nrows -= 1

        fig,ax = plt.subplots(nrows, 4, figsize=(5+5*nrows, 5*nrows))

        # turn off all axes
        for i in range(nrows*4):
            ri = int(i/4)
            ci = i%4
            if ax.ndim == 1:
                ax[i].axis('off')
            else:
                ax[ri,ci].axis('off')

        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])

        # plot observations
        for i in range(len(self.lc_data)):
            ri = int(i/4)
            ci = i%4
            ncolor = next(colors)
            nmarker = next(markers)

            model = transit(self.lc_data[i]['time'], self.lc_data[i]['priors'])
            airmass = np.exp(self.lc_data[i]['airmass']*self.lc_data[i]['priors']['a2'])
            detrend = self.lc_data[i]['flux']/(model*airmass)

            if ax.ndim == 1:
                ax[i].axis('on')
                ax[i].errorbar(self.lc_data[i]['time'], self.lc_data[i]['flux']/airmass/detrend.mean(), yerr=self.lc_data[i]['ferr']/airmass/detrend.mean(), 
                                ls='none', marker=nmarker, color=ncolor, alpha=0.5, zorder=1)
                
                ax[i].plot(self.lc_data[i]['time_upsample'], self.lc_data[i]['transit_upsample'], 'r-', zorder=2)
                ax[i].set_xlabel("Time [BJD]", fontsize=14)
                ax[i].set_ylabel("Relative Flux", fontsize=14)
                ax[i].set_title(f"{self.lc_data[i].get('name','')}", fontsize=16)
            else:
                ax[ri,ci].axis('on')
                ax[ri,ci].errorbar(self.lc_data[i]['time'], self.lc_data[i]['flux']/airmass/detrend.mean(), yerr=self.lc_data[i]['ferr']/airmass/detrend.mean(), 
                                   ls='none', marker=nmarker, color=ncolor, alpha=0.5, zorder=1)
                ax[ri,ci].plot(self.lc_data[i]['time_upsample'], self.lc_data[i]['transit_upsample'], 'r-', zorder=2)
                ax[ri,ci].set_xlabel("Time[BJD]", fontsize=14)
                ax[ri,ci].set_ylabel("Relative Flux", fontsize=14)
                ax[ri,ci].set_title(f"{self.lc_data[i].get('name','')}", fontsize=16)

        plt.tight_layout()
        return fig

    def plot_bestfit(self, title="", bin_dt=30./(60*24), alpha=0.05, ylim_sigma=5, phase_limits='median', show_legend=True, limit_legend=False, show_individual_fits=False):
        """
        Plot the best fit model and residuals

        Parameters
        ----------
        title : str
            Title for the plot

        bin_dt : float
            Bin size for plotting the residuals

        alpha : float
            Alpha value for plotting the data

        ylim_sigma : float
            Number of sigma to plot the residuals

        phase_limits : str
            'median' or 'all' to set the phase limits

        show_legend : bool
            Show the legend

        limit_legend : bool
            Limit the legend to 3 entries
        """
        f = plt.figure(figsize=(15,12))
        f.subplots_adjust(top=0.92,bottom=0.09,left=0.1,right=0.98, hspace=0)
        ax_lc = plt.subplot2grid((4,5), (0,0), colspan=5,rowspan=3)
        ax_res = plt.subplot2grid((4,5), (3,0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title, fontsize=18)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True,ls='--')

        try:
            rprs2 = self.parameters['rprs']**2
            rprs2err = 2*self.parameters['rprs']*self.errors['rprs']
        except:
            rprs2 = self.lc_data[0]['priors']['rprs']**2
            rprs2err = 2*self.lc_data[0]['priors']['rprs']*self.lc_data[0]['errors']['rprs']

        lclabel1 = r"$R^{2}_{p}/R^{2}_{s}$ = %s $\pm$ %s" %(
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )
        
        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
            str(round_to_2(self.errors.get('tmid',0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2
        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])

        alldata = {
            'time': [],
            'flux': [],
            'detrend': [],
            'ferr': [],
            'residuals': [],
        }

        for n in range(len(self.lc_data)):
            ncolor = next(colors)
            nmarker = next(markers)
            alldata['time'].extend(self.lc_data[n]['time'].tolist())
            alldata['detrend'].extend(self.lc_data[n]['detrend'].tolist())
            alldata['flux'].extend(self.lc_data[n]['flux'].tolist())
            alldata['ferr'].extend(self.lc_data[n]['ferr'].tolist())
            alldata['residuals'].extend(self.lc_data[n]['residuals'].tolist())
            
            phase = get_phase(self.lc_data[n]['time'], self.parameters['per'], self.lc_data[n]['priors']['tmid'])
            si = np.argsort(phase)
            #bt2, br2, _ = time_bin(phase[si]*self.parameters['per'], self.lc_data[n]['residuals'][si]/np.median(self.lc_data[n]['flux'])*1e2, bin_dt)

            # plot data
            axs[0].errorbar(phase, self.lc_data[n]['detrend'], yerr=np.std(self.lc_data[n]['residuals'])/np.median(self.lc_data[n]['flux']), 
                            ls='none', marker=nmarker, color=ncolor, zorder=1, alpha=alpha)

            # plot residuals
            axs[1].plot(phase, self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])*1e2, color=ncolor, marker=nmarker, ls='none',
                         alpha=0.2)

            # plot binned data
            bt2, bf2, bs = time_bin(phase[si]*self.lc_data[n]['priors']['per'], self.lc_data[n]['detrend'][si], bin_dt)

            if limit_legend:
                axs[0].errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker)
            else:
                axs[0].errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker,
                                label=r'{}: {:.2f} %'.format(self.lc_data[n].get('name',''),np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])*1e2)))

            # replace min and max for upsampled lc model
            minp = min(minp, min(phase))
            maxp = max(maxp, max(phase))
            min_std = min(min_std, np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])))

            # plot individual best fit models
            if show_individual_fits:
                axs[0].plot(self.lc_data[n]['phase_upsample'], self.lc_data[n]['transit_upsample'], color=ncolor, zorder=3, alpha=0.5)

        # create binned plot for all the data
        for k in alldata.keys():
            alldata[k] = np.array(alldata[k])
            
        phase = get_phase(alldata['time'], self.parameters['per'], self.lc_data[n]['priors']['tmid'])
        si = np.argsort(phase)
        bt, br, _ = time_bin(phase[si]*self.parameters['per'], alldata['residuals'][si]/np.median(alldata['flux']), 2*bin_dt)
        bt, bf, bs = time_bin(phase[si]*self.parameters['per'], alldata['detrend'][si], 2*bin_dt)

        axs[0].errorbar(bt/self.parameters['per'],bf,yerr=bs,alpha=1,zorder=2,color='white',ls='none',marker='o',ms=15,
                        markeredgecolor='black',
                        ecolor='black',
                        label=r'Binned Data: {:.2f} %'.format(np.std(br)*1e2))

        axs[1].plot(bt/self.parameters['per'],br*1e2,color='white',ls='none',marker='o',ms=11,markeredgecolor='black')

        # best fit model
        self.time_upsample = np.linspace(minp*self.parameters['per']+self.parameters['tmid'], 
                                         maxp*self.parameters['per']+self.parameters['tmid'], 10000)
        self.transit_upsample = transit(self.time_upsample, self.parameters)
        self.phase_upsample = get_phase(self.time_upsample, self.parameters['per'], self.parameters['tmid'])
        sii = np.argsort(self.phase_upsample)
        axs[0].plot(self.phase_upsample[sii], self.transit_upsample[sii], 'r-', zorder=3, label=lclabel, lw=3)

        # set up axes limits
        axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[0].set_xlabel("Phase ", fontsize=14)
        axs[0].set_ylim([1-self.parameters['rprs']**2-ylim_sigma*min_std, 1+ylim_sigma*min_std])
        axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[1].set_xlabel("Phase", fontsize=14)
        axs[1].set_ylim([-5*min_std*1e2, 5*min_std*1e2])

        # compute average min and max for all the data
        mins = []; maxs = []
        for n in range(len(self.lc_data)):
            mins.append(min(self.lc_data[n]['phase']))
            maxs.append(max(self.lc_data[n]['phase']))

        # set up phase limits
        if isinstance(phase_limits, str):
            if phase_limits == "minmax":
                axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
                axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            elif phase_limits == "median":
                axs[0].set_xlim([np.median(mins), np.median(maxs)])
                axs[1].set_xlim([np.median(mins), np.median(maxs)])
            else:
                axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
                axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        elif isinstance(phase_limits, list):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        elif isinstance(phase_limits, tuple):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        else:
            axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])

        axs[0].get_xaxis().set_visible(False)
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True,ls='--',axis='y')
    
        if show_legend:
            axs[0].legend(loc='best',ncol=len(self.lc_data)//7+1)
    
        return f,axs

    def plot_stack(self, title="", bin_dt=30./(60*24), dy=0.02):
        f, ax = plt.subplots(1,figsize=(9,12))
        
        ax.set_title(title)
        ax.set_ylabel("Relative Flux", fontsize=14)
        ax.grid(True,ls='--')

        rprs2 = self.parameters['rprs']**2
        rprs2err = 2*self.parameters['rprs']*self.errors['rprs']
        lclabel1 = r"$R^{2}_{p}/R^{2}_{s}$ = %s $\pm$ %s" %(
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )
        
        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
            str(round_to_2(self.errors.get('tmid',0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2
        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])
        for n in range(len(self.lc_data)):
            ncolor = next(colors)
            nmarker = next(markers)

            phase = get_phase(self.lc_data[n]['time'], self.parameters['per'], self.lc_data[n]['priors']['tmid'])
            si = np.argsort(phase)
            bt2, br2, _ = time_bin(phase[si]*self.parameters['per'], self.lc_data[n]['residuals'][si]/np.median(self.lc_data[n]['flux'])*1e2, bin_dt)
            
            # plot data
            ax.errorbar(phase, self.lc_data[n]['detrend']-n*dy, yerr=np.std(self.lc_data[n]['residuals'])/np.median(self.lc_data[n]['flux']), 
                            ls='none', marker=nmarker, color=ncolor, zorder=1, alpha=0.25)
        
            # plot binned data
            bt2, bf2, bs = time_bin(phase[si]*self.lc_data[n]['priors']['per'], self.lc_data[n]['detrend'][si]-n*dy, bin_dt)
            ax.errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker)

            # replace min and max for upsampled lc model
            minp = min(minp, min(phase))
            maxp = max(maxp, max(phase))
            min_std = min(min_std, np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])))

            # best fit model
            self.time_upsample = np.linspace(minp*self.parameters['per']+self.parameters['tmid'], 
                                            maxp*self.parameters['per']+self.parameters['tmid'], 10000)
            self.transit_upsample = transit(self.time_upsample, self.parameters)
            self.phase_upsample = get_phase(self.time_upsample, self.parameters['per'], self.parameters['tmid'])
            sii = np.argsort(self.phase_upsample)
            ax.plot(self.phase_upsample[sii], self.transit_upsample[sii]-n*dy, ls='-', color=ncolor, zorder=3, label=self.lc_data[n].get('name',''))

        ax.set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        ax.set_xlabel("Phase ", fontsize=14)
        ax.set_ylim([1-self.parameters['rprs']**2-5*min_std-n*dy, 1+5*min_std])
        ax.get_xaxis().set_visible(False)
        ax.legend(loc='best')
        return f,ax


if __name__ == "__main__":

    prior = {
        'rprs': 0.02,  # Rp/Rs
        'ars': 14.25,  # a/Rs
        'per': 3.33,  # Period [day]
        'inc': 88.5,  # Inclination [deg]
        'u0': 0, 'u1': 0, 'u2': 0, 'u3': 0,  # limb darkening (nonlinear)
        'ecc': 0.5,  # Eccentricity
        'omega': 120,  # Arg of periastron
        'tmid': 0.75,  # Time of mid transit [day],
        'a1': 0,  # Airmass coefficients
        'a2': 0.,  # trend = a1 * np.exp(a2 * airmass) or a1 + a2 * airmass

        'teff': 5000,
        'tefferr': 50,
        'met': 0,
        'meterr': 0,
        'logg': 3.89,
        'loggerr': 0.01
    }

    # example generating LD coefficients
    from pylightcurve import exotethys

    u0, u1, u2, u3 = exotethys(prior['logg'], prior['teff'], prior['met'], 'TESS', method='claret',
                               stellar_model='phoenix')

    prior['u0'], prior['u1'], prior['u2'], prior['u3'] = u0, u1, u2, u3

    time = np.linspace(0.7, 0.8, 1000)  # [day]

    # simulate extinction from airmass
    stime = time - time[0]
    alt = 90 * np.cos(4 * stime - np.pi / 6)
    airmass = 1./np.cos(np.deg2rad(90-alt))
    #airmass = np.zeros(time.shape[0])

    # GENERATE NOISY DATA
    data = transit(time, prior) * prior['a1'] * np.exp(prior['a2'] * airmass)
    data += np.random.normal(0, prior['a1'] * 250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time)) + np.random.normal(300e-6, 50e-6, len(time))

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, 0.1],
        'tmid': [prior['tmid'] - 0.01, prior['tmid'] + 0.01],
        'ars': [13, 15],
        'a2': [-0.1,0.1] # uncomment if you want to fit for airmass
        # never list 'a1' in bounds, it is perfectly correlated to exp(a2*airmass)
        # and is solved for during the fit
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    fig, axs = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()

    fig = myfit.plot_triangle()
    plt.tight_layout()
    plt.show()
