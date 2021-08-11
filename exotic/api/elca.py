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
# ########################################################################### #
# Exoplanet light curve analysis
#
# Authors: Kyle Pearson, Marlee Smith, Rob Zellem
# Supplemental Code: Gael Roudier and Jason Eastman
# ########################################################################### #

import copy
from numba import njit
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
# from scipy.ndimage import gaussian_filter as norm_kde
from scipy.signal import savgol_filter
# from scipy.stats import gaussian_kde
from ultranest import ReactiveNestedSampler


def tldlc(z, rprs, g1=0, g2=0, g3=0, g4=0, nint=int(2**3)):
    """
    G. ROUDIER: Light curve model
    """
    ldlc = np.zeros(z.size)
    xin = z.copy() - rprs
    xin[xin < 0e0] = 0e0
    xout = z.copy() + rprs
    xout[xout > 1e0] = 1e0
    select = xin > 1e0
    if True in select:
        ldlc[select] = 1e0
    xint = np.linspace(1e0, 0e0, nint)
    znot = z.copy()[~select]
    xinnot = np.arccos(xin[~select])
    xoutnot = np.arccos(xout[~select])
    xrs = np.array([xint]).T*(xinnot - xoutnot) + xoutnot
    xrs = np.cos(xrs)
    diffxrs = np.diff(xrs, axis=0)
    extxrs = np.zeros((xrs.shape[0]+1, xrs.shape[1]))
    extxrs[1:-1, :] = xrs[1:, :] - diffxrs/2.
    extxrs[0, :] = xrs[0, :] - diffxrs[0]/2.
    extxrs[-1, :] = xrs[-1, :] + diffxrs[-1]/2.
    occulted = vecoccs(znot, extxrs, rprs)
    diffocc = np.diff(occulted, axis=0)
    si = vecistar(xrs, g1, g2, g3, g4)
    drop = np.sum(diffocc*si, axis=0)
    inldlc = 1. - drop
    ldlc[~select] = np.array(inldlc)
    return ldlc


def vecistar(xrs, g1, g2, g3, g4):
    """
    G. ROUDIER: Stellar surface extinction model
    """
    ldnorm = (-g1/10e0 - g2/6e0 - 3e0*g3/14e0 - g4/4e0 + 5e-1)*2e0*np.pi
    select = xrs < 1e0
    mu = np.zeros(xrs.shape)
    mu[select] = (1e0 - xrs[select]**2)**(1e0/4e0)
    s1 = g1*(1e0 - mu)
    s2 = g2*(1e0 - mu**2)
    s3 = g3*(1e0 - mu**3)
    s4 = g4*(1e0 - mu**4)
    outld = (1e0 - (s1+s2+s3+s4))/ldnorm
    return outld


def vecoccs(z, xrs, rprs):
    """
    G. ROUDIER: Stellar surface occulation model
    """
    out = np.zeros(xrs.shape)
    vecxrs = xrs.copy()
    selx = vecxrs > 0e0
    veczsel = np.array([z.copy()]*xrs.shape[0])
    veczsel[veczsel < 0e0] = 0e0
    select1 = (vecxrs <= rprs - veczsel) & selx
    select2 = (vecxrs >= rprs + veczsel) & selx
    select = (~select1) & (~select2) & selx
    zzero = veczsel == 0e0
    if True in select1 & zzero:
        out[select1 & zzero] = np.pi*(np.square(vecxrs[select1 & zzero]))
        pass
    if True in select2 & zzero:
        out[select2 & zzero] = np.pi*(rprs**2)
    if True in select & zzero:
        out[select & zzero] = np.pi*(rprs**2)
    if True in select1 & ~zzero:
        out[select1 & ~zzero] = np.pi*(np.square(vecxrs[select1 & ~zzero]))
        pass
    if True in select2:
        out[select2 & ~zzero] = np.pi*(rprs**2)
    if True in select & ~zzero:
        redxrs = vecxrs[select & ~zzero]
        redz = veczsel[select & ~zzero]
        s1 = (np.square(redz) + np.square(redxrs) - rprs**2)/(2e0*redz*redxrs)
        s1[s1 > 1e0] = 1e0
        s2 = (np.square(redz) + rprs**2 - np.square(redxrs))/(2e0*redz*(rprs+0.0001))
        s2[s2 > 1e0] = 1e0
        fs3 = -redz + redxrs + rprs
        ss3 = redz + redxrs - rprs
        ts3 = redz - redxrs + rprs
        os3 = redz + redxrs + rprs
        s3 = fs3*ss3*ts3*os3
        zselect = s3 < 0e0
        if True in zselect:
            s3[zselect] = 0e0
        out[select & ~zzero] = (np.square(redxrs)*np.arccos(s1) +
                                (rprs**2)*np.arccos(s2) - (5e-1)*np.sqrt(s3))
        pass
    return out


@njit(cache=True)
def time2z(time, ipct, tknot, sma, orbperiod, ecc, tperi=None, epsilon=1e-5):
    """
    G. ROUDIER: Time samples in [Days] to separation in [R*]
    """
    tperii = 0
    if tperi is not None:
        ft0 = (tperi - tknot) % orbperiod
        ft0 /= orbperiod
        if ft0 > 0.5:
            ft0 += -1e0
        M0 = 2e0*np.pi*ft0
        E0 = solveme(M0, ecc, epsilon)
        realf = np.sqrt(1e0 - ecc)*np.cos(E0/2e0)
        imagf = np.sqrt(1e0 + ecc)*np.sin(E0/2e0)
        w = np.angle(np.complex(realf, imagf))
        if abs(ft0) < epsilon:
            w = np.pi/2e0
            tperii = tknot
            pass
        pass
    else:
        w = np.pi/2e0
        tperii = tknot
        pass

    ft = (time - tperii) % orbperiod
    ft /= orbperiod
    sft = np.copy(ft)
    sft[(sft > 0.5)] += -1e0
    M = 2e0*np.pi*ft
    E = solveme(M, ecc, epsilon)
    realf = np.sqrt(1. - ecc)*np.cos(E/2e0)
    imagf = np.sqrt(1. + ecc)*np.sin(E/2e0)
    st = np.vstack((realf, imagf))
    f = np.array([2e0 * np.angle(np.complex(st[0][i], st[1][i])) for i in range(len(imagf))])
    r = sma*(1e0 - ecc**2)/(1e0 + ecc*np.cos(f))
    z = r*np.sqrt(1e0**2 - (np.sin(w+f)**2)*(np.sin(ipct*np.pi/180e0))**2)
    z[sft < 0] *= -1e0
    return z, sft


@njit(cache=True)
def solveme(M, e, eps):
    """
    G. ROUDIER: Newton Raphson solver for true anomaly
    M is a numpy array
    """
    E = M.copy()  # numba only allows copy() with no arguments
    for i in np.arange(M.shape[0]):
        while abs(E[i] - e*np.sin(E[i]) - M[i]) > eps:
            num = E[i] - e*np.sin(E[i]) - M[i]
            den = 1. - e*np.cos(E[i])
            E[i] = E[i] - num/den
            pass
        pass
    return E


def transit(times, values):
    sep, phase = time2z(times, values['inc'], values['tmid'], values['ars'], values['per'], values['ecc'])
    model = tldlc(abs(sep), values['rprs'], values['u0'], values['u1'], values['u2'], values['u3'])
    return model


def getPhase(curTime, pPeriod, tMid):
    phase = (curTime - tMid) / pPeriod
    return phase - int(np.nanmin(phase))


@njit(fastmath=True)
def mc_a1(m_a2, sig_a2, transit, airmass, data, n=10000):
    a2 = np.random.normal(m_a2, sig_a2, n)
    model = transit * np.exp(np.mean(a2) * airmass)
    detrend = data / model
    a1 = np.random.normal(np.mean(detrend), np.std(detrend), n)
    return np.mean(a1), np.std(a1)


# average data into bins of dt from start to finish
def time_bin(time, flux, dt):
    bins = int(np.floor((max(time) - min(time))/dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time)+i*dt)) & (time < (min(time)+(i+1)*dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
    zmask = (bflux == 0) | (btime == 0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask]


# Function that bins an array
def binner(arr, n, err=''):
    if len(err) == 0:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant', constant_values=np.NaN).reshape(-1, n)
        arr = np.nanmean(ecks, axis=1)
        return arr
    else:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant', constant_values=np.NaN).reshape(-1, n)
        why = np.pad(err.astype(float), (0, ((n - err.size % n) % n)), mode='constant', constant_values=np.NaN).reshape(-1, n)
        weights = 1./(why**2.)
        # Calculate the weighted average
        arr = np.nansum(ecks * weights, axis=1) / np.nansum(weights, axis=1)
        err = np.array([np.sqrt(1. / np.nansum(1. / (np.array(i) ** 2.))) for i in why])
        return arr, err


class lc_fitter(object):

    def __init__(self, time, data, dataerr, airmass, prior, bounds, log_dir=None, mode='ns', verbose=True):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.airmass = airmass
        self.prior = prior
        self.bounds = bounds
        self.max_ncalls = 2e5
        self.log_dir = log_dir
        self.verbose = verbose
        if mode == "lm":
            self.fit_LM()
        elif mode == "ns":
            self.fit_nested()

    def fit_LM(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])

        def lc2min(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= np.exp(self.prior['a2']*self.airmass)
            detrend = self.data / model
            model *= np.median(detrend)
            return ((self.data-model)/self.dataerr)**2

        try:
            res = least_squares(lc2min, x0=[self.prior[k] for k in freekeys],
                                bounds=[boundarray[:, 0], boundarray[:, 1]], jac='3-point', loss='linear')
        except Exception as e:
            print(f"{e} \nbounded light curve fitting failed...check priors "
                  "(e.g. estimated mid-transit time + orbital period)")

            for i, k in enumerate(freekeys):
                if not boundarray[i, 0] < self.prior[k] < boundarray[i, 1]:
                    print(f"bound: [{boundarray[i, 0]}, {boundarray[i, 1]}] prior: {self.prior[k]}")

            print("removing bounds and trying again...")
            res = least_squares(lc2min, x0=[self.prior[k] for k in freekeys], method='lm', jac='3-point', loss='linear')

        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}

        for i, k in enumerate(freekeys):
            self.parameters[k] = res.x[i]
            self.errors[k] = 0

        self.create_fit_variables()

    def create_fit_variables(self):
        self.phase = (self.time - self.parameters['tmid'] + 0.25 * self.parameters['per']) / self.parameters['per'] % 1 - 0.25
        self.transit = transit(self.time, self.parameters)
        self.parameters['a1'], self.errors['a1'] = mc_a1(self.parameters['a2'], self.errors['a2'],
                                                         self.transit, self.airmass, self.data)
        self.airmass_model = self.parameters['a1'] * np.exp(self.parameters['a2'] * self.airmass)
        self.model = self.transit * self.airmass_model
        self.detrended = self.data / self.airmass_model
        self.detrendederr = self.dataerr / self.airmass_model
        self.residuals = self.data - self.model
        self.chi2 = np.sum(self.residuals ** 2 / self.dataerr ** 2)
        self.bic = len(self.bounds) * np.log(len(self.time)) - 2 * np.log(self.chi2)

        # compare fit chi2 to smoothed data chi2
        dt = np.diff(np.sort(self.time)).mean()
        si = np.argsort(self.time)
        try:
            self.sdata = savgol_filter(self.data[si], 1+2*int(0.5/24/dt), 2)
        except:
            self.sdata = np.ones(len(self.time))

        schi2 = np.sum((self.data[si] - self.sdata)**2/self.dataerr[si]**2)
        self.quality = schi2/self.chi2

        # measured duration
        tdur = (self.transit < 1).sum() * np.median(np.diff(np.sort(self.time)))

        # test for partial transit
        newtime = np.linspace(self.parameters['tmid']-0.2, self.parameters['tmid']+0.2, 10000)
        newtran = transit(newtime, self.parameters)
        masktran = newtran < 1
        newdur = np.diff(newtime).mean()*masktran.sum()

        self.duration_measured = tdur
        self.duration_expected = newdur

    def fit_nested(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)

        def loglike(pars):
            # chi-squared
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= np.exp(self.prior['a2'] * self.airmass)
            detrend = self.data / model  # used to estimate a1
            model *= np.median(detrend)
            return -0.5 * np.sum(((self.data - model) / self.dataerr) ** 2)

        @njit(fastmath=True)
        def prior_transform(upars):
            # transform unit cube to prior volume
            return boundarray[:, 0] + bounddiff * upars

        # TODO vectorize loglike
        test = ReactiveNestedSampler(freekeys, loglike, prior_transform, log_dir=self.log_dir)
        self.results = test.run(max_ncalls=int(self.max_ncalls))

        # plots
        test.plot()

        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # final model
        self.create_fit_variables()

    def plot_bestfit(self, nbins=10, phase=True, title=""):
        # import pdb; pdb.set_trace()

        f, (ax_lc, ax_res) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})

        if phase:
            ax_res.set_xlabel('Phase')

            ecks = self.phase

        else:
            ax_res.set_xlabel('Time [day]')

            ecks = self.time

        # clip plot to get rid of white space
        ax_res.set_xlim([min(ecks), max(ecks)])
        ax_lc.set_xlim([min(ecks), max(ecks)])

        # making borders and tick labels black
        ax_lc.spines['bottom'].set_color('black')
        ax_lc.spines['top'].set_color('black')
        ax_lc.spines['right'].set_color('black')
        ax_lc.spines['left'].set_color('black')
        ax_lc.tick_params(axis='x', colors='black')
        ax_lc.tick_params(axis='y', colors='black')

        ax_res.spines['bottom'].set_color('black')
        ax_res.spines['top'].set_color('black')
        ax_res.spines['right'].set_color('black')
        ax_res.spines['left'].set_color('black')
        ax_res.tick_params(axis='x', colors='black')
        ax_res.tick_params(axis='y', colors='black')

        # residual plot
        ax_res.errorbar(ecks, self.residuals / np.median(self.data), yerr=self.detrendederr, color='gray',
                        marker='o', markersize=5, linestyle='None', mec='None', alpha=0.75)
        ax_res.plot(ecks, np.zeros(len(ecks)), 'r-', lw=2, alpha=1, zorder=100)
        ax_res.set_ylabel('Residuals')
        ax_res.set_ylim([-3 * np.nanstd(self.residuals / np.median(self.data)),
                         3 * np.nanstd(self.residuals / np.median(self.data))])

        correctedSTD = np.std(self.residuals / np.median(self.data))
        ax_lc.errorbar(ecks, self.detrended, yerr=self.detrendederr, ls='none',
                       marker='o', color='gray', markersize=5, mec='None', alpha=0.75)
        ax_lc.plot(ecks, self.transit, 'r', zorder=1000, lw=2)

        ax_lc.set_ylabel('Relative Flux')
        ax_lc.get_xaxis().set_visible(False)

        ax_res.errorbar(binner(ecks, len(self.residuals) // 10),
                        binner(self.residuals / np.median(self.data), len(self.residuals) // 10),
                        yerr=
                        binner(self.residuals / np.median(self.data), len(self.residuals) // 10, self.detrendederr)[
                            1],
                        fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)
        ax_lc.errorbar(binner(ecks, len(ecks) // 10),
                       binner(self.detrended, len(self.detrended) // 10),
                       yerr=
                       binner(self.residuals / np.median(self.data), len(self.residuals) // 10, self.detrendederr)[
                           1],
                       fmt='s', ms=5, mfc='b', mec='None', ecolor='b', zorder=10)

        # remove vertical whitespace
        f.subplots_adjust(hspace=0)

        return f, (ax_lc, ax_res)


if __name__ == "__main__":

    prior = {
        'rprs': 0.03,                                   # Rp/Rs
        'ars': 14.25,                                   # a/Rs
        'per': 3.336817,                                # Period [day]
        'inc': 87.5,                                    # Inclination [deg]
        'u0': 1.8, 'u1': -3.3, 'u2': 3.9, 'u3': -1.5,   # limb darkening (nonlinear)
        'ecc': 0,                                       # Eccentricity
        'omega': 0,                                     # Arg of periastron
        'tmid': 0.75,                                   # Time of mid transit [day],
        # 'a1': 50,                                     # Airmass coefficients
        'a2': 0.25
    }

    time = np.linspace(0.65, 0.85, 150)  # [day]

    # simulate extinction from airmass
    stime = time-time[0]
    alt = 90 * np.cos(4*stime-np.pi/6)
    airmass = 1./np.cos(np.deg2rad(90-alt))

    # GENERATE NOISY DATA
    data = transit(time, prior)*prior['a1']*np.exp(prior['a2']*airmass)
    data += np.random.normal(0, prior['a1']*250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time))

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, 0.1],
        'tmid': [min(time), max(time)],
        'ars': [13, 15],
        # 'a1': [25, 75],
        'a2': [0, 0.3]
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    fig, axs = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()
