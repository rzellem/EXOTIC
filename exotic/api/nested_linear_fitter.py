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
# Various functions for fitting a linear model to data, including nested
# sampling, linear least squares. Includes residual plotting, posteriors and 
# a periodogram analysis.
# 
# ########################################################################### #
import numpy as np
from itertools import cycle
import statsmodels.api as sm
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle
try:
    from ultranest import ReactiveNestedSampler
except ImportError:
    print("Warning: ultranest not installed. Nested sampling will not work.")

try:
    from plotting import corner
except ImportError:
    from .plotting import corner

class linear_fitter(object):

    def __init__(self, data, dataerr, bounds=None, prior=None, labels=None, verbose=True):
        """
        Fit a linear model to data using nested sampling.

        Parameters
        ----------
        data : array
            Data to fit.
        dataerr : array
            Error on data.
        bounds : dict, optional
            Bounds on parameters. Dictionary of tuples/list
        prior : dict, optional
            Prior on parameters for slope and intercept. Dictionary of tuples/list
        """
        self.data = data
        self.dataerr = dataerr
        self.bounds = bounds
        self.labels = np.array(labels)
        self.prior = prior.copy()  # dict {'m':(0.1,0.5), 'b':(0,1)}
        self.verbose = verbose
        if bounds is None:
            # use +- 3 sigma prior as bounds
            self.bounds = {
                'm': [prior['m'][0] - 3 * prior['m'][1], prior['m'][0] + 3 * prior['m'][1]],
                'b': [prior['b'][0] - 3 * prior['b'][1], prior['b'][0] + 3 * prior['b'][1]]
            }
        self.results = None
        self.fit_nested()

    def fit_nested(self):
        """ Fit a linear model to data using nested sampling. """

        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)
        self.epochs = np.round((self.data - np.mean(self.bounds['b'])) / np.mean(self.bounds['m']))

        def loglike(pars):
            # chi-squared
            model = pars[0] * self.epochs + pars[1]
            return -0.5 * np.sum(((self.data - model) / self.dataerr) ** 2)

        def prior_transform(upars):
            # transform unit cube to prior volume
            return (boundarray[:, 0] + bounddiff * upars)

        # estimate slope and intercept
        noop = lambda *args, **kwargs: None
        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=4e5,
                                                                                         min_num_live_points=420,
                                                                                         show_status=True)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=4e5,
                                                                                         min_num_live_points=420,
                                                                                         show_status=False,
                                                                                         viz_callback=noop)
        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = {}

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # final model
        self.model = self.epochs * self.parameters['m'] + self.parameters['b']
        self.residuals = self.data - self.model

    def plot_oc(self, savefile=None, ylim='none', show_2sigma=False, prior_name="Prior"):
        """ Plot the data in the form of residuals vs. time

        Parameters
        ----------
        savefile : str, optional
            Save the figure to a file.
        ylim : str, optional
            Set the y-axis limits. Default is 'none'. Can be prior, average, or none.
        show_2sigma : bool, optional
            Show a fill between using the 2 sigma limits. Default is False (aka 1 sigma)
        """

        # set up the figure        
        fig, ax = plt.subplots(1, figsize=(9, 6))

        # check if labels are not None
        if self.labels.size:
            # find unique set of labels
            ulabels = np.unique(self.labels)
            # set up a color/marker cycle
            markers = cycle(['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
            colors = cycle(['black', 'blue', 'green', 'orange', 'purple', 'grey', 'magenta', 'cyan', 'lime'])

            # plot each label separately
            for i, ulabel in enumerate(ulabels):
                # find where the label matches
                mask = self.labels == ulabel
                # plot the data/residuals
                ax.errorbar(self.epochs[mask], self.residuals[mask] * 24 * 60, yerr=self.dataerr[mask] * 24 * 60,
                            ls='none', marker=next(markers), color=next(colors), label=ulabel)
        else:
            # plot the data/residuals
            ax.errorbar(self.epochs, self.residuals * 24 * 60, yerr=self.dataerr * 24 * 60,
                        ls='none', marker='o', color='black')

        ylower = (self.residuals.mean() - 3 * np.std(self.residuals)) * 24 * 60
        yupper = (self.residuals.mean() + 3 * np.std(self.residuals)) * 24 * 60

        # upsample data
        epochs = (np.linspace(self.data.min() - 7, self.data.max() + 7, 1000) -
                  self.parameters['b']) / self.parameters['m']

        # set the y-axis limits
        depoch = self.epochs.max() - self.epochs.min()
        ax.set_xlim([self.epochs.min() - depoch * 0.01, self.epochs.max() + depoch * 0.01])

        # best fit solution
        model = epochs * self.parameters['m'] + self.parameters['b']

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['m'], self.errors['m'], size=10000)
        mc_b = np.random.normal(self.parameters['b'], self.errors['b'], size=10000)
        mc_model = np.expand_dims(epochs, -1) * mc_m + mc_b

        # create a fill between area for uncertainty of new ephemeris
        diff = mc_model.T - model

        if show_2sigma:
            ax.fill_between(epochs, np.percentile(diff, 2, axis=0) * 24 * 60, np.percentile(diff, 98, axis=0) * 24 * 60,
                            alpha=0.2, color='k', label=r'Uncertainty ($\pm$ 2$\sigma$)')
        else:
            # show 1 sigma
            ax.fill_between(epochs, np.percentile(diff, 36, axis=0) * 24 * 60,
                            np.percentile(diff, 64, axis=0) * 24 * 60, alpha=0.2, color='k',
                            label=r'Uncertainty ($\pm$ 1$\sigma$)')

        # duplicate axis and plot days since mid-transit
        ax2 = ax.twiny()
        ax2.set_xlabel(f"Time [BJD - {self.parameters['b']:.1f}]", fontsize=14)
        ax2.set_xlim(ax.get_xlim())
        xticks = ax.get_xticks()
        dt = np.round(xticks * self.parameters['m'], 1)
        # ax2.set_xticks(dt)
        ax2.set_xticklabels(dt)

        if ylim == 'diff':
            ax.set_ylim([min(np.percentile(diff, 1, axis=0) * 24 * 60),
                         max(np.percentile(diff, 99, axis=0) * 24 * 60)])

        # overlay the prior ephemeris
        if self.prior is not None:
            # create fill between area for uncertainty of old/prior ephemeris
            epochs_p = (np.linspace(self.data.min() - 7, self.data.max() + 7, 1000) -
                        self.prior['b'][0]) / self.prior['m'][0]
            prior = epochs_p * self.prior['m'][0] + self.prior['b'][0]
            mc_m_p = np.random.normal(self.prior['m'][0], self.prior['m'][1], size=10000)
            mc_b_p = np.random.normal(self.prior['b'][0], self.prior['b'][1], size=10000)
            mc_model_p = np.expand_dims(epochs_p, -1) * mc_m_p + mc_b_p
            diff_p = mc_model_p.T - model

            # plot an invisible line so the 2nd axes are happy
            ax2.plot(epochs, (model - prior) * 24 * 60, ls='--', color='r', alpha=0)

            # why is this so small!?!?!? consistent to within machine precision?
            # ax.plot(epochs, (model-prior)*24*60, ls='--', color='r')

            if show_2sigma:
                ax.fill_between(epochs, np.percentile(diff_p, 2, axis=0) * 24 * 60,
                                np.percentile(diff_p, 98, axis=0) * 24 * 60, alpha=0.1, color='r',
                                label=rf'{prior_name} ($\pm$ 2$\sigma$)')
            else:
                # show ~1 sigma
                ax.fill_between(epochs, np.percentile(diff_p, 36, axis=0) * 24 * 60,
                                np.percentile(diff_p, 64, axis=0) * 24 * 60, alpha=0.1, color='r',
                                label=rf'{prior_name} ($\pm$ 1$\sigma$)')

            if ylim == 'prior':
                ax.set_ylim([min(np.percentile(diff_p, 1, axis=0) * 24 * 60),
                             max(np.percentile(diff_p, 99, axis=0) * 24 * 60)])
            elif ylim == 'average':
                ax.set_ylim([0.5 * (min(np.percentile(diff, 1, axis=0) * 24 * 60) +
                                    min(np.percentile(diff_p, 1, axis=0) * 24 * 60)),
                             0.5 * (max(np.percentile(diff, 99, axis=0) * 24 * 60) +
                                    max(np.percentile(diff_p, 99, axis=0) * 24 * 60))])

        ax.axhline(0, color='black', alpha=0.5, ls='--',
                   label="Period: {:.7f}+-{:.7f} days\nT_mid: {:.7f}+-{:.7f} BJD".format(self.parameters['m'],
                                                                                         self.errors['m'],
                                                                                         self.parameters['b'],
                                                                                         self.errors['b']))

        # TODO sig figs
        # lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
        #    str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
        #    str(round_to_2(self.errors.get('tmid',0)))
        # )

        ax.legend(loc='best')
        ax.set_xlabel("Epoch [number]", fontsize=14)
        ax.set_ylabel("Residuals [min]", fontsize=14)
        ax.grid(True, ls='--')
        return fig, ax

    def plot_triangle(self):
        """ Create a posterior triangle plot of the results. """
        ranges = []
        mask1 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        mask2 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        mask3 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        titles = []
        labels = []
        flabels = {
            'm': 'Period [day]',
            'b': 'T_mid [JD]',
        }
        for i, key in enumerate(self.quantiles):
            labels.append(flabels.get(key, key))
            titles.append(f"{self.parameters[key]:.7f} +-\n {self.errors[key]:.7f}")

            # set the axes limits for the plots
            ranges.append([
                self.parameters[key] - 5 * self.errors[key],
                self.parameters[key] + 5 * self.errors[key]
            ])

            if key == 'a2' or key == 'a1':
                continue

            # create masks for contouring on sigma bounds
            mask3 = (mask3 &
                     (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - 3 * self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 3 * self.errors[key])))

            mask1 = (mask1 &
                     (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + self.errors[key])))

            mask2 = (mask2 &
                     (self.results['weighted_samples']['points'][:, i] > ( self.parameters[key] - 2 * self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 2 * self.errors[key])))

        chi2 = self.results['weighted_samples']['logl'] * -2
        fig = corner(self.results['weighted_samples']['points'],
                     labels=labels,
                     bins=int(np.sqrt(self.results['samples'].shape[0])),
                     range=ranges,
                     figsize=(10, 10),
                     # quantiles=(0.1, 0.84),
                     plot_contours=True,
                     levels=[np.percentile(chi2[mask1], 95), np.percentile(chi2[mask2], 95),
                             np.percentile(chi2[mask3], 95)],
                     plot_density=False,
                     titles=titles,
                     data_kwargs={
                         'c': chi2,  # color code by chi2
                         'vmin': np.percentile(chi2[mask3], 1),
                         'vmax': np.percentile(chi2[mask3], 95),
                         'cmap': 'viridis',
                     },
                     label_kwargs={'labelpad': 50, },
                     hist_kwargs={'color': 'black', }
                     )
        return fig

    def plot_periodogram(self, minper=0, maxper=0, maxper2=50):
        """ Search the residuals for periodic signals. """

        ########################################
        # create basis vectors for Tn = T0 + n*P
        basis = np.ones((2, len(self.epochs)))
        basis[1] = self.epochs
        res_linear = sm.WLS(self.data, basis.T, weights=1.0 / self.dataerr ** 2).fit()
        coeffs_linear = res_linear.params  # retrieve the slope and intercept
        y_bestfit_linear = np.dot(basis.T, coeffs_linear)  # reconstruct signal
        residuals = self.data - y_bestfit_linear
        ########################################

        # compute a period range based on nyquist frequency
        si = np.argsort(self.epochs)

        if minper == 0:
            minper = max(3, 2 * np.diff(self.epochs[si]).min())
        if maxper == 0:
            maxper = (np.max(self.epochs) - np.min(self.epochs)) * 3.

        # recompute on new grid
        ls = LombScargle(self.epochs, residuals, dy=self.dataerr)
        freq, power = ls.autopower(maximum_frequency=1. / minper, minimum_frequency=1. / maxper, nyquist_factor=2)

        # Phase fold data at max peak
        mi = np.argmax(power)
        per = 1. / freq[mi]
        newphase = self.epochs / per % 1
        self.periods = 1. / freq
        self.power = power

        ########################################
        # create basis vectors for Tn = T0 + n*P + Asin(wn) + Bcos(wn)
        basis = np.ones((4, len(self.epochs)))
        basis[1] = self.epochs
        basis[2] = np.sin(2 * np.pi * self.epochs / per)
        basis[3] = np.cos(2 * np.pi * self.epochs / per)

        # perform the weighted least squares regression
        res_first_order = sm.WLS(self.data, basis.T, weights=1.0 / self.dataerr ** 2).fit()
        coeffs_first_order = res_first_order.params
        y_bestfit_first_order = np.dot(basis.T, coeffs_first_order)  # reconstruct signal
        residuals_first_order = self.data - y_bestfit_first_order

        # reconstruct signal with first two terms, i.e. linear solution (T0 + n*P)
        y_bestfit_linear_first_order = np.dot(basis[:2].T, coeffs_first_order[:2])
        residuals_linear = self.data - y_bestfit_linear_first_order
        ########################################

        ########################################
        # subtract first order solution from data and recompute periodogram
        maxper = maxper2

        # recompute on new grid
        ls2 = LombScargle(self.epochs, residuals_linear, dy=self.dataerr)
        # freq2,power2 = ls.autopower(maximum_frequency=1./(1.01*per),
        #                             minimum_frequency=1./maxper, nyquist_factor=2)
        freq2, power2 = ls.autopower(maximum_frequency=1. / (1 + minper),
                                     minimum_frequency=1. / maxper, nyquist_factor=2)

        # find max period
        mi2 = np.argmax(power2)
        per2 = 1. / freq2[mi2]

        # create basis vectors for second order solution
        basis = np.ones((6, len(self.epochs)))
        basis[1] = self.epochs
        basis[2] = np.sin(2 * np.pi * self.epochs / per)
        basis[3] = np.cos(2 * np.pi * self.epochs / per)
        basis[4] = np.sin(4 * np.pi * self.epochs / per2)
        basis[5] = np.cos(4 * np.pi * self.epochs / per2)

        # perform the weighted least squares regression
        res_second_order = sm.WLS(self.data, basis.T, weights=1.0 / self.dataerr ** 2).fit()
        coeffs_second_order = res_second_order.params
        y_bestfit_second_order = np.dot(basis.T, coeffs_second_order)  # reconstruct signal
        ########################################

        # find the best bic
        bics = [res_linear.bic, res_first_order.bic, res_second_order.bic]
        best_bic = np.argmin(bics)

        ########################################
        # create plot
        fig, ax = plt.subplots(4, figsize=(10, 14))

        # periodogram plot for residuals
        ax[0].semilogx(self.periods, self.power, 'k-', label='Data', zorder=5)
        ax[0].set_xlabel("Period [epoch]", fontsize=14)
        ax[0].set_ylabel('Power', fontsize=14)
        ax[0].axvline(per, color='red', label=f'Period: {per:.2f} epochs', alpha=0.75, zorder=10)
        # find power at closest period
        idx = np.argmin(np.abs(self.periods - per))
        per_power1 = self.power[idx]

        # find power at closest period
        idx = np.argmin(np.abs(self.periods - per2))
        per_power2 = self.power[idx]

        ax[0].set_title("Lomb-Scargle Periodogram", fontsize=18)
        ax[0].set_xlim([minper, (np.max(self.epochs) - np.min(self.epochs)) * 3.])

        # plot false alarm probability on lomb-scargle periodogram
        fp = ls.false_alarm_probability(power.max(), method='davies')
        fp_levels = ls.false_alarm_level([0.01, 0.05, 0.1], method='davies')

        # set upper y-limit on plot
        ax[0].set_ylim([0, self.power.max()])

        # plot as horizontal line
        ax[0].axhline(fp_levels[0], color='red', ls='--', label=f'99% FAP (Power = {fp_levels[0]:.1f})')

        # plot lomb-scargle for detrended data
        ax[0].semilogx(1. / freq2, power2 * 0.99, 'g-', alpha=0.5, label='Data - Fourier Fit 1', zorder=7)

        # plot false alarm probability on lomb-scargle periodogram
        fp = ls2.false_alarm_probability(power2.max(), method='davies')
        fp_levels2 = ls2.false_alarm_level([0.01, 0.05, 0.1], method='davies')

        # best period + false alarm for second order solution
        ax[0].axvline(per2, color='cyan', alpha=0.5, label=f'Period: {per2:.2f} epochs', zorder=10)
        ax[0].axhline(fp_levels2[0], color='cyan', ls='--', label=f'99% FAP (Power = {fp_levels2[0]:.1f})')

        # add horizontal dotted line at zero
        linear_label = f"Linear Fit (BIC: {res_linear.bic:.2f})"
        if best_bic == 0:
            linear_label += " best"
        ax[1].axhline(0, color='black', ls='--', label=linear_label)

        # super sample fourier solution for first order
        xnew = np.linspace(self.epochs.min(), self.epochs.max(), 1000)
        basis_new = np.ones((2, len(xnew)))
        basis_new[0] = np.sin(2 * np.pi * xnew / per)
        basis_new[1] = np.cos(2 * np.pi * xnew / per)
        y_bestfit_new = np.dot(basis_new.T, coeffs_first_order[2:])  # reconstruct signal

        # plot first order fourier solution
        fourier1_label = f'Fourier Fit 1st (BIC: {res_first_order.bic:.2f})'
        if per_power1 < fp_levels[0]:
            fourier1_label += " below 99% FAP"
        ax[1].plot(xnew, y_bestfit_new * 24 * 60, 'r-', label=fourier1_label, alpha=0.75)

        # set up ax labels
        ax[1].set_xlabel(f"Epochs", fontsize=14)
        ax[1].set_ylabel("O-C [min]", fontsize=14)
        ax[1].grid(True, ls='--')
        depoch = self.epochs.max() - self.epochs.min()
        ax[1].set_xlim([self.epochs.min() - depoch * 0.01, self.epochs.max() + depoch * 0.01])

        # o-c time series with fourier solution
        ax[1].errorbar(self.epochs, residuals_linear * 24 * 60,
                       yerr=self.dataerr * 24 * 60, ls='none',
                       marker='.', color='black', label=f'Data')

        # super sample fourier solution for second order
        xnew = np.linspace(self.epochs.min(), self.epochs.max(), 1000)
        basis_new = np.ones((4, len(xnew)))
        basis_new[0] = np.sin(2 * np.pi * xnew / per)
        basis_new[1] = np.cos(2 * np.pi * xnew / per)
        basis_new[2] = np.sin(4 * np.pi * xnew / per2)
        basis_new[3] = np.cos(4 * np.pi * xnew / per2)
        y_bestfit_new2 = np.dot(basis_new.T, coeffs_second_order[2:])  # reconstruct signal

        # plot first order fourier solution
        fourier2_label = f'Fourier Fit 2nd (BIC: {res_second_order.bic:.2f})'
        if per_power2 < fp_levels2[0]:  # should be fp_levels2?
            fourier2_label += " below 99% FAP"
        ax[1].plot(xnew, y_bestfit_new2 * 24 * 60, 'c-', label=fourier2_label, alpha=0.75)

        # set up ax labels
        ax[1].set_xlabel(f"Epochs", fontsize=14)
        ax[1].set_ylabel("O-C [min]", fontsize=14)
        ax[1].grid(True, ls='--')
        ax[1].legend(loc='best')
        ax[0].legend(loc='best')

        # plot phase folded signal for first order solution
        xnew = np.linspace(0, per, 1000)
        basis_new = np.ones((2, len(xnew)))
        basis_new[0] = np.sin(2 * np.pi * xnew / per)
        basis_new[1] = np.cos(2 * np.pi * xnew / per)
        y_bestfit_new = np.dot(basis_new.T, coeffs_first_order[2:])  # reconstruct signal
        xnewphase = xnew / per % 1
        si = np.argsort(xnewphase)

        # use uncertainty to derive fill between region
        std_dev = np.sqrt(np.diagonal(res_first_order.normalized_cov_params))
        samples_1st = []
        samples_2nd = []

        for i in range(1000):
            coeffs_1st = res_first_order.params + np.random.normal(0, 1, res_first_order.params.shape[0]) * std_dev
            samples_1st.append(np.dot(basis_new.T, coeffs_1st[2:]))  # TODO add in offset + period error
        samples_1st = np.array(samples_1st)

        # fill between region +/- 1 sigma
        ax[2].fill_between(xnewphase[si],
                           np.percentile(samples_1st, 16, axis=0)[si] * 24 * 60,
                           np.percentile(samples_1st, 84, axis=0)[si] * 24 * 60,
                           color='red', alpha=0.3, label='1-sigma region')

        # fill for 3 sigma
        ax[2].fill_between(xnewphase[si],
                           np.percentile(samples_1st, 0.15, axis=0)[si] * 24 * 60,
                           np.percentile(samples_1st, 99.85, axis=0)[si] * 24 * 60,
                           color='red', alpha=0.1, label='3-sigma region')

        # sort data in phase
        si = np.argsort(newphase)

        # plot data
        ax[2].errorbar(self.epochs / per % 1, residuals_linear * 24 * 60,
                       yerr=self.dataerr * 24 * 60, ls='none',
                       marker='.', color='black', label=f'Data')

        ax[2].set_xlabel(f"Phase (Period: {per:.2f} epochs)", fontsize=14)
        ax[2].set_ylabel("O-C [min]", fontsize=14)
        ax[2].grid(True, ls='--')
        ax[2].legend(loc='best')

        # find best fit signal with 2 periods
        # construct basis vectors with sin and cos
        basis2 = np.ones((5, len(self.epochs)))
        basis2[0] = np.sin(2 * np.pi * self.epochs / per)
        basis2[1] = np.cos(2 * np.pi * self.epochs / per)
        basis2[2] = np.sin(2 * np.pi * self.epochs / per2)
        basis2[3] = np.cos(2 * np.pi * self.epochs / per2)

        # perform the weighted least squares regression to find second order fourier solution
        res = sm.WLS(residuals_first_order, basis2.T, weights=1.0 / self.dataerr ** 2).fit()
        coeffs = res.params  # retrieve the slope and intercept of the fit from res
        y_bestfit = np.dot(basis2.T, coeffs)

        # super sample fourier solution
        xnew = np.linspace(self.epochs.min(), self.epochs.max(), 1000)
        basis_new = np.ones((5, len(xnew)))
        basis_new[1] = np.sin(2 * np.pi * xnew / per)
        basis_new[2] = np.cos(2 * np.pi * xnew / per)
        basis_new[3] = np.sin(2 * np.pi * xnew / per2)
        basis_new[4] = np.cos(2 * np.pi * xnew / per2)
        y_bestfit_new = np.dot(basis_new.T, coeffs)
        xnewphase = xnew / per2 % 1
        si = np.argsort(xnewphase)

        newphase = self.epochs / per2 % 1
        ax[3].errorbar(newphase, residuals_first_order * 24 * 60,
                       yerr=self.dataerr * 24 * 60, ls='none',
                       marker='.', color='black', label=f'Data - Fourier Fit 1')

        # create single sine wave from detrended data
        basis_new = np.ones((3, len(xnew)))
        basis_new[1] = np.sin(2 * np.pi * xnew / per)
        basis_new[2] = np.cos(2 * np.pi * xnew / per)
        y_best_single = np.dot(basis_new.T, coeffs[:3])

        # create best double sine wave from detrended data
        basis_new = np.ones((5, len(xnew)))
        basis_new[1] = np.sin(2 * np.pi * xnew / per)
        basis_new[2] = np.cos(2 * np.pi * xnew / per)
        basis_new[3] = np.sin(2 * np.pi * xnew / per2)
        basis_new[4] = np.cos(2 * np.pi * xnew / per2)
        y_best_double = np.dot(basis_new.T, coeffs)

        # use uncertainty to derive fill between region
        std_dev = np.sqrt(np.diagonal(res.normalized_cov_params))
        samples = []
        for i in range(1000):
            coeffs = res.params + np.random.normal(0, std_dev)
            ynew = np.dot(basis_new.T, coeffs)
            ynew -= y_best_single
            samples.append(ynew)
        samples = np.array(samples)

        y_bestfit_new = y_best_double - y_best_single

        ax[3].set_xlabel(f"Phase (Period: {per2:.2f} epochs)", fontsize=14)
        ax[3].set_ylabel("Residuals [min]", fontsize=14)
        ax[3].grid(True, ls='--')

        # fill between region +/- 1 sigma
        ax[3].fill_between(xnewphase[si],
                           np.percentile(samples, 16, axis=0)[si] * 24 * 60,
                           np.percentile(samples, 84, axis=0)[si] * 24 * 60,
                           color='cyan', alpha=0.4, label='1-sigma region')

        # fill for 3 sigma
        ax[3].fill_between(xnewphase[si],
                           np.percentile(samples, 0.15, axis=0)[si] * 24 * 60,
                           np.percentile(samples, 99.85, axis=0)[si] * 24 * 60,
                           color='cyan', alpha=0.2, label='3-sigma region')

        ax[3].legend(loc='best')

        return fig, ax

class non_linear_fitter(object):

    def __init__(self, data, dataerr, bounds=None, prior=None, labels=None, verbose=True):
        """
        Fit a Non-linear model to data using nested sampling.
        https://arxiv.org/pdf/2211.05646.pdf eq. 3

        Parameters
        ----------
        data : array
            Data to fit.
        dataerr : array
            Error on data.
        bounds : dict, optional
            Bounds on parameters. Dictionary of tuples/list
        prior : dict, optional
            Prior on parameters for slope and intercept. Dictionary of tuples/list
        """
        self.data = data
        self.dataerr = dataerr
        self.bounds = bounds
        self.labels = np.array(labels)
        self.prior = prior.copy()  # dict {'m':(0.1,0.5), 'b':(0,1)}
        self.verbose = verbose
        if bounds is None:
            # use +- 3 sigma prior as bounds
            # prior is list of tuples (mean, std)
            self.bounds = {
                'P': [prior['P'][0] - 3 * prior['P'][1], prior['P'][0] + 3 * prior['P'][1]],
                't0': [prior['t0'][0] - 3 * prior['t0'][1], prior['t0'][0] + 3 * prior['t0'][1]],
                'dPdN': [prior['dPdN'][0] - 3 * prior['dPdN'][1], prior['dPdN'][0] + 3 * prior['dPdN'][1]]
            }
        self.results = None
        self.fit_nested()

    def fit_nested(self):
        """ Fit a linear model to data using nested sampling. """

        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)

        # orbital epochs (N)
        self.epochs = np.round((self.data - np.mean(self.bounds['t0'])) / np.mean(self.bounds['P']))

        def loglike(pars):
            # chi-squared
            # tmid = t0 + N*P + 0.5*dPdN*N**2 (eq 3 from paper)
            model = pars[0] * self.epochs + pars[1] + 0.5 * pars[2] * self.epochs ** 2
            return -0.5 * np.sum(((self.data - model) / self.dataerr) ** 2)

        def prior_transform(upars):
            # transform unit cube to prior volume
            return (boundarray[:, 0] + bounddiff * upars)

        # estimate slope and intercept
        noop = lambda *args, **kwargs: None
        if self.verbose:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=4e5,
                                                                                         min_num_live_points=420,
                                                                                         show_status=True)
        else:
            self.results = ReactiveNestedSampler(freekeys, loglike, prior_transform).run(max_ncalls=4e5,
                                                                                         min_num_live_points=420,
                                                                                         show_status=False,
                                                                                         viz_callback=noop)
        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = {}

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # final model
        self.model = self.parameters['t0'] + self.epochs * self.parameters['P'] #+ 0.5 * self.parameters['dPdN'] * self.epochs ** 2
        self.residuals = self.data - self.model

    def plot_oc(self, savefile=None, ylim='none', show_2sigma=False, prior_name="Prior"):
        """ Plot the data in the form of residuals vs. time

        Parameters
        ----------
        savefile : str, optional
            Save the figure to a file.
        ylim : str, optional
            Set the y-axis limits. Default is 'none'. Can be prior, average, or none.
        show_2sigma : bool, optional
            Show a fill between using the 2 sigma limits. Default is False (aka 1 sigma)
        """

        # set up the figure        
        fig, ax = plt.subplots(1, figsize=(9, 6))

        # check if labels are not None
        if self.labels.size:
            # find unique set of labels
            ulabels = np.unique(self.labels)
            # set up a color/marker cycle
            markers = cycle(['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
            colors = cycle(['black', 'blue', 'green', 'orange', 'purple', 'grey', 'magenta', 'cyan', 'lime'])

            # plot each label separately
            for i, ulabel in enumerate(ulabels):
                # find where the label matches
                mask = self.labels == ulabel
                # plot the data/residuals
                ax.errorbar(self.epochs[mask], self.residuals[mask] * 24 * 60, yerr=self.dataerr[mask] * 24 * 60,
                            ls='none', marker=next(markers), color=next(colors), label=ulabel)
        else:
            # plot the data/residuals
            ax.errorbar(self.epochs, self.residuals * 24 * 60, yerr=self.dataerr * 24 * 60,
                        ls='none', marker='o', color='black')

        ylower = (self.residuals.mean() - 3 * np.std(self.residuals)) * 24 * 60
        yupper = (self.residuals.mean() + 3 * np.std(self.residuals)) * 24 * 60

        # upsample data
        epochs = (np.linspace(self.data.min() - 7, self.data.max() + 365*5, 1000) -
                  self.parameters['t0']) / self.parameters['P']

        # set the y-axis limits
        depoch = self.epochs.max() - self.epochs.min()
        #ax.set_xlim([self.epochs.min() - depoch * 0.01, self.epochs.max() + depoch * 0.01])

        # best fit solution
        model = epochs * self.parameters['P'] + self.parameters['t0']
        # nl_model = epochs * self.parameters['P'] + self.parameters['t0'] + 0.5 * self.parameters['dPdN'] * epochs ** 2

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['P'], self.errors['P'], size=10000)
        mc_b = np.random.normal(self.parameters['t0'], self.errors['t0'], size=10000)
        mc_dPdN = np.random.normal(self.parameters['dPdN'], self.errors['dPdN'], size=10000)
        mc_model = np.expand_dims(epochs, -1) * mc_m + mc_b + 0.5 * mc_dPdN * np.expand_dims(epochs, -1) ** 2

        # create a fill between area for uncertainty of new ephemeris
        diff = mc_model.T - model

        if show_2sigma:
            ax.fill_between(epochs, np.percentile(diff, 2, axis=0) * 24 * 60, np.percentile(diff, 98, axis=0) * 24 * 60,
                            alpha=0.2, color='k', label=r'Uncertainty ($\pm$ 2$\sigma$)')
        else:
            # show 1 sigma
            ax.fill_between(epochs, np.percentile(diff, 36, axis=0) * 24 * 60,
                            np.percentile(diff, 64, axis=0) * 24 * 60, alpha=0.2, color='k',
                            label=r'Uncertainty ($\pm$ 1$\sigma$)')

        # duplicate axis and plot days since mid-transit
        ax2 = ax.twiny()
        ax2.set_xlabel(f"Time [BJD - {self.parameters['t0']:.1f}]", fontsize=14)
        ax2.set_xlim(ax.get_xlim())
        xticks = ax.get_xticks()
        dt = np.round(xticks * self.parameters['P'], 1)
        # ax2.set_xticks(dt)
        ax2.set_xticklabels(dt)

        if ylim == 'diff':
            ax.set_ylim([min(np.percentile(diff, 1, axis=0) * 24 * 60),
                         max(np.percentile(diff, 99, axis=0) * 24 * 60)])

        # overlay the prior ephemeris
        if self.prior is not None:
            # create fill between area for uncertainty of old/prior ephemeris
            epochs_p = (np.linspace(self.data.min() - 7, self.data.max() + 365*5, 1000) -
                        self.prior['t0'][0]) / self.prior['P'][0]
            prior = epochs_p * self.prior['P'][0] + self.prior['t0'][0]
            mc_m_p = np.random.normal(self.prior['P'][0], self.prior['P'][1], size=10000)
            mc_b_p = np.random.normal(self.prior['t0'][0], self.prior['t0'][1], size=10000)
            #mc_dPdN_p = np.random.normal(self.prior['dPdN'][0], self.prior['dPdN'][1], size=10000)
            mc_model_p = np.expand_dims(epochs_p, -1) * mc_m_p + mc_b_p #+ 0.5 * mc_dPdN_p * np.expand_dims(epochs_p, -1) ** 2
            diff_p = mc_model_p.T - model

            # plot an invisible line so the 2nd axes are happy
            ax2.plot(epochs, (model - prior) * 24 * 60, ls='--', color='r', alpha=0)

            # why is this so small!?!?!? consistent to within machine precision?
            # ax.plot(epochs, (model-prior)*24*60, ls='--', color='r')

            if show_2sigma:
                ax.fill_between(epochs, np.percentile(diff_p, 2, axis=0) * 24 * 60,
                                np.percentile(diff_p, 98, axis=0) * 24 * 60, alpha=0.1, color='r',
                                label=rf'{prior_name} ($\pm$ 2$\sigma$)')
            else:
                # show ~1 sigma
                ax.fill_between(epochs, np.percentile(diff_p, 36, axis=0) * 24 * 60,
                                np.percentile(diff_p, 64, axis=0) * 24 * 60, alpha=0.1, color='r',
                                label=rf'{prior_name} ($\pm$ 1$\sigma$)')

            if ylim == 'prior':
                ax.set_ylim([min(np.percentile(diff_p, 1, axis=0) * 24 * 60),
                             max(np.percentile(diff_p, 99, axis=0) * 24 * 60)])
            elif ylim == 'average':
                ax.set_ylim([0.5 * (min(np.percentile(diff, 1, axis=0) * 24 * 60) +
                                    min(np.percentile(diff_p, 1, axis=0) * 24 * 60)),
                             0.5 * (max(np.percentile(diff, 99, axis=0) * 24 * 60) +
                                    max(np.percentile(diff_p, 99, axis=0) * 24 * 60))])

        
        # constant period model
        ax.axhline(0, color='black', alpha=0.5, ls='--',
                   label="Period: {:.7f}+-{:.7f} days\nT_mid: {:.7f}+-{:.7f} BJD".format(self.parameters['P'],
                                                                                         self.errors['P'],
                                                                                         self.parameters['t0'],
                                                                                         self.errors['t0']))
        
        # dPdN model
        model_diff = 0.5 * self.parameters['dPdN'] * epochs ** 2
        ax.plot(epochs, model_diff * 24 * 60, color='cyan', alpha=0.5, ls='--',
                label="dPdN: {:.2e}+-{:.2e} days/epoch".format(self.parameters['dPdN'], self.errors['dPdN']))


        ax.legend(loc='best')
        ax.set_xlabel("Epoch [number]", fontsize=14)
        ax.set_ylabel("Residuals [min]", fontsize=14)
        ax.grid(True, ls='--')
        return fig, ax

    def plot_triangle(self):
        """ Create a posterior triangle plot of the results. """
        ranges = []
        mask1 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        mask2 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        mask3 = np.ones(len(self.results['weighted_samples']['logl']), dtype=bool)
        titles = []
        labels = []
        flabels = {
            'P': 'Period [day]',
            't0': 'T_mid [JD]',
            'dPdN': 'dPdN [day/epoch]'
        }
        for i, key in enumerate(self.quantiles):
            labels.append(flabels.get(key, key))
            titles.append(f"{self.parameters[key]:.7f} +-\n {self.errors[key]:.7f}")

            # set the axes limits for the plots
            ranges.append([
                self.parameters[key] - 5 * self.errors[key],
                self.parameters[key] + 5 * self.errors[key]
            ])

            if key == 'a2' or key == 'a1':
                continue

            # create masks for contouring on sigma bounds
            mask3 = (mask3 &
                     (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - 3 * self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 3 * self.errors[key])))

            mask1 = (mask1 &
                     (self.results['weighted_samples']['points'][:, i] > (self.parameters[key] - self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + self.errors[key])))

            mask2 = (mask2 &
                     (self.results['weighted_samples']['points'][:, i] > ( self.parameters[key] - 2 * self.errors[key])) &
                     (self.results['weighted_samples']['points'][:, i] < (self.parameters[key] + 2 * self.errors[key])))

        chi2 = self.results['weighted_samples']['logl'] * -2
        fig = corner(self.results['weighted_samples']['points'],
                     labels=labels,
                     bins=int(np.sqrt(self.results['samples'].shape[0])),
                     range=ranges,
                     figsize=(10, 10),
                     # quantiles=(0.1, 0.84),
                     plot_contours=True,
                     levels=[np.percentile(chi2[mask1], 95), np.percentile(chi2[mask2], 95),
                             np.percentile(chi2[mask3], 95)],
                     plot_density=False,
                     titles=titles,
                     data_kwargs={
                         'c': chi2,  # color code by chi2
                         'vmin': np.percentile(chi2[mask3], 1),
                         'vmax': np.percentile(chi2[mask3], 95),
                         'cmap': 'viridis',
                     },
                     label_kwargs={'labelpad': 50, },
                     hist_kwargs={'color': 'black', }
                     )
        return fig
