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


class ephemeris_fitter(object):

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
        if labels is not None:
            self.labels = np.array(labels)
        else:
            self.labels = np.array(["Data"]*len(data))
        self.prior = prior.copy()  # dict {'P':(0.1,0.5), 'T0':(0,1), (value, uncertainty)}
        self.verbose = verbose

        # if no bounds use prior
        if bounds is None:
            # use +- 5 sigma prior as bounds
            self.bounds = {
                'P': [prior['P'][0] - 5 * prior['P'][1], prior['P'][0] + 5 * prior['P'][1]],
                'T0': [prior['T0'][0] - 5 * prior['T0'][1], prior['T0'][0] + 5 * prior['T0'][1]]
            }
        self.results = None
        self.fit_nested()

    def fit_nested(self):
        """ Fit a linear model to data using nested sampling. """
        # order matters :/
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray, 1).reshape(-1)
        self.epochs = np.round((self.data - np.mean(self.bounds['T0'])) / np.mean(self.bounds['P']))

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
        self.model = self.epochs * self.parameters['P'] + self.parameters['T0']
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
        if self.labels is not None:
            # find unique set of labels
            ulabels = np.unique(self.labels)
            # set up a color/marker cycle
            markers = cycle(['o', 'v', '^', '<', '>', 's', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
            colors = cycle(['black', 'blue', 'green', 'orange', 'purple', 'grey', 'magenta', 'cyan', 'lime'])

            # plot each label separately
            for i, ulabel in enumerate(ulabels):
                # find where the label matches
                mask = self.labels == ulabel
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
                  self.parameters['T0']) / self.parameters['P']

        # set the y-axis limits
        depoch = self.epochs.max() - self.epochs.min()
        ax.set_xlim([self.epochs.min() - depoch * 0.01, self.epochs.max() + depoch * 0.01])

        # best fit solution
        model = epochs * self.parameters['P'] + self.parameters['T0']

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['P'], self.errors['P'], size=10000)
        mc_b = np.random.normal(self.parameters['T0'], self.errors['T0'], size=10000)
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
        ax2.set_xlabel(f"Time [BJD - {self.parameters['T0']:.1f}]", fontsize=14)
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
            epochs_p = (np.linspace(self.data.min() - 7, self.data.max() + 7, 1000) -
                        self.prior['T0'][0]) / self.prior['P'][0]
            prior = epochs_p * self.prior['P'][0] + self.prior['T0'][0]
            mc_m_p = np.random.normal(self.prior['P'][0], self.prior['P'][1], size=10000)
            mc_b_p = np.random.normal(self.prior['T0'][0], self.prior['T0'][1], size=10000)
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
                   label="Period: {:.7f}+-{:.7f} days\nT_mid: {:.7f}+-{:.7f} BJD".format(self.parameters['P'],
                                                                                         self.errors['P'],
                                                                                         self.parameters['T0'],
                                                                                         self.errors['T0']))

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
            'P': 'Period [day]',
            'T0': 'T_mid [JD]',
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

    def plot_periodogram(self, minper=0, maxper=0, minper2=0, maxper2=0):
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
            minper = max(3, 2. * np.diff(self.epochs[si]).min())
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
        basis[2] = np.sin(2. * np.pi * self.epochs / per)
        basis[3] = np.cos(2. * np.pi * self.epochs / per)

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
        if minper2 == 0:
            minper2 = per*2.1
        if maxper2 == 0:
            maxper2 = (np.max(self.epochs) - np.min(self.epochs)) * 2.

        # reset minper2 if it is greater than maxper2
        if minper2 > maxper2:
            minper2 = 3.1

        # recompute on new grid
        ls2 = LombScargle(self.epochs, residuals_linear, dy=self.dataerr)

        # search for periods greater than first order
        freq2,power2 = ls.autopower(maximum_frequency=1./(minper2+1),
                                    minimum_frequency=1./maxper2, nyquist_factor=2)

        # TODO do a better job defining period grid, ignoring harmonics of first order solution +- 0.25 day

        # find max period
        try:
            mi2 = np.argmax(power2)
        except:
            import pdb; pdb.set_trace()

        per2 = 1. / freq2[mi2]

        # create basis vectors for second order solution
        basis = np.ones((6, len(self.epochs)))
        basis[1] = self.epochs
        basis[2] = np.sin(2. * np.pi * self.epochs / per)
        basis[3] = np.cos(2. * np.pi * self.epochs / per)
        basis[4] = np.sin(2. * np.pi * self.epochs / per2)
        basis[5] = np.cos(2. * np.pi * self.epochs / per2)

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
        basis_new[0] = np.sin(2. * np.pi * xnew / per)
        basis_new[1] = np.cos(2. * np.pi * xnew / per)
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
        basis_new[0] = np.sin(2. * np.pi * xnew / per)
        basis_new[1] = np.cos(2. * np.pi * xnew / per)
        basis_new[2] = np.sin(2. * np.pi * xnew / per2)
        basis_new[3] = np.cos(2. * np.pi * xnew / per2)
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
        basis_new[0] = np.sin(2. * np.pi * xnew / per)
        basis_new[1] = np.cos(2. * np.pi * xnew / per)
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
        basis2[0] = np.sin(2. * np.pi * self.epochs / per)
        basis2[1] = np.cos(2. * np.pi * self.epochs / per)
        basis2[2] = np.sin(2. * np.pi * self.epochs / per2)
        basis2[3] = np.cos(2. * np.pi * self.epochs / per2)

        # perform the weighted least squares regression to find second order fourier solution
        res = sm.WLS(residuals_first_order, basis2.T, weights=1.0 / self.dataerr ** 2).fit()
        coeffs = res.params  # retrieve the slope and intercept of the fit from res
        self.coeffs = res.params
        y_bestfit = np.dot(basis2.T, coeffs)

        # super sample fourier solution
        xnew = np.linspace(self.epochs.min(), self.epochs.max(), 1000)
        basis_new = np.ones((5, len(xnew)))
        basis_new[1] = np.sin(2. * np.pi * xnew / per)
        basis_new[2] = np.cos(2. * np.pi * xnew / per)
        basis_new[3] = np.sin(2. * np.pi * xnew / per2)
        basis_new[4] = np.cos(2. * np.pi * xnew / per2)
        y_bestfit_new = np.dot(basis_new.T, coeffs)
        xnewphase = xnew / per2 % 1
        si = np.argsort(xnewphase)

        newphase = self.epochs / per2 % 1
        ax[3].errorbar(newphase, residuals_first_order * 24 * 60,
                       yerr=self.dataerr * 24 * 60, ls='none',
                       marker='.', color='black', label=f'Data - Fourier Fit 1')

        # create single sine wave from detrended data
        basis_new = np.ones((3, len(xnew)))
        basis_new[1] = np.sin(2. * np.pi * xnew / per)
        basis_new[2] = np.cos(2. * np.pi * xnew / per)
        y_best_single = np.dot(basis_new.T, coeffs[:3])

        # create best double sine wave from detrended data
        basis_new = np.ones((5, len(xnew)))
        basis_new[1] = np.sin(2. * np.pi * xnew / per)
        basis_new[2] = np.cos(2. * np.pi * xnew / per)
        basis_new[3] = np.sin(2. * np.pi * xnew / per2)
        basis_new[4] = np.cos(2. * np.pi * xnew / per2)
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

        self.best_periods = np.array([per, per2])
        coeff1 = np.sqrt(coeffs[1]**2 + coeffs[1]**2)
        coeff2 = np.sqrt(coeffs[3]**2 + coeffs[4]**2)
        self.amplitudes = np.array([coeff1, coeff2])

        return fig, ax

class decay_fitter(object):

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
        if labels is not None:
            self.labels = np.array(labels)
        else:
            self.labels = ["Data"]*len(data)
        self.prior = prior.copy()  # dict {'m':(0.1,0.5), 'b':(0,1)}
        self.verbose = verbose
        if bounds is None:
            # use +- 3 sigma prior as bounds
            # prior is list of tuples (mean, std)
            self.bounds = {
                'P': [prior['P'][0] - 3 * prior['P'][1], prior['P'][0] + 3 * prior['P'][1]],
                'T0': [prior['T0'][0] - 3 * prior['T0'][1], prior['T0'][0] + 3 * prior['T0'][1]],
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
        self.epochs = np.round((self.data - np.mean(self.bounds['T0'])) / np.mean(self.bounds['P']))

        def loglike(pars):
            # chi-squared
            # tmid = T0 + N*P + 0.5*dPdN*N**2 (eq 3 from paper)
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
        self.model = self.parameters['T0'] + self.epochs * self.parameters['P'] #+ 0.5 * self.parameters['dPdN'] * self.epochs ** 2
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
        if self.labels is not None:
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
        # SET limits on x-axis here
        epochs = (np.linspace(self.data.min() - 7, self.data.max() + 365*5, 1000) -
                  self.parameters['T0']) / self.parameters['P']

        # set the y-axis limits
        depoch = self.epochs.max() - self.epochs.min()
        #ax.set_xlim([self.epochs.min() - depoch * 0.01, self.epochs.max() + depoch * 0.01])

        # best fit solution
        model = epochs * self.parameters['P'] + self.parameters['T0']
        # nl_model = epochs * self.parameters['P'] + self.parameters['T0'] + 0.5 * self.parameters['dPdN'] * epochs ** 2

        # MonteCarlo the new ephemeris for uncertainty
        mc_m = np.random.normal(self.parameters['P'], self.errors['P'], size=10000)
        mc_b = np.random.normal(self.parameters['T0'], self.errors['T0'], size=10000)
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
        ax2.set_xlabel(f"Time [BJD - {self.parameters['T0']:.1f}]", fontsize=14)
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
                        self.prior['T0'][0]) / self.prior['P'][0]
            prior = epochs_p * self.prior['P'][0] + self.prior['T0'][0]
            mc_m_p = np.random.normal(self.prior['P'][0], self.prior['P'][1], size=10000)
            mc_b_p = np.random.normal(self.prior['T0'][0], self.prior['T0'][1], size=10000)
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
                                                                                         self.parameters['T0'],
                                                                                         self.errors['T0']))
        
        # dPdN model
        model_diff = 0.5 * self.parameters['dPdN'] * epochs ** 2
        ax.plot(epochs, model_diff * 24 * 60, color='cyan', alpha=0.5, ls='--',
                label="dP/dN: {:.2e}+-{:.2e} days/epoch".format(self.parameters['dPdN'], self.errors['dPdN']))

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
            'T0': 'T_mid [JD]',
            'dPdN': 'dPdN [day/epoch]'
        }
        for i, key in enumerate(self.quantiles):
            labels.append(flabels.get(key, key))
            if key == 'dPdN':
                titles.append(f"{self.parameters[key]:.2e} +-\n {self.errors[key]:.2e}")
            else:
                titles.append(f"{self.parameters[key]:.8f} +-\n {self.errors[key]:.8f}")

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

if __name__ == "__main__":

    Tc = np.array([ # measured mid-transit times
        #Exoplanet Watch transit times array
        2457025.7867,2457347.7655,2457694.8263,2458112.8495,2458492.6454,2459532.7834,2459604.81036,2459604.8137,2459614.6365,2459616.8209,
        2459651.7466,2459903.8588,2459914.7783,2459915.8684,2459925.6921,2459939.8793,2459949.7047,2459959.5249,2459962.8037,2459973.7129,
        2459974.798,2459986.8057,2459994.4489,2460009.7312,2458867.01587,2459501.1295,

        #ExoClock transit times array
        2454508.977,2454515.525,2454801.477,2454840.769,2455123.447,2455140.91,2455147.459,2455147.459,2455148.551,2455148.552,2455158.372,
        2455158.373,2455159.464,2455159.467,2455160.556,2455163.831,2455172.562,2455192.205,2455209.669,2455210.762,2455215.129,2455230.407,
        2455238.046,2455254.419,2455265.331,2455494.53,2455494.531,2455509.81,2455510.902,2455530.547,2455542.553,2455548.011,2455565.472,
        2455566.563,2455575.296,2455575.297,2455578.569,2455589.483,2455590.576,2455598.216,2455600.397,2455600.398,2455600.398,2455601.49,
        2455601.49,2455601.49,2455603.673,2455612.404,2455623.318,2455623.319,2455624.41,2455624.41,2455663.702,2455842.693,2455842.694,
        2455876.528,2455887.442,2455887.442,2455887.442,2455888.533,2455888.534,2455888.534,2455889.625,2455890.716,2455901.63,2455903.814,
        2455903.814,2455920.185,2455923.459,2455926.733,2455939.829,2455945.287,2455946.378,2455946.379,2455947.47,2455948.561,2455951.835,
        2455952.927,2455959.475,2455960.567,2455970.39,2455971.481,2455981.304,2455982.395,2455983.487,2455984.578,2455985.67,2455996.584,
        2456000.95,2456005.315,2456006.406,2456014.046,2456175.577,2456245.427,2456249.794,2456273.805,2456282.536,2456284.719,2456297.816,
        2456302.182,2456305.455,2456319.644,2456328.376,2456329.467,2456604.505,2456605.596,2456606.688,2456607.779,2456629.607,2456630.699,
        2456654.71,2456659.076,2456662.35,2456663.441,2456664.533,2456674.356,2456677.63,2456688.544,2456694.002,2456703.824,2456711.464,
        2456719.104,2456721.287,2456722.378,2456986.502,2457010.513,2457012.696,2457033.434,2457045.438,2457046.53,2457059.627,2457060.718,
        2457067.267,2457068.358,2457103.284,2457345.579,2457368.5,2457378.321,2457390.327,2457391.418,2457426.343,2457426.344,2457427.435,
        2457451.446,2457474.364,2457486.372,2457671.913,2457691.559,2457703.564,2457706.838,2457726.484,2457727.575,2457760.318,2457765.775,
        2457766.866,2457772.324,2457773.415,2457776.689,2457786.512,2457788.695,2457800.7,2457808.34,2457809.432,2457809.434,2457810.523,
        2457832.348,2457832.351,2458026.624,2458050.635,2458060.459,2458072.462,2458073.555,2458074.647,2458077.921,2458123.76,2458124.852,
        2458132.491,2458134.675,2458136.858,2458155.41,2458155.412,2458156.503,2458159.778,2458166.326,2458178.331,2458214.347,2458411.895,
        2458467.557,2458471.923,2458478.472,2458479.564,2458490.477,2458494.843,2458501.39,2458501.392,2458506.848,2458536.315,2458537.408,
        2458871.382,2458880.113,2458883.384,2458893.209,2458903.032,2459088.575,2459148.602,2459156.242,2459157.334,2459190.076,2459194.441,
        2459194.443,2459196.623,2459197.717,2459204.264,

        #ETD mid transit times array
        2454508.976854,2454836.403393,2454840.768603,2454840.771193,2454908.437992,2454931.358182,2455164.923956,2455164.924316,2455172.561816,2455172.562786,
        2455198.756735,2455230.406730,2455254.418870,2455256.596033,2455265.334053,2455506.533285,2455509.808915,2455519.632074,2455532.729964,2455541.461084,
        2455576.387242,2455577.475172,2455577.476782,2455577.477112,2455578.568422,2455593.850192,2455600.398282,2455601.489211,2455603.671161,2455612.400801,
        2455612.401321,2455615.675911,2455647.328760,2455857.973783,2455857.973783,2455867.797832,2455873.253612,2455877.618482,2455877.619162,2455889.623302,
        2455899.448441,2455899.451251,2455900.539781,2455904.905271,2455905.995241,2455912.544561,2455912.544561,2455913.636131,2455914.726741,2455922.370121,
        2455928.913321,2455936.552560,2455936.561520,2455937.645430,2455938.738730,2455946.383030,2455947.468480,2455949.650900,2455949.650900,2455950.742510,
        2455957.287600,2455957.294400,2455957.296780,2455957.296780,2455961.659130,2455962.749930,2455970.392369,2455982.394919,2455993.309869,2455993.313119,
        2455994.398719,2455994.399359,2456003.129689,2456003.130169,2456005.315328,2456222.509695,2456246.516835,2456249.793845,2456269.425195,2456271.623265,
        2456296.724384,2456304.364534,2456304.365804,2456304.367454,2456304.367784,2456329.469294,2456364.392544,2456584.860173,2456595.772813,2456596.866383,
        2456604.503073,2456608.869564,2456616.511354,2456626.333324,2456627.429034,2456628.515004,2456628.517074,2456630.699754,2456639.430344,2456639.432894,
        2456640.522554,2456641.610274,2456650.341584,2456663.443614,2456688.544734,2456722.372535,2456722.382615,2456734.383505,2456734.388675,2456746.385515,
        2456746.391885,2456927.546138,2456950.485119,2456960.306869,2456963.578439,2456963.578959,2456963.582479,2456986.501140,2456998.508180,2456998.508180,
        2457032.333311,2457033.433561,2457069.448912,2457080.367572,2457080.367572,2457092.368353,2457344.486723,2457344.488693,2457345.578813,2457345.579833,
        2457356.493733,2457357.585013,2457357.586493,2457368.498964,2457368.499654,2457369.589474,2457379.414424,2457379.414424,2457414.334726,2457426.343356,
        2457451.447237,2457474.364078,2457474.364348,2457486.371248,2457668.637786,2457726.484138,2457749.405049,2457750.494949,2457751.585649,2457751.586939,
        2457757.043281,2457758.134851,2457760.317471,2457760.318171,2457761.408951,2457770.140361,2457772.327372,2457773.415162,2457774.506522,2457809.431273,
        2457809.435043,2457832.350874,2458087.745843,2458096.474134,2458109.571734,2458131.400315,2458396.615162,2458396.617562,2458411.896052,2458455.554593,
        2458457.734253,2458479.568743,2458489.383454,2458489.387804,2458490.477904,2458490.482824,2458501.391674,2458501.392014,2458501.396544,2458513.402254,
        2458537.407634,2458538.507734,2458837.550156,2458848.461686,2458859.375956,2458860.465876,2458860.469616,2458871.379136,2458871.382556,2458872.470796,
        2458883.387165,2458883.389555,2458884.478175,2458884.478425,2458884.479455,2458906.306925,2458930.332215,2459147.510251,2459171.522921,2459172.608741,
        2459172.617851,2459183.531910,2459193.352480,2459194.442340,2459194.443610,2459196.623950,2459197.719250,2459220.632699,2459242.461499,2459265.384258,
        2459265.384318,2459505.495609,2459553.514627,2459553.515137,2459564.427746,2459565.520376,2459566.615006,2459576.438096])

    Tc_error = np.array([
        0.0039,0.0058,0.0054,0.0022,0.0065,0.0024,0.00091,0.00097,0.0022,0.0018,
        0.002,0.0021,0.0065,0.0036,0.0071,0.0031,0.0044,0.0021,0.0028,0.0023,
        0.0039,0.0071,0.0015,0.0055,0.00092,0.00092,

        #ExoClock error array
        0.0002,0.00013,0.00041,0.00047,0.0007,0.000417,0.0009,0.00042,0.001,0.0013,0.00057,
        0.0013,0.00095,0.00061,0.0007,0.000324,0.00044,0.0005,0.000463,0.000405,0.00082,0.00011,
        0.00066,0.00014,0.0016,0.00048,0.00063,0.00037,0.000313,0.00077,0.00029,0.00024,0.0011,
        0.00022,0.00063,0.00079,0.0008,0.00094,0.001,0.00036,0.00083,0.00026,0.00078,0.000289,
        0.00017,0.00039,0.001,0.00088,0.0009,0.00054,0.00086,0.0014,0.00063,0.00093,0.0013,
        0.0003,0.00038,0.00062,0.00075,0.00059,0.0003,0.00038,0.0009,0.00023,0.00083,0.00093,
        0.000324,0.00046,0.0002,0.00092,0.0019,0.00079,0.00036,0.00034,0.00022,0.00028,0.00011,
        0.0001,0.00018,0.0004,0.00029,0.00029,0.00094,0.00047,0.00029,0.000324,0.000417,0.00037,
        0.0004,0.00038,0.0004,0.00023,0.00033,0.00033,0.000394,0.000301,0.0003,0.000301,0.000301,
        0.00046,0.00026,0.000382,0.00027,0.00029,0.0002,0.0003,0.00034,0.000706,0.00019,0.00043,
        0.000336,0.00034,0.00019,0.00019,0.00032,0.00028,0.000324,0.00041,0.00029,0.00029,0.00026,
        0.00034,0.00034,0.00046,0.00043,0.00039,0.000486,0.0005,0.00049,0.00049,0.000347,0.000359,
        0.00022,0.00021,0.0003,0.00042,0.0004,0.0013,0.00034,0.00033,0.00055,0.0006,0.00023,0.00021,
        0.0007,0.0013,0.00035,0.00025,0.00034,0.00037,0.00028,0.00023,0.0006,0.00028,0.00039,
        0.00024,0.00022,0.00029,0.00026,0.00048,0.00032,0.0004,0.00018,0.0009,0.00021,0.0006,
        0.0006,0.00056,0.00023,0.0003,0.0003,0.00022,0.00034,0.00028,0.00027,0.00035,0.00031,
        0.00032,0.00033,0.0005,0.00031,0.00032,0.00091,0.00034,0.00038,0.0017,0.0004,0.0005,
        0.00026,0.0006,0.0006,0.0008,0.0003,0.0009,0.0003,0.00044,0.0008,0.0007,0.0009,
        0.0003,0.0007,0.0005,0.0003,0.0013,0.0007,0.0003,0.0004,0.0003,0.0012,0.0006,
        0.0005,0.0013,0.0004,

        #ETD error array
        0.000200,0.000600,0.000470,0.001000,0.001000,0.000980,0.001490,0.001290,0.000440,0.000140,
        0.001410,0.000110,0.000140,0.000590,0.001290,0.001120,0.000870,0.000700,0.000350,0.000500,
        0.001350,0.000380,0.000690,0.000850,0.000820,0.000470,0.000420,0.001090,0.000940,0.000830,
        0.000780,0.000540,0.001730,0.000710,0.000710,0.000750,0.001260,0.000920,0.001290,0.000730,
        0.000630,0.000570,0.000380,0.001030,0.001350,0.000570,0.000570,0.000780,0.000470,0.000900,
        0.000730,0.000910,0.001230,0.001190,0.000720,0.000770,0.001020,0.000590,0.000590,0.000660,
        0.001100,0.001040,0.000570,0.000570,0.001070,0.001320,0.000860,0.001160,0.000600,0.000760,
        0.000680,0.000760,0.000630,0.000600,0.000440,0.000810,0.000740,0.000670,0.000900,0.000550,
        0.000520,0.001460,0.000890,0.001560,0.000580,0.001640,0.001170,0.000510,0.000960,0.000510,
        0.000920,0.000710,0.000900,0.000510,0.001050,0.000970,0.000880,0.000440,0.000740,0.000680,
        0.000820,0.000800,0.001040,0.000760,0.000540,0.000890,0.001080,0.001120,0.000730,0.001390,
        0.001410,0.001640,0.000740,0.000570,0.000930,0.000890,0.000610,0.000510,0.001450,0.001450,
        0.001130,0.000460,0.000510,0.001190,0.001190,0.000600,0.000650,0.001070,0.001090,0.000490,
        0.000610,0.000520,0.000430,0.000580,0.000460,0.000340,0.000890,0.000890,0.001310,0.000650,
        0.000860,0.000770,0.000550,0.001350,0.000820,0.000550,0.000840,0.000530,0.000940,0.000720,
        0.000370,0.000520,0.000670,0.001030,0.000630,0.000330,0.000990,0.000550,0.000620,0.000720,
        0.000940,0.000580,0.000350,0.000570,0.001380,0.000640,0.001180,0.000700,0.000700,0.000710,
        0.000550,0.000580,0.000680,0.001030,0.000860,0.000850,0.000300,0.000760,0.000680,0.000820,
        0.000610,0.001450,0.000730,0.000700,0.001350,0.000820,0.000670,0.000940,0.000490,0.001130,
        0.000540,0.000540,0.000700,0.000790,0.000840,0.000600,0.000520,0.000730,0.000640,0.001020,
        0.000780,0.000610,0.001330,0.000770,0.000610,0.000520,0.001130,0.001130,0.000530,0.000780,
        0.000420,0.001250,0.000380,0.000720,0.000860,0.000470,0.000950,0.000540])

    labels = np.full(len(Tc), 'Data')

    P = 1.0914203  # orbital period for your target

    Tc_norm = Tc - Tc.min()  # normalize the data to the first observation
    orbit = np.rint(Tc_norm / P)  # number of orbits since first observation (rounded to nearest integer)

    # make a n x 2 matrix with 1's in the first column and values of orbit in the second
    A = np.vstack([np.ones(len(Tc)), orbit]).T

    # perform the weighted least squares regression
    res = sm.WLS(Tc, A, weights=1.0 / Tc_error ** 2).fit()
    # use sm.WLS for weighted LS, sm.OLS for ordinary LS, or sm.GLS for general LS

    params = res.params  # retrieve the slope and intercept of the fit from res
    std_dev = np.sqrt(np.diagonal(res.normalized_cov_params))

    slope = params[1]
    slope_std_dev = std_dev[1]
    intercept = params[0]
    intercept_std_dev = std_dev[0]

    # 3 sigma clip based on residuals
    calculated = orbit * slope + intercept
    residuals = (Tc - calculated) / Tc_error
    mask = np.abs(residuals) < 3
    Tc = Tc[mask]
    Tc_error = Tc_error[mask]
    labels = labels[mask]

    # print(res.summary())
    # print("Params =",params)
    # print("Error matrix =",res.normalized_cov_params)
    # print("Standard Deviations =",std_dev)

    print("Weighted Linear Least Squares Solution")
    print("T0 =", intercept, "+-", intercept_std_dev)
    print("P =", slope, "+-", slope_std_dev)

    # min and max values to search between for fitting
    bounds = {
        'P': [slope - 0.1, slope + 0.1],  # orbital period
        'T0': [intercept - 0.1, intercept + 0.1]  # mid-transit time
    }

    # used to plot red overlay in O-C figure
    prior = {
        'P': [slope, slope_std_dev],  # value from WLS (replace with literature value)
        'T0': [intercept, intercept_std_dev]  # value from WLS (replace with literature value)
    }

    lf = ephemeris_fitter(Tc, Tc_error, bounds, prior=prior, labels=labels)

    lf.plot_triangle()
    plt.subplots_adjust(top=0.9, hspace=0.2, wspace=0.2)
    plt.savefig("posterior.png")
    plt.close()
    print("image saved to: posterior.png")

    fig, ax = lf.plot_oc()
    plt.tight_layout()
    plt.savefig("oc.png")
    plt.show()
    plt.close()
    print("image saved to: oc.png")

    fig, ax = lf.plot_periodogram()
    plt.tight_layout()
    plt.savefig("periodogram.png")
    plt.close()
    print("image saved to: periodogram.png")

    # min and max values to search between for fitting
    bounds = {
        'P': [P - 0.1, P + 0.1],                   # orbital period
        'T0': [intercept - 0.1, intercept + 0.1],  # mid-transit time
        'dPdN': [-0.00001, 0.00001]  # dPdN
    }

    # used to plot red overlay in O-C figure
    prior = {
        'P': [slope, slope_std_dev],  # value from WLS (replace with literature value)
        'T0': [intercept, intercept_std_dev],  # value from WLS (replace with literature value)
        'dPdN': [0, 0.0001]  # best guess
    }

    nlf = decay_fitter(Tc, Tc_error, bounds, prior=prior, labels=labels)

    nlf.plot_triangle()
    plt.subplots_adjust(top=0.9, hspace=0.2, wspace=0.2)
    plt.savefig("nl_posterior.png")
    plt.close()
    print("image saved to: nl_posterior.png")

    nfig, ax = nlf.plot_oc()
    plt.tight_layout()
    plt.savefig("nl_oc.png")
    plt.show()
    plt.close()
    print("image saved to: nl_oc.png")

    # for each free key print parameter and error
    for key in nlf.parameters:
        print(f"Parameter {key} = {nlf.parameters[key]:.2e} +- {nlf.errors[key]:.2e}")

    # TODO BIC Values