import copy
import numpy as np
import matplotlib.pyplot as plt

import dynesty
from dynesty import plotting
from dynesty.utils import resample_equal

from scipy.ndimage import gaussian_filter as norm_kde
from scipy.stats import gaussian_kde

# Computes Hasting's polynomial approximation for the complete
# elliptic integral of the first (ek) and second (kk) kind
def ellke(k):
    m1=1.-k**2
    np.logm1 = np.log(m1)

    a1=0.44325141463
    a2=0.06260601220
    a3=0.04757383546
    a4=0.01736506451
    b1=0.24998368310
    b2=0.09200180037
    b3=0.04069697526
    b4=0.00526449639
    ee1=1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-np.logm1)
    ek = ee1+ee2

    a0=1.38629436112
    a1=0.09666344259
    a2=0.03590092383
    a3=0.03742563713
    a4=0.01451196212
    b0=0.5
    b1=0.12498593597
    b2=0.06880248576
    b3=0.03328355346
    b4=0.00441787012
    ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*np.logm1
    kk = ek1-ek2

    return [ek,kk]

# Computes the complete elliptical integral of the third kind using
# the algorithm of Bulirsch (1965):
def ellpic_bulirsch(n,k):
    kc=np.sqrt(1.-k**2); la=n+1.
    if(min(la) < 0.):
        print('Negative l')
    m0=1.; c=1.; la=np.sqrt(la); d=1./la; e=kc
    while 1:
        f = c; c = d/la+c; g = e/la; d = 2.*(f*g+d)
        la = g + la; g = m0; m0 = kc + m0
        if max(abs(1.-kc/g)) > 1.e-8:
            kc = 2*np.sqrt(e); e=kc*m0
        else:
            return 0.5*np.pi*(c*m0+d)/(m0*(m0+la))

#   Python translation of IDL code.
#   This routine computes the lightcurve for occultation of a
#   quadratically limb-darkened source without microlensing.  Please
#   cite Mandel & Agol (2002) and Eastman & Agol (2008) if you make use
#   of this routine in your research.  Please report errors or bugs to
#   jdeast@astronomy.ohio-state.edu
def occultquad(z,u1,u2,p0):

    nz = np.size(z)
    lambdad = np.zeros(nz)
    etad = np.zeros(nz)
    lambdae = np.zeros(nz)
    omega=1.-u1/3.-u2/6.

    ## tolerance for double precision equalities
    ## special case integrations
    tol = 1e-14

    p = abs(p0)

    z = np.where(abs(p-z) < tol,p,z)
    z = np.where(abs((p-1)-z) < tol,p-1.,z)
    z = np.where(abs((1-p)-z) < tol,1.-p,z)
    z = np.where(z < tol,0.,z)

    x1=(p-z)**2.
    x2=(p+z)**2.
    x3=p**2.-z**2.

    ## trivial case of no planet
    if p <= 0.:
        muo1 = np.zeros(nz) + 1.
        mu0  = np.zeros(nz) + 1.
        return [muo1,mu0]

    ## Case 1 - the star is unocculted:
    ## only consider points with z lt 1+p
    notusedyet = np.where( z < (1. + p) )
    notusedyet = notusedyet[0]
    if np.size(notusedyet) == 0:
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+ \
                  u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]

    # Case 11 - the  source is completely occulted:
    if p >= 1.:
        occulted = np.where(z[notusedyet] <= p-1.)#,complement=notused2)
        if np.size(occulted) != 0:
            ndxuse = notusedyet[occulted]
            etad[ndxuse] = 0.5 # corrected typo in paper
            lambdae[ndxuse] = 1.
            # lambdad = 0 already
            notused2 = np.where(z[notusedyet] > p-1)
            if np.size(notused2) == 0:
                muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.* \
                                                 (p > z))+u2*etad)/omega
                mu0=1.-lambdae
                return [muo1,mu0]
            notusedyet = notusedyet[notused2]

    # Case 2, 7, 8 - ingress/egress (uniform disk only)
    inegressuni = np.where((z[notusedyet] >= abs(1.-p)) & (z[notusedyet] < 1.+p))
    if np.size(inegressuni) != 0:
        ndxuse = notusedyet[inegressuni]
        tmp = (1.-p**2.+z[ndxuse]**2.)/2./z[ndxuse]
        tmp = np.where(tmp > 1.,1.,tmp)
        tmp = np.where(tmp < -1.,-1.,tmp)
        kap1 = np.arccos(tmp)
        tmp = (p**2.+z[ndxuse]**2-1.)/2./p/z[ndxuse]
        tmp = np.where(tmp > 1.,1.,tmp)
        tmp = np.where(tmp < -1.,-1.,tmp)
        kap0 = np.arccos(tmp)
        tmp = 4.*z[ndxuse]**2-(1.+z[ndxuse]**2-p**2)**2
        tmp = np.where(tmp < 0,0,tmp)
        lambdae[ndxuse] = (p**2*kap0+kap1 - 0.5*np.sqrt(tmp))/np.pi
        # eta_1
        etad[ndxuse] = 1./2./np.pi*(kap1+p**2*(p**2+2.*z[ndxuse]**2)*kap0- \
           (1.+5.*p**2+z[ndxuse]**2)/4.*np.sqrt((1.-x1[ndxuse])*(x2[ndxuse]-1.)))

    # Case 5, 6, 7 - the edge of planet lies at origin of star
    ocltor = np.where(z[notusedyet] == p)#, complement=notused3)
    t = np.where(z[notusedyet] == p)
    if np.size(ocltor) != 0:
        ndxuse = notusedyet[ocltor]
        if p < 0.5:
            # Case 5
            q=2.*p  # corrected typo in paper (2k -> 2p)
            Ek,Kk = ellke(q)
            # lambda_4
            lambdad[ndxuse] = 1./3.+2./9./np.pi*(4.*(2.*p**2-1.)*Ek+\
                                              (1.-4.*p**2)*Kk)
            # eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)
            lambdae[ndxuse] = p**2 # uniform disk
        elif p > 0.5:
            # Case 7
            q=0.5/p # corrected typo in paper (1/2k -> 1/2p)
            Ek,Kk = ellke(q)
            # lambda_3
            lambdad[ndxuse] = 1./3.+16.*p/9./np.pi*(2.*p**2-1.)*Ek-\
                              (32.*p**4-20.*p**2+3.)/9./np.pi/p*Kk
            # etad = eta_1 already
        else:
            # Case 6
            lambdad[ndxuse] = 1./3.-4./np.pi/9.
            etad[ndxuse] = 3./32.
        notused3 = np.where(z[notusedyet] != p)
        if np.size(notused3) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                      (lambdad+2./3.*(p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused3]

    # Case 2, Case 8 - ingress/egress (with limb darkening)
    inegress = np.where( ((z[notusedyet] > 0.5+abs(p-0.5)) & \
                       (z[notusedyet] < 1.+p))  | \
                      ( (p > 0.5) & (z[notusedyet] > abs(1.-p)) & \
                        (z[notusedyet] < p)) )#, complement=notused4)
    if np.size(inegress) != 0:

        ndxuse = notusedyet[inegress]
        q=np.sqrt((1.-x1[ndxuse])/(x2[ndxuse]-x1[ndxuse]))
        Ek,Kk = ellke(q)
        n=1./x1[ndxuse]-1.

        # lambda_1:
        lambdad[ndxuse]=2./9./np.pi/np.sqrt(x2[ndxuse]-x1[ndxuse])*\
                         (((1.-x2[ndxuse])*(2.*x2[ndxuse]+x1[ndxuse]-3.)-\
                           3.*x3[ndxuse]*(x2[ndxuse]-2.))*Kk+(x2[ndxuse]-\
                           x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-\
                          3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        notused4 = np.where( ( (z[notusedyet] <= 0.5+abs(p-0.5)) | \
                            (z[notusedyet] >= 1.+p) ) & ( (p <= 0.5) | \
                            (z[notusedyet] <= abs(1.-p)) | \
                            (z[notusedyet] >= p) ))
        if np.size(notused4) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*\
                                                     (p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused4]

    # Case 3, 4, 9, 10 - planet completely inside star
    if p < 1.:
        inside = np.where(z[notusedyet] <= (1.-p))#, complement=notused5)
        if np.size(inside) != 0:
            ndxuse = notusedyet[inside]

            ## eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)

            ## uniform disk
            lambdae[ndxuse] = p**2

            ## Case 4 - edge of planet hits edge of star
            edge = np.where(z[ndxuse] == 1.-p)#, complement=notused6)
            if np.size(edge[0]) != 0:
                ## lambda_5
                lambdad[ndxuse[edge]] = 2./3./pi*np.arccos(1.-2.*p)-\
                                      4./9./pi*np.sqrt(p*(1.-p))*(3.+2.*p-8.*p**2)
                if p > 0.5:
                    lambdad[ndxuse[edge]] -= 2./3.
                notused6 = np.where(z[ndxuse] != 1.-p)
                if np.size(notused6) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                              (lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused6[0]]

            ## Case 10 - origin of planet hits origin of star
            origin = np.where(z[ndxuse] == 0)#, complement=notused7)
            if np.size(origin) != 0:
                ## lambda_6
                lambdad[ndxuse[origin]] = -2./3.*(1.-p**2)**1.5
                notused7 = np.where(z[ndxuse] != 0)
                if np.size(notused7) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*\
                              (lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused7[0]]

            q=np.sqrt((x2[ndxuse]-x1[ndxuse])/(1.-x1[ndxuse]))
            n=x2[ndxuse]/x1[ndxuse]-1.
            Ek,Kk = ellke(q)

            ## Case 3, Case 9 - anynp.where in between
            ## lambda_2
            lambdad[ndxuse] = 2./9./np.pi/np.sqrt(1.-x1[ndxuse])*\
                              ((1.-5.*z[ndxuse]**2+p**2+x3[ndxuse]**2)*Kk+\
                               (1.-x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-\
                               3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        ## if there are still unused elements, there's a bug in the code
        ## (please report it)
        notused5 = np.where(z[notusedyet] > (1-p))
        if notused5[0].shape[0] != 0:
            print("ERROR: the following values of z didn't fit into a case:")
            return [-1,-1]

        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+\
                  u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]

def time2z(time, ipct, tknot, sma, orbperiod, ecc, tperi=None, epsilon=1e-5):
    '''
    G. ROUDIER: Time samples in [Days] to separation in [R*]
    '''
    if tperi is not None:
        ft0 = (tperi - tknot) % orbperiod
        ft0 /= orbperiod
        if ft0 > 0.5: ft0 += -1e0
        M0 = 2e0*np.pi*ft0
        E0 = solveme(M0, ecc, epsilon)
        realf = np.sqrt(1e0 - ecc)*np.cos(E0/2e0)
        imagf = np.sqrt(1e0 + ecc)*np.sin(E0/2e0)
        w = np.angle(np.complex(realf, imagf))
        if abs(ft0) < epsilon:
            w = np.pi/2e0
            tperi = tknot
            pass
        pass
    else:
        w = np.pi/2e0
        tperi = tknot
        pass
    ft = (time - tperi) % orbperiod
    ft /= orbperiod
    sft = np.copy(ft)
    sft[(sft > 0.5)] += -1e0
    M = 2e0*np.pi*ft
    E = solveme(M, ecc, epsilon)
    realf = np.sqrt(1. - ecc)*np.cos(E/2e0)
    imagf = np.sqrt(1. + ecc)*np.sin(E/2e0)
    f = []
    for r, i in zip(realf, imagf):
        cn = np.complex(r, i)
        f.append(2e0*np.angle(cn))
        pass
    f = np.array(f)
    r = sma*(1e0 - ecc**2)/(1e0 + ecc*np.cos(f))
    z = r*np.sqrt(1e0**2 - (np.sin(w+f)**2)*(np.sin(ipct*np.pi/180e0))**2)
    z[sft < 0] *= -1e0
    return z, sft

def solveme(M, e, eps):
    '''
    G. ROUDIER: Newton Raphson solver for true anomaly
    M is a numpy array
    '''
    E = np.copy(M)
    for i in np.arange(M.shape[0]):
        while abs(E[i] - e*np.sin(E[i]) - M[i]) > eps:
            num = E[i] - e*np.sin(E[i]) - M[i]
            den = 1. - e*np.cos(E[i])
            E[i] = E[i] - num/den
            pass
        pass
    return E

def transit(time, values):
    sep,phase = time2z(time, values['inc'], values['tmid'], values['ars'], values['per'], values['ecc'])
    model, _ = occultquad(abs(sep), values['u1'], values['u2'], values['rprs'])
    return model


def getPhase(curTime, pPeriod, tMid):
    phase = ((curTime - tMid) / pPeriod) % 1
    mask = phase >= 0.5
    phase[mask] = -1 * (1-phase[mask])
    return phase

# Creates binned data points for Final Light Curve
def time_bin(time, flux, dt):
    bins = int(np.floor((max(time) - min(time))/dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time)+i*dt)) & (time < (min(time)+(i+1)*dt))
        if mask.sum() > 0:
            bflux[i] = np.nanmean(flux[mask])
            btime[i] = np.nanmean(time[mask])
    zmask = (bflux==0) | (btime==0) | np.isnan(bflux) | np.isnan(btime)
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

    def __init__(self, time, data, dataerr, airmass, prior, bounds):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.airmass = airmass
        self.prior = prior
        self.bounds = bounds
        self.fit_nested()

    def fit_nested(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        bounddiff = np.diff(boundarray,1).reshape(-1)

        def loglike(pars):
            # chi-squared
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= self.prior['a1']*np.exp(self.prior['a2']*self.airmass)
            return -0.5 * np.sum( ((self.data-model)/self.dataerr)**2 )

        def prior_transform(upars):
            # transform unit cube to prior volume
            return (boundarray[:,0] + bounddiff*upars)

        dsampler = dynesty.DynamicNestedSampler(
            loglike, prior_transform,
            ndim=len(freekeys), bound='multi', sample='unif',
            maxiter_init=5000, dlogz_init=1, dlogz=0.05,
            maxiter_batch=100, maxbatch=10, nlive_batch=100
        )
        dsampler.run_nested()
        self.results = dsampler.results

        # alloc data for best fit + error
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        tests = [copy.deepcopy(self.prior) for i in range(6)]

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
            self.errors[freekeys[i]] = cov[i,i]**0.5
            tests[0][freekeys[i]] = mean[i]
            tests[1][freekeys[i]] = mean2[i]

            counts, bins = np.histogram(samples[:,i], bins=100, weights=weights)
            mi = np.argmax(counts)
            tests[5][freekeys[i]] = bins[mi] + 0.5*np.mean(np.diff(bins))

            # finds median and +- 2sigma, will vary from mode if non-gaussian
            self.quantiles[freekeys[i]] = dynesty.utils.quantile(self.results.samples[:,i], [0.025, 0.5, 0.975], weights=weights)
            tests[2][freekeys[i]] = self.quantiles[freekeys[i]][1]

        # find minimum near weighted mean
        mask = (samples[:,0] < self.parameters[freekeys[0]]+2*self.errors[freekeys[0]]) & (samples[:,0] > self.parameters[freekeys[0]]-2*self.errors[freekeys[0]])
        bi = np.argmin(self.weights[mask])

        for i in range(len(freekeys)):
            tests[3][freekeys[i]] = samples[mask][bi,i]
            tests[4][freekeys[i]] = np.average(samples[mask][:,i],weights=self.weights[mask],axis=0)

        # find best fit
        chis = []
        for i in range(len(tests)):
            lightcurve = transit(self.time, tests[i])
            airmass = tests[i]['a1']*np.exp(tests[i]['a2']*self.airmass)
            residuals = self.data - (lightcurve*airmass)
            chis.append( np.sum(residuals**2) )
            # chi2 of binned data in transit + small Oot
            # btime, br = time_bin(self.time, residuals)
            # blc = transit(btime, tests[i])
            # mask = blc < 1
            # duration = btime[mask].max() - btime[mask].min()
            # tmask = ((btime - tests[i]['tmid']) < duration) & ((btime - tests[i]['tmid']) > -1*duration)
            # chis.append( np.mean(br[tmask]**2) )

        mi = np.argmin(chis)
        self.parameters = copy.deepcopy(tests[mi])

        # final model
        self.phase = getPhase(self.time, self.parameters['per'], self.parameters['tmid'])
        self.transit = transit(self.time, self.parameters)
        self.airmass_model = self.parameters['a1']*np.exp(self.parameters['a2']*self.airmass)
        self.model = self.transit * self.airmass_model
        self.detrended = self.data / self.airmass_model
        self.residuals = self.data - self.model

    def plot_bestfit(self, nbins=10, phase=True):

        f = plt.figure( figsize=(12,7) )
        #f.subplots_adjust(top=0.94,bottom=0.08,left=0.07,right=0.96)
        ax_lc = plt.subplot2grid( (4,5), (0,0), colspan=5,rowspan=3 )
        ax_res = plt.subplot2grid( (4,5), (3,0), colspan=5, rowspan=1, sharex=ax_lc )
        axs = [ax_lc, ax_res]

        dt = (max(self.time) - min(self.time))/nbins

        if phase == True:
            axs[0].errorbar(self.phase, self.detrended, yerr=self.dataerr, ls='none', marker='o', color='gray', markersize=5, zorder=1)
            axs[0].plot(self.phase, self.transit, 'r-', zorder=2)
            axs[0].set_ylabel("Relative Flux")
            axs[0].grid(True,ls='--')

            axs[1].plot(self.phase, 1e6*self.residuals/np.median(self.data), marker='o', color='gray', markersize=5, ls='none')
            axs[1].plot(self.phase, np.zeros(len(self.phase)), 'r-', lw=2, alpha=1, zorder=100)   ###maybe
            axs[1].set_xlabel("Phase")
            axs[1].set_ylabel("Residuals [PPM]")

            ax_res.errorbar(binner(self.phase, len(self.residuals) // 10), binner(1e6*self.residuals/np.median(self.data), len(self.residuals) // 10),
                            yerr=binner(1e6*self.residuals/np.median(self.data), len(self.residuals) // 10, self.dataerr / self.airmass_model)[1],
                            fmt='s', mfc='b', mec='b', ecolor='b', zorder=10)
            ax_lc.errorbar(binner(self.phase, len(self.phase) // 10),
                            binner(self.detrended / self.airmass_model, len(self.phase) // 10),
                            yerr=binner(1e6*self.residuals/np.median(self.data), len(self.residuals) // 10, self.dataerr / self.airmass_model)[1],
                            fmt='s', mfc='b', mec='b', ecolor='b', zorder=10)

        #    binnedPhase,
        ##    binnedPhase2,
        #    binnedResids = binner(self.phase, 1e6*self.residuals/np.median(self.data), dt/self.parameters['per'])
            #binnedPhase, binnedFlux = time_bin(self.phase, self.detrended, dt/self.parameters['per'])
            #binnedPhase2, binnedResids = time_bin(self.phase, 1e6*self.residuals/np.median(self.data), dt/self.parameters['per'])
        #    axs[0].plot(binnedPhase, binnedFlux, marker='s', color='blue', markersize=6, ls='none')
        #    axs[1].plot(binnedPhase2, binnedResids, marker='s', color='blue', markersize=6, ls='none')

        #COPY INTO BELOW!

        else:
            axs[0].errorbar(self.time, self.detrended, yerr=self.dataerr, ls='none', marker='o', color='gray', markersize=5, zorder=1)
            axs[0].plot(self.time, self.transit, 'r-', zorder=2)
            axs[0].set_ylabel("Relative Flux")
            axs[0].grid(True,ls='--')

            axs[1].plot(self.time, 1e6*self.residuals/np.median(self.data), marker='o', color='gray', markersize=5, ls='none')
            axs[1].plot(phase, np.zeros(len(self.phase)), 'r-', lw=2, alpha=1, zorder=100)   ###maybe
            axs[1].set_xlabel("Time [day]")
            axs[1].set_ylabel("Residuals [PPM]")

            ax_res.errorbar(binner(self.time, len(self.residuals) // 10), binner(self.residuals, len(self.residuals) // 10),
                            yerr=binner(self.residuals, len(self.residuals) // 10, self.dataerr / self.airmass_model)[1],
                            fmt='s', mfc='b', mec='b', ecolor='b', zorder=10)
            ax_lc.errorbar(binner(self.time, len(self.time) // 10),
                            binner(self.detrended / self.airmass_model, len(self.phase) // 10),
                            yerr=binner(self.residuals, len(self.residuals) // 10, self.dataerr / self.airmass_model)[1],
                            fmt='s', mfc='b', mec='b', ecolor='b', zorder=10)
#
        #    binnedTime,
#            binnedFlux = binner(self.time, self.detrended, dt)
#            binnedTime2,
#            binnedResids = binner(self.time, 1e6*self.residuals/np.median(self.data), dt)
            #binnedTime, binnedFlux = time_bin(self.time, self.detrended, dt)
            #binnedTime2, binnedResids = time_bin(self.time, 1e6*self.residuals/np.median(self.data), dt)
#            axs[0].plot(binnedTime, binnedFlux, marker='s', color='blue', markersize=6, ls='none')
#            axs[1].plot(binnedTime2, binnedResids, marker='s', color='blue', markersize=6, ls='none')

        plt.tight_layout()

        return f,axs

if __name__ == "__main__":

    prior = {
        'rprs':0.03,        # Rp/Rs
        'ars':14.25,        # a/Rs
        'per':3.336817,     # Period [day]
        'inc':87.5,        # Inclination [deg]
        'u1': 0.3, 'u2': 0.1, # limb darkening (linear, quadratic)
        'ecc':0,            # Eccentricity
        'omega':0,          # Arg of periastron
        'tmid':0.75,         # time of mid transit [day],

        'a1':50,            # airmass coeffcients
        'a2':0.25
    }

    time = np.linspace(0.65,0.85,150) # [day]

    # simulate extinction from airmass
    stime = time-time[0]
    alt = 90* np.cos(4*stime-np.pi/6)
    airmass = 1./np.cos( np.deg2rad(90-alt))

    # GENERATE NOISY DATA
    data = transit(time, prior)*prior['a1']*np.exp(prior['a2']*airmass)
    data += np.random.normal(0, prior['a1']*250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time))


    # add bounds for free parameters only
    mybounds = {
        'rprs':[0,0.1],
        'tmid':[min(time),max(time)],
        'ars':[13,15],

        'a1':[25,75],
        'a2':[0,0.3]
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds)

    for k in myfit.bounds.keys():
        print("{:.6f} +- {}".format( myfit.parameters[k], myfit.errors[k]))

    fig,axs = myfit.plot_bestfit()

    # triangle plot
    fig,axs = dynesty.plotting.cornerplot(myfit.results, labels=list(mybounds.keys()), quantiles_2d=[0.4,0.85], smooth=0.015, show_titles=True,use_math_text=True, title_fmt='.2e',hist2d_kwargs={'alpha':1,'zorder':2,'fill_contours':False})
    dynesty.plotting.cornerpoints(myfit.results, labels=list(mybounds.keys()), fig=[fig,axs[1:,:-1]],plot_kwargs={'alpha':0.1,'zorder':1,} )
    plt.tight_layout()
    plt.show()
    #plt.savefig("temp.png")
