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

import logging
import numpy as np
import lmfit as lm
import matplotlib.pyplot as plt
import ldtk
from ldtk import LDPSetCreator, BoxcarFilter
from ldtk.ldmodel import LinearModel, QuadraticModel, NonlinearModel

log = logging.getLogger(__name__)


class LDPSet(ldtk.LDPSet):
    """
    A. NIESSNER: INLINE HACK TO ldtk.LDPSet
    """
    @staticmethod
    def is_mime(): return True

    @property
    def profile_mu(self): return self._mu
    pass


setattr(ldtk, 'LDPSet', LDPSet)
setattr(ldtk.ldtk, 'LDPSet', LDPSet)


def createldgrid(minmu, maxmu, orbp,
                 ldmodel='nonlinear', phoenixmin=1e-1,
                 segmentation=int(10), verbose=False):
    """
    G. ROUDIER: Wrapper around LDTK downloading tools
    LDTK: Parviainen et al. https://github.com/hpparvi/ldtk
    """
    tstar = orbp['T*']
    terr = np.sqrt(abs(orbp['T*_uperr']*orbp['T*_lowerr']))
    fehstar = orbp['FEH*']
    feherr = np.sqrt(abs(orbp['FEH*_uperr']*orbp['FEH*_lowerr']))
    loggstar = orbp['LOGG*']
    loggerr = np.sqrt(abs(orbp['LOGG*_uperr']*orbp['LOGG*_lowerr']))
    log.warning(f">-- Temperature: {tstar} +/- {terr}")
    log.warning(f">-- Metallicity: {fehstar} +/- {feherr}")
    log.warning(f">-- Surface Gravity: {loggstar} +/- {loggerr}")
    niter = int(len(minmu)/segmentation) + 1
    allcl = None
    allel = None
    out = {}
    avmu = [np.mean([mm, xm]) for mm, xm in zip(minmu, maxmu)]
    for i in np.arange(niter):
        loweri = i*segmentation
        upperi = (i+1)*segmentation
        if i == (niter-1):
            upperi = len(avmu)
        munm = 1e3*np.array(avmu[loweri:upperi])
        munmmin = 1e3*np.array(minmu[loweri:upperi])
        munmmax = 1e3*np.array(maxmu[loweri:upperi])
        filters = [BoxcarFilter(str(mue), mun, mux)
                   for mue, mun, mux in zip(munm, munmmin, munmmax)]
        sc = LDPSetCreator(teff=(tstar, terr), logg=(loggstar, loggerr),
                           z=(fehstar, feherr), filters=filters)
        ps = sc.create_profiles(nsamples=int(1e4))
        itpfail = False
        for testprof in ps.profile_averages:
            if np.all(~np.isfinite(testprof)):
                itpfail = True
            pass
        nfail = 1e0
        while itpfail:
            nfail *= 2
            sc = LDPSetCreator(teff=(tstar, nfail*terr), logg=(loggstar, loggerr),
                               z=(fehstar, feherr), filters=filters)
            ps = sc.create_profiles(nsamples=int(1e4))
            itpfail = False
            for testprof in ps.profile_averages:
                if np.all(~np.isfinite(testprof)):
                    itpfail = True
                pass
            pass
        cl, el = ldx(ps.profile_mu, ps.profile_averages, ps.profile_uncertainties,
                     mumin=phoenixmin, debug=verbose, model=ldmodel)
        if allcl is None:
            allcl = cl
        else:
            allcl = np.concatenate((allcl, cl), axis=0)
        if allel is None:
            allel = el
        else:
            allel = np.concatenate((allel, el), axis=0)
        pass
    allel[allel > 1.] = 0.
    allel[~np.isfinite(allel)] = 0.
    out['MU'] = avmu
    out['LD'] = allcl.T
    out['ERR'] = allel.T
    for i in range(0, len(allcl.T)):
        log.warning(f">-- LD{int(i)}: {float(allcl.T[i])} +/- {float(allel.T[i])}")
        pass
    return out


def ldx(psmu, psmean, psstd, mumin=1e-1, debug=False, model='nonlinear'):
    """
    G. ROUDIER: Limb darkening coefficient retrievial on PHOENIX GRID models,
    OPTIONAL mumin set up on HAT-P-41 stellar class
    """
    mup = np.array(psmu).copy()
    prfs = np.array(psmean).copy()
    sprfs = np.array(psstd).copy()
    nwave = prfs.shape[0]
    select = (mup > mumin)
    fitmup = mup[select]
    fitprfs = prfs[:, select]
    fitsprfs = sprfs[:, select]
    cl = []
    el = []
    params = lm.Parameters()
    params.add('gamma1', value=1e-1)
    params.add('gamma2', value=5e-1)
    params.add('gamma3', value=1e-1)
    params.add('gamma4', expr='1 - gamma1 - gamma2 - gamma3')

    if debug:
        plt.figure()
    for iwave in np.arange(nwave):
        select = fitsprfs[iwave] == 0e0
        if True in select:
            fitsprfs[iwave][select] = 1e-10
        if model == 'linear':
            params['gamma1'].value = 0
            params['gamma1'].vary = False
            params['gamma3'].value = 0
            params['gamma3'].vary = False
            params['gamma4'].value = 0
            params['gamma4'].vary = False
            out = lm.minimize(ln_ldx, params,
                              args=(fitmup, fitprfs[iwave], fitsprfs[iwave]))
            cl.append([out.params['gamma1'].value])
            el.append([out.params['gamma1'].stderr])
            pass
        if model == 'quadratic':
            params['gamma1'].value = 0
            params['gamma1'].vary = False
            params['gamma3'].value = 0
            params['gamma3'].vary = False
            out = lm.minimize(qd_ldx, params,
                              args=(mup, prfs[iwave], sprfs[iwave]))
            cl.append([out.params['gamma2'].value, out.params['gamma4'].value])
            el.append([out.params['gamma2'].stderr, out.params['gamma4'].stderr])
            pass
        if model == 'nonlinear':
            out = lm.minimize(nl_ldx, params,
                              args=(fitmup, fitprfs[iwave], fitsprfs[iwave]))
            cl.append([out.params['gamma1'].value, out.params['gamma2'].value,
                       out.params['gamma3'].value, out.params['gamma4'].value])
            el.append([out.params['gamma1'].stderr, out.params['gamma2'].stderr,
                       out.params['gamma3'].stderr, out.params['gamma4'].stderr])
            pass
        if debug:
            plt.plot(mup, prfs[iwave], 'k^')
            plt.errorbar(fitmup, fitprfs[iwave], yerr=fitsprfs[iwave], ls='None')
            if model == 'linear':
                plt.plot(fitmup, ln_ldx(out.params, fitmup))
            if model == 'quadratic':
                plt.plot(fitmup, qd_ldx(out.params, fitmup))
            if model == 'nonlinear':
                plt.plot(fitmup, nl_ldx(out.params, fitmup))
            pass
        pass
    if debug:
        plt.ylabel('$I(\\mu)$')
        plt.xlabel('$\\mu$')
        plt.title(model)
        plt.show()
        pass
    return np.array(cl), np.array(el)


def ln_ldx(params, x, data=None, weights=None):
    """
    G. ROUDIER: Linear law
    """
    gamma1 = params['gamma1'].value
    model = LinearModel.evaluate(x, np.array(gamma1))
    if data is None:
        return model
    if weights is None:
        return data - model
    return (data - model)/weights


def qd_ldx(params, x, data=None, weights=None):
    """
    G. ROUDIER: Quadratic law
    """
    gamma1 = params['gamma2'].value
    gamma2 = params['gamma4'].value
    model = QuadraticModel.evaluate(x, np.array((gamma1, gamma2)))
    if data is None:
        return model
    if weights is None:
        return data - model
    return (data - model)/weights


def nl_ldx(params, x, data=None, weights=None):
    """
    G. ROUDIER: Non Linear law
    """
    gamma1 = params['gamma1'].value
    gamma2 = params['gamma2'].value
    gamma3 = params['gamma3'].value
    gamma4 = params['gamma4'].value
    model = NonlinearModel.evaluate(x, np.array((gamma1, gamma2, gamma3, gamma4)))
    if data is None:
        return model
    if weights is None:
        return data - model
    return (data - model)/weights

