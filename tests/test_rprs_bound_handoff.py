"""Search preparation must retain support already opened by the LM scout."""
import copy
from types import SimpleNamespace

import numpy as np
import pytest

import exotic.exotic as ex


@pytest.fixture
def wasp_window(monkeypatch):
    monkeypatch.setattr(ex,'RPRS_RANGE_RESTRICTION_ENABLED',True)
    monkeypatch.setattr(ex,'RPRS_RANGE_RESTRICTION_PERCENTAGE',10.)
    monkeypatch.setattr(ex,'RPRS_SEARCH_BOUND_MAX',.5)
    center=.11064687659396447
    return dict(rprs=center,rprs_data_uncertainty=(.15210345102939868-center)/3.)


def test_real_wasp_scout_range_survives_data_uncertainty_widening(wasp_window):
    bounds=dict(rprs=[.05053160490809812,.21233694594729358],ars=[3.232143276808586,26.895530418528125])
    original=copy.deepcopy(bounds)
    result=ex.widen_rprs_bounds_to_data_uncertainty_window(bounds,wasp_window)
    assert result==original
    assert bounds==original
    assert result['rprs'] is not bounds['rprs']


@pytest.mark.parametrize('interval',[[.09,.13],[.04,.14],[.09,.24],[.04,.24]])
def test_data_window_expands_each_side_without_shrinking_existing_support(wasp_window,interval):
    window=ex.rprs_prior_centered_bounds(wasp_window)
    result=ex.widen_rprs_bounds_to_data_uncertainty_window({'rprs':interval},wasp_window)
    assert result['rprs']==[min(interval[0],window[0]),max(interval[1],window[1])]


def test_sampler_handoff_keeps_scout_expansion_but_honors_explicit_restriction(wasp_window):
    bounds={'rprs':[.05053160490809812,.21233694594729358]}
    expanded=ex.apply_configured_prior_search_restrictions(bounds,wasp_window,allow_rprs_expansion=True)
    assert expanded==bounds
    restricted=ex.apply_configured_prior_search_restrictions(bounds,wasp_window,allow_rprs_expansion=False)
    np.testing.assert_allclose(restricted['rprs'],ex.rprs_prior_centered_bounds(wasp_window),rtol=0,atol=0)


@pytest.mark.parametrize('malformed',[None,[np.nan,np.nan],['bad']])
def test_existing_data_window_can_still_supply_missing_numeric_endpoints(wasp_window,malformed):
    result=ex.widen_rprs_bounds_to_data_uncertainty_window({'rprs':malformed},wasp_window)
    assert result['rprs']==ex.rprs_prior_centered_bounds(wasp_window)


def test_full_scout_to_sampler_handoff_retains_recorded_range(monkeypatch,wasp_window):
    opened={'rprs':[.05053160490809812,.21233694594729358]}
    def scout(times,flux,unc,airmass,prior,bounds,**kwargs):
        return dict(prior),copy.deepcopy(opened),dict(
            enabled=True,applied=True,adjusted=True,expanded_keys=['rprs'],bounds=opened)
    captured=[]
    def fitter(times,flux,unc,airmass,prior,bounds,**kwargs):
        captured.append(copy.deepcopy(bounds))
        return SimpleNamespace(parameters=dict(prior),errors={},bounds=bounds,results=None)
    monkeypatch.setattr(ex,'prepare_ultranest_bounds_with_lm_boundary_scout',scout)
    monkeypatch.setattr(ex,'lc_fitter',fitter)
    t=np.linspace(-.1,.1,20)
    fit=ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        t,np.ones(20),np.full(20,.01),np.ones(20),
        dict(wasp_window,per=3.12,tmid=0.),{'rprs':[.09,.13]},
        search_restriction_prior=wasp_window,use_lm_boundary_scout=True,
        use_prior_rprs_when_posterior_pinned=False)
    assert captured==[opened]
    assert fit.bounds==fit.lm_boundary_scout_final_bounds==fit.posterior_refit_final_bounds
