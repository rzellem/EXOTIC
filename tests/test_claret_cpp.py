"""Native precision=3 Claret model against the local Python reference."""
import numpy as np
import pytest

from exotic.api import elca


def parameters():
    return dict(rprs=.1,per=3.2,ars=9.,ecc=0.,inc=89.,omega=54.,tmid=2460000.,
                u0=1.772483643097537,u1=-2.9680582066134003,
                u2=3.4417486511535746,u3=-1.2461740876377116)


@pytest.mark.parametrize('ecc',[0.,.018,.12,.7])
@pytest.mark.parametrize('omega',[-90.,0.,54.,90.,270.])
@pytest.mark.parametrize('exposure',[None,15.,120.])
def test_claret_orbit_and_exposure_matches_reference(ecc,omega,exposure):
    p=parameters();p.update(ecc=ecc,omega=omega)
    t=p['tmid']+np.r_[np.linspace(-.2,.2,401),np.linspace(-7.36,7.36,181)]
    if exposure is not None:
        p['_exposure_time_days']=np.linspace(exposure/2,exposure*2,len(t))/86400.
    old=elca.transit(t,{**p,'_claret_backend':'python'})
    new=elca.transit(t,{**p,'_claret_backend':'cpp'})
    np.testing.assert_allclose(new,old,rtol=0,atol=2.e-12,equal_nan=True)


@pytest.mark.parametrize('radius',[0.,.001,.1,.5,.9999,1.,1.0001,2.])
def test_claret_contacts_and_grazing_transits(radius):
    p=parameters();p.update(rprs=radius,tmid=0.,inc=90.)
    z=np.r_[np.linspace(0,1+radius,500),abs(1-radius),radius,1+radius,np.sqrt(1+radius*radius)]
    t=np.arcsin(z/p['ars'])*p['per']/(2*np.pi)
    old=elca.transit(t,{**p,'_claret_backend':'python'})
    new=elca.transit(t,{**p,'_claret_backend':'cpp'})
    np.testing.assert_allclose(new[:-1],old[:-1],rtol=0,atol=2.e-12,equal_nan=True)
    # The last point is the double-root tangent z**2 - radius**2 == 1.
    # Its sector endpoint contains sqrt(radius**2 - (z*sin(theta))**2).
    # An O(eps) libm/NumPy difference inside that cancellation becomes
    # O(sqrt(eps)) in the endpoint; Ubuntu differs by 8.2e-10 in flux.
    # Allow that roundoff only at this tangent, keeping every other contact
    # and the ordinary orbit/exposure comparisons at their strict tolerance.
    np.testing.assert_allclose(new[-1:],old[-1:],rtol=0,
                               atol=np.sqrt(np.finfo(float).eps),equal_nan=True)


def test_claret_controls_shapes_and_invalid_exposure_preserve_reference():
    p=parameters()
    for t in [np.asarray(p['tmid']),np.array([]),p['tmid']+np.linspace(-.2,.2,60).reshape(6,10)[:,::2]]:
        for exposure in [None,30./86400.,[1.,2.],'invalid']:
            p['_exposure_time_days']=exposure
            # The reference PyLightcurve orbit expects a one-dimensional time
            # vector, while the public adapter preserves the caller's shape.
            expected=elca.transit(np.asarray(t).reshape(-1),{**p,'_claret_backend':'python'}).reshape(t.shape)
            actual=elca.transit(t,p)
            np.testing.assert_allclose(actual,expected,rtol=0,atol=2.e-12,equal_nan=True)


def test_native_model_does_not_call_python_claret(monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('native model used PyLightcurve')
    monkeypatch.setattr(elca,'pytransit',forbidden)
    actual=elca.transit(parameters()['tmid']+np.linspace(-.1,.1,75),parameters())
    assert np.isfinite(actual).all()
    assert actual.min()<1.


@pytest.mark.parametrize('mode',['lm','ns'])
def test_claret_build_failure_precedes_optimizer_or_sampler(monkeypatch,mode):
    from exotic.api import limbdarkfast_cpp
    def fail():
        raise RuntimeError('native build unavailable')
    monkeypatch.setattr(limbdarkfast_cpp,'_library',fail)
    t=parameters()['tmid']+np.linspace(-.1,.1,30)
    with pytest.raises(RuntimeError,match='native build unavailable'):
        elca.lc_fitter(t,np.ones(30),np.full(30,.001),np.ones(30),parameters(),
                       {'rprs':[.05,.2]},mode=mode,verbose=False)
