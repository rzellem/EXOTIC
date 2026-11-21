"""Native kernel parity and explicit native-execution failures."""
import numpy as np
import pytest

from exotic.api import elca, limbdarkfast_cpp
from exotic.api.limbdarkfast import hybrid_flux


@pytest.fixture
def cpp():
    if not limbdarkfast_cpp.library_path().is_file():
        pytest.skip('Build the optional native library with scripts/build_hybrid4_cpp.py')
    return limbdarkfast_cpp.hybrid4_flux


@pytest.mark.parametrize('radius', [0.0, 0.01, 0.1, 0.5, 0.9999, 1.0, 1.0001, 2.0])
def test_cpp_matches_numpy_across_contacts_and_switches(cpp, radius):
    weights = [0.3709047608728925, 0.00711438176464726, 0.3162293789763985, 0.01888876666845794]
    contacts = np.array([0.0, radius, abs(1.0-radius), 1.0+radius])
    distances = np.concatenate([np.linspace(0, 1+radius, 400),
                                np.maximum(0, (contacts[:,None] + np.array([-1e-13, 0, 1e-13]))).ravel()])
    np.testing.assert_allclose(cpp(distances, radius, weights),
                               hybrid_flux(distances, radius, weights, 'hybrid4'), rtol=0, atol=1e-12)


def test_cpp_broadcasts_and_accepts_noncontiguous_inputs_and_negative_weights(cpp):
    z = np.linspace(0, 2, 60)[::2].reshape(5, 6)
    radius = np.linspace(0.01, 0.2, 6)
    weights = [0.5, -0.2, 0.3, 0.6]
    np.testing.assert_allclose(cpp(z, radius, weights), hybrid_flux(z, radius, weights, 'hybrid4'),
                               rtol=0, atol=1e-12)
    assert cpp(0.3, 0.1, weights).shape == ()
    assert cpp([], 0.1, weights).shape == (0,)


def test_cpp_reports_native_domain_error_without_crossing_the_c_abi(cpp):
    with pytest.raises(ValueError, match='r >= 0'):
        cpp([0.0, 0.5], -0.1, [0.3,0.2,0.1,0.05])


def test_missing_cpp_library_is_an_explicit_error(tmp_path, monkeypatch):
    limbdarkfast_cpp._library.cache_clear()
    with monkeypatch.context() as patch:
        patch.setattr(limbdarkfast_cpp, 'library_path', lambda: tmp_path/'absent.dll')
        # Exercise the load error independently of the source-build tests. A
        # synthetic missing-library path must not start a real compiler here.
        patch.setattr(limbdarkfast_cpp, 'ensure_hybrid4_native', lambda **kwargs: None)
        with pytest.raises(RuntimeError, match='scripts/build_hybrid4_cpp.py'):
            limbdarkfast_cpp.hybrid4_flux([0.2],0.1,[0.3,0.2,0.1,0.05])
    limbdarkfast_cpp._library.cache_clear()


def test_cpp_model_retains_orbit_and_variable_exposure_integration(cpp):
    prior = dict(rprs=0.1, per=3.0, ars=10.0, inc=88.0, ecc=0.4, omega=120.0,
                 tmid=2460000.0, u0=0.0,u1=0.0,u2=0.0,u3=0.0,
                 _limb_darkening_law='hybrid4', h0=0.3,h1=0.2,h2=0.1,h3=0.05)
    times=2460000.0+np.linspace(-0.07,0.07,400)
    prior['_exposure_time_days']=np.linspace(10,120,len(times))/86400.0
    numpy_flux=elca.transit(times,{**prior,'_hybrid4_backend':'numpy'})
    cpp_flux=elca.transit(times,{**prior,'_hybrid4_backend':'cpp'})
    np.testing.assert_allclose(cpp_flux,numpy_flux,rtol=0,atol=1e-12)


def test_missing_cpp_library_fails_before_sampling(monkeypatch):
    def fail():
        raise RuntimeError('native library unavailable')
    monkeypatch.setattr(limbdarkfast_cpp, '_library', fail)
    p={'_limb_darkening_law':'hybrid4','_hybrid4_backend':'cpp'}
    with pytest.raises(RuntimeError, match='native library unavailable'):
        elca.lc_fitter(np.arange(3),np.ones(3),np.ones(3),np.ones(3),p,{'rprs':[0,1]},mode='ns')
