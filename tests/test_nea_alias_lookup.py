import requests

from exotic.api.nea import NEA_ALIAS_LOOKUP_URL, resolve_planet_alias


class DummyResponse:
    def __init__(self, payload, status_ok=True):
        self._payload = payload
        self._status_ok = status_ok

    def raise_for_status(self):
        if not self._status_ok:
            raise requests.exceptions.HTTPError('503')

    def json(self):
        return self._payload


def alias_payload(requested, resolved, status, planets):
    return {
        'manifest': {
            'requested_name': requested,
            'resolved_name': resolved,
            'lookup_status': status,
        },
        'system': {
            'objects': {
                'planet_set': {
                    'item_count': len(planets),
                    'planets': {name: {} for name in planets},
                }
            }
        },
    }


def make_getter(payload, status_ok=True, calls=None):
    def getter(url, params=None, timeout=None):
        if calls is not None:
            calls.append((url, params, timeout))
        return DummyResponse(payload, status_ok=status_ok)
    return getter


def test_planet_alias_resolves_to_default_planet_name():
    calls = []
    payload = alias_payload('HAT-P-10 b', 'WASP-11 b', 'OK', ['WASP-11 b'])
    assert resolve_planet_alias('HAT-P-10 b', getter=make_getter(payload, calls=calls)) == 'WASP-11 b'
    assert calls[0][0] == NEA_ALIAS_LOOKUP_URL
    assert calls[0][1] == {'objname': 'HAT-P-10 b'}


def test_star_alias_resolves_when_system_has_one_planet():
    payload = alias_payload('HAT-P-10', 'WASP-11', 'OK', ['WASP-11 b'])
    assert resolve_planet_alias('HAT-P-10', getter=make_getter(payload)) == 'WASP-11 b'


def test_star_alias_is_ambiguous_in_multi_planet_system():
    payload = alias_payload('TRAPPIST-1', 'TRAPPIST-1', 'OK',
                            ['TRAPPIST-1 b', 'TRAPPIST-1 c', 'TRAPPIST-1 d'])
    assert resolve_planet_alias('TRAPPIST-1', getter=make_getter(payload)) is None


def test_unknown_name_returns_none():
    payload = alias_payload('Nonsense-99 b', None, 'System Not Found', [])
    assert resolve_planet_alias('Nonsense-99 b', getter=make_getter(payload)) is None


def test_service_failure_returns_none():
    payload = alias_payload('HAT-P-10 b', 'WASP-11 b', 'OK', ['WASP-11 b'])
    assert resolve_planet_alias('HAT-P-10 b', getter=make_getter(payload, status_ok=False)) is None

    def broken_getter(url, params=None, timeout=None):
        raise requests.exceptions.ConnectionError('offline')

    assert resolve_planet_alias('HAT-P-10 b', getter=broken_getter) is None


def test_malformed_payload_returns_none():
    assert resolve_planet_alias('HAT-P-10 b', getter=make_getter({'manifest': 'not a dict'})) is None
    assert resolve_planet_alias('HAT-P-10 b', getter=make_getter([])) is None
