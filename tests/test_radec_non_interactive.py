import pytest

from exotic import exotic as exotic_module


def test_invalid_target_coordinates_use_nasa_archive_fallback_without_prompt(monkeypatch):
    messages = []
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("non-interactive coordinate resolution must not prompt"),
    )
    monkeypatch.setattr(
        exotic_module,
        'log_info',
        lambda message, **kwargs: messages.append((message, kwargs)),
    )

    ra, dec = exotic_module.radec_hours_to_degree(
        'not-an-ra',
        '+20:00:00',
        non_interactive_run=True,
        archive_ra=123.456,
        archive_dec=-45.678,
        target_name='Example b',
    )

    assert ra == pytest.approx(123.456)
    assert dec == pytest.approx(-45.678)
    assert len(messages) == 1
    assert "Using NASA Exoplanet Archive coordinates" in messages[0][0]
    assert messages[0][1] == {'warn': True}


def test_invalid_target_and_archive_coordinates_abort_without_prompt(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("non-interactive coordinate resolution must not prompt"),
    )

    with pytest.raises(
        ValueError,
        match=(
            r"Non-interactive run cancelled for target Example b: .*"
            r"NASA Exoplanet Archive coordinates .* are also unusable"
        ),
    ):
        exotic_module.radec_hours_to_degree(
            'not-an-ra',
            '+20:00:00',
            non_interactive_run=True,
            archive_ra='also-not-an-ra',
            archive_dec='also-not-a-dec',
            target_name='Example b',
        )


def test_invalid_target_coordinates_abort_when_archive_coordinates_unavailable(monkeypatch):
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("non-interactive coordinate resolution must not prompt"),
    )

    with pytest.raises(
        ValueError,
        match=(
            r"Non-interactive run cancelled for target Example b: .*"
            r"NASA Exoplanet Archive coordinates are unavailable"
        ),
    ):
        exotic_module.radec_hours_to_degree(
            'not-an-ra',
            '+20:00:00',
            non_interactive_run=True,
            target_name='Example b',
        )


def test_convert_jd_to_bjd_coerces_sexagesimal_strings():
    """Regression: with -ov -nea the main flow used to deliver sexagesimal
    RA/Dec STRINGS to convert_jd_to_bjd. barycorrpy then raised TypeError and
    the astropy fallback parsed "19:27:06.50" as 19.45 DEGREES instead of
    19h27m = 291.78 degrees, silently shifting BJD_TDB by hundreds of seconds
    (-538 s for these CoRoT-2 coordinates). String and numeric inputs must
    produce the identical, correct conversion."""
    site = {'lat': 31.68, 'long': -110.88, 'elev': 1268.0}
    epochs = [2461222.65]

    from_strings = exotic_module.convert_jd_to_bjd(
        epochs, {'ra': '19:27:06.50', 'dec': '+01:23:01.5'}, site)
    from_degrees = exotic_module.convert_jd_to_bjd(
        epochs, {'ra': 291.777083, 'dec': 1.383750}, site)

    import numpy as np
    a = float(np.atleast_1d(from_strings)[0])
    b = float(np.atleast_1d(from_degrees)[0])
    assert a == pytest.approx(b, abs=1e-8)
    # The correct value is ~ +522.5 s after the input JD_UTC (TDB-UTC ~ +69 s
    # plus ~ +453 s light travel time for this geometry); the broken fallback
    # produced a value ~16 s BEFORE it. Guard the sign and scale.
    assert (a - epochs[0]) * 86400 == pytest.approx(522.5, abs=2.0)
