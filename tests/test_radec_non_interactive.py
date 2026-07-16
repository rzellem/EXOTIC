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
