import csv

from exotic.plate_status import PlateStatus


def test_out_of_frame_warning_reports_only_fits_basename():
    messages = []
    status = PlateStatus(lambda message, **kwargs: messages.append(message))
    status.setCurrentFilename('/mnt/data/session/TIC 13510052901-R-20240403-030-063313_out.fits')

    status.outOfFrameWarning(11)

    assert messages[0] == (
        'Comparison star #11 is beyond the edge of file '
        'TIC 13510052901-R-20240403-030-063313_out.fits'
    )
    assert 'repeated frame-level star warnings are aggregated' in messages[1]


def test_out_of_frame_warning_preserves_fits_fz_basename():
    messages = []
    status = PlateStatus(lambda message, **kwargs: messages.append(message))
    status.setCurrentFilename(r'C:\data\session\compressed-frame.fits.fz')

    status.outOfFrameWarning(0)

    assert messages[0] == 'Target star is beyond the edge of file compressed-frame.fits.fz'


def test_tracked_vsx_label_is_used_for_every_star_warning_type():
    messages = []
    status = PlateStatus(lambda message, **kwargs: messages.append(message))
    status.setComparisonStarLabels({18: 'Tracked VSX variable DI Her'})

    status.setCurrentFilename('/mnt/data/frame-1.fits.fz')
    status.outOfFrameWarning(18)
    status.setCurrentFilename('/mnt/data/frame-2.fits.fz')
    status.lowFluxAmplitudeWarning(18, 123.4, 234.5)
    status.setCurrentFilename('/mnt/data/frame-3.fits.fz')
    status.overexposedWarning(18, 124.4, 235.5, 58981.5)
    status.setCurrentFilename('/mnt/data/frame-4.fits.fz')
    status.skyBackgroundWarning(18, 125.4, 236.5)

    warning_messages = [message for message in messages if 'file frame-' in message]
    assert len(warning_messages) == 4
    assert all('Tracked VSX variable DI Her' in message for message in warning_messages)
    assert all('Comparison star' not in message for message in warning_messages)
    assert all('/mnt/data/' not in message for message in warning_messages)


def test_repeated_frame_warnings_are_aggregated_with_progress_and_summary():
    messages = []
    status = PlateStatus(lambda message, **kwargs: messages.append(message))

    for frame_index in range(205):
        status.setCurrentFilename(f'/mnt/data/frame-{frame_index:04d}.fits')
        status.overexposedWarning(1, 100.0, 200.0, 50000.0)

    status.logAggregatedWarningSummary()

    assert sum(
        'Comparison star #1 is overexposed in file' in message
        for message in messages
    ) == 1
    assert any(
        'Comparison star #1: overexposed in 100 frame(s) so far.' in message
        for message in messages
    )
    assert any(
        'Comparison star #1: overexposed in 200 frame(s) so far.' in message
        for message in messages
    )
    assert messages[-1] == '>-- Comparison star #1: overexposed in 205 frame(s).'
    assert sum(
        'overexposed_comp1' in frame_status
        for frame_status in status.statusByFilename.values()
    ) == 205

    message_count = len(messages)
    status.logAggregatedWarningSummary()
    assert len(messages) == message_count


def test_aggregated_warnings_preserve_exact_per_frame_csv_flags(tmp_path):
    messages = []
    filenames = [f'/mnt/data/frame-{frame_index:04d}.fits' for frame_index in range(3)]
    status = PlateStatus(lambda message, **kwargs: messages.append(message))
    status.initializeFilenames(filenames)
    status.initializeComparisonStarCount(1)

    for filename in filenames:
        status.setCurrentFilename(filename)
        status.overexposedWarning(1, 100.0, 200.0, 50000.0)

    output_path = tmp_path / 'PlateStatus.csv'
    status.writePlateStatus(output_path)

    with output_path.open(newline='') as handle:
        rows = list(csv.reader(handle))
    header = rows[0]
    overexposed_column = header.index('overexposed_comp1')
    assert [row[overexposed_column] for row in rows[1:]] == ['True', 'True', 'True']
    assert sum(
        'Comparison star #1 is overexposed in file' in message
        for message in messages
    ) == 1
    assert messages[-1] == '>-- Comparison star #1: overexposed in 3 frame(s).'
