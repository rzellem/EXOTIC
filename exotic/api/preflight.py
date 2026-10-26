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
"""Pre-flight checks for an initialization file, run before any photometry.

Every check here is cheap (FITS headers, one image, one archive query, at most
one plate solve) and targets a mistake that -ov otherwise lets through until an
hour of reduction has been spent on it:

  * planetary parameters copied from a template and never edited
    (archive agreement: period, Rp/Rs, a/Rs, inclination)
  * an ephemeris whose transit does not fall in the observing window
  * a target pixel that is not on the target (checked against a WCS)
  * comparison stars off the frame, saturated, or within noise of blank sky

The functions in this module take plain values (dicts, arrays, a WCS) so they
can be tested without frames or network; ``exotic.py`` wires them to the
initialization file, the NASA Exoplanet Archive and the plate solvers.
"""
import numpy as np

PASS = 'pass'
FAIL = 'fail'
LOOK = 'look'
SKIP = 'skip'

_STATUS_MARK = {PASS: 'ok  ', FAIL: 'FAIL', LOOK: 'look', SKIP: '    '}

# Relative tolerances for archive agreement. The period is held tight because
# it is the parameter a template copy gets wrong by a lot (a wrong period puts
# the predicted transit anywhere), while the geometry may legitimately differ
# between the archive default row and a newer literature source.
ARCHIVE_TOLERANCES = {
    'pPer': ('Orbital Period (days)', 1e-4),
    'rprs': ('Ratio of Planet to Stellar Radius (Rp/Rs)', 0.15),
    'aRs': ('Ratio of Distance to Stellar Radius (a/Rs)', 0.15),
    'inc': ('Orbital Inclination (deg)', 0.02),
}

TARGET_PIXEL_TOLERANCE_PX = 8.0
COMPARISON_MIN_SNR = 5.0
COMPARISON_BRIGHTNESS_RATIO_RANGE = (0.1, 20.0)
PEAK_SEARCH_HALF_WIDTH_PX = 5


class PreflightReport:
    """Ordered list of (status, label, detail) with a hard-failure verdict."""

    def __init__(self):
        self.checks = []

    def add(self, status, label, detail=''):
        self.checks.append((status, label, detail))
        return status

    def passed(self, ok, label, detail=''):
        return self.add(PASS if ok else FAIL, label, detail)

    def look(self, ok, label, detail=''):
        """Report without failing: a judgment call, not a rule."""
        return self.add(PASS if ok else LOOK, label, detail)

    def skip(self, label, detail=''):
        return self.add(SKIP, label, detail)

    @property
    def failures(self):
        return [label for status, label, _ in self.checks if status == FAIL]

    @property
    def ok(self):
        return not self.failures

    def lines(self):
        for status, label, detail in self.checks:
            yield f"  [{_STATUS_MARK[status]}] {label}" + (f"  {detail}" if detail else '')

    def summary(self):
        if self.failures:
            return f"{len(self.failures)} pre-flight check(s) failed. Do not reduce with this file yet."
        return "All pre-flight checks passed."


def _finite(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def compare_archive_parameters(planet_dict, archive_dict, report, tolerances=ARCHIVE_TOLERANCES):
    """Check the initialization file's planetary parameters against the archive."""
    for key, (label, tolerance) in tolerances.items():
        mine = _finite(planet_dict.get(key))
        theirs = _finite(archive_dict.get(key)) if archive_dict else None
        if theirs is None or theirs == 0:
            report.skip(label, 'no archive value to compare against')
            continue
        if mine is None:
            report.passed(False, label, f"missing in the initialization file (archive {theirs})")
            continue
        relative = abs(mine - theirs) / abs(theirs)
        report.passed(relative <= tolerance, label,
                      f"inits {mine}  archive {theirs}  ({relative * 100:.2f}% off, tolerance {tolerance * 100:g}%)")
    return report


def predicted_mid_transit(jd_start, jd_end, mid_transit, period):
    """Epoch of the published ephemeris nearest the middle of the observing window."""
    cycle = np.round(((jd_start + jd_end) / 2.0 - mid_transit) / period)
    return float(mid_transit + cycle * period), int(cycle)


def check_transit_window(jd_start, jd_end, mid_transit, period, duration_days, report):
    """Fail if no part of the predicted transit falls in the observing window.

    A partial transit is reported for the observer's judgment but is not a
    failure: the pipeline fits partials. The failure this exists for is the
    template-copy ephemeris whose transit is nowhere near the night.
    """
    jd_start, jd_end = _finite(jd_start), _finite(jd_end)
    mid_transit, period = _finite(mid_transit), _finite(period)
    if jd_start is None or jd_end is None or jd_end < jd_start:
        report.passed(False, 'observing window from the FITS headers', 'could not read start and end times')
        return report
    if mid_transit is None or period is None or period <= 0:
        report.passed(False, 'ephemeris usable', f"Tmid {mid_transit} period {period}")
        return report

    tmid, cycle = predicted_mid_transit(jd_start, jd_end, mid_transit, period)
    hours = (jd_end - jd_start) * 24.0
    offset_hours = (tmid - jd_start) * 24.0
    report.add(SKIP, 'observing window', f"{jd_start:.5f} -> {jd_end:.5f}  ({hours:.2f} h, epoch {cycle})")

    duration = _finite(duration_days)
    if duration is None or duration <= 0:
        report.passed(jd_start <= tmid <= jd_end, 'predicted mid-transit inside the window',
                      f"Tmid {tmid:.5f} ({offset_hours:+.2f} h from start); no duration available")
        return report

    ingress, egress = tmid - duration / 2.0, tmid + duration / 2.0
    overlaps = egress >= jd_start and ingress <= jd_end
    report.passed(overlaps, 'predicted transit overlaps the window',
                  f"Tmid {tmid:.5f} ({offset_hours:+.2f} h from start), duration {duration * 24:.2f} h")
    if overlaps:
        report.look(ingress >= jd_start and egress <= jd_end, 'full transit inside the window',
                    f"ingress {(ingress - jd_start) * 24:+.2f} h  egress {(egress - jd_start) * 24:+.2f} h from start"
                    + ('' if ingress >= jd_start and egress <= jd_end else '  (partial transit)'))
    return report


def check_target_pixel(wcs, ra_deg, dec_deg, target_xy, report, tolerance_px=TARGET_PIXEL_TOLERANCE_PX,
                       frame_label='first frame'):
    """Project the target's RA/Dec through a WCS and compare with the seed pixel."""
    ra_deg, dec_deg = _finite(ra_deg), _finite(dec_deg)
    if ra_deg is None or dec_deg is None:
        report.skip('target pixel on the target', 'no target RA/Dec available')
        return report
    try:
        expected = np.asarray(wcs.all_world2pix(ra_deg, dec_deg, 0), dtype=float)
        seed = np.asarray(target_xy, dtype=float)
        separation = float(np.hypot(*(expected - seed)))
    except Exception as exc:
        report.skip('target pixel on the target', f"WCS projection failed: {exc}")
        return report
    report.passed(separation <= tolerance_px, f"target pixel on the target ({frame_label})",
                  f"WCS puts the target at ({expected[0]:.1f}, {expected[1]:.1f}); "
                  f"inits pixel ({seed[0]:.0f}, {seed[1]:.0f}) is {separation:.1f} px away (tolerance {tolerance_px:g})")
    return report


def peak_above_background(image, x, y, background, half_width=PEAK_SEARCH_HALF_WIDTH_PX):
    """Brightest pixel in a small box around (x, y), minus the background; None if off the frame."""
    height, width = image.shape[:2]
    if not (0 <= x < width and 0 <= y < height):
        return None
    x0, x1 = max(0, int(x) - half_width), int(x) + half_width + 1
    y0, y1 = max(0, int(y) - half_width), int(y) + half_width + 1
    return float(np.nanmax(image[y0:y1, x0:x1])) - background


def check_comparison_stars(image, target_xy, comparison_xy, saturation_value, report,
                           overexposure_fraction=0.9):
    """Comparison stars must be on the frame and below the overexposure threshold.

    Brightness relative to the target is a judgment call reported without failing:
    one frame is a small sample, and on a clouded night the first frame can be the
    worst one. Off-frame and overexposed are rules; the pipeline would reject them
    later anyway, and this says so before the photometry is run.
    """
    image = np.asarray(image, dtype=float)
    background = float(np.nanmedian(image))
    quiet = image[image < np.nanpercentile(image, 95)]
    noise = float(np.nanstd(quiet)) if quiet.size else float('nan')
    target_peak = peak_above_background(image, target_xy[0], target_xy[1], background)
    threshold = saturation_value * overexposure_fraction if _finite(saturation_value) else None

    if target_peak is None:
        report.passed(False, 'target pixel inside the frame', f"({target_xy[0]}, {target_xy[1]}) is off the frame")
    else:
        report.add(SKIP, 'first-frame photometry context',
                   f"background {background:.0f}, noise {noise:.0f}, target peak {target_peak:+.0f} ADU above background")
        if threshold is not None:
            report.passed(target_peak + background < threshold, 'target unsaturated',
                          f"peak {target_peak + background:.0f} vs threshold {threshold:.0f}")

    if not comparison_xy:
        report.skip('comparison stars', 'none given in the initialization file')
        return report

    for comp in comparison_xy:
        try:
            cx, cy = float(comp[0]), float(comp[1])
        except (TypeError, ValueError, IndexError):
            report.passed(False, f"comp {comp!r}", 'not an [x, y] pair')
            continue
        label = f"comp ({cx:.0f}, {cy:.0f})"
        peak = peak_above_background(image, cx, cy, background)
        if peak is None:
            report.passed(False, f"{label} inside the frame", 'outside the frame')
            continue
        if threshold is not None:
            report.passed(peak + background < threshold, f"{label} unsaturated",
                          f"peak {peak + background:.0f} vs threshold {threshold:.0f}")
        if target_peak and target_peak > 0:
            ratio = peak / target_peak
            low, high = COMPARISON_BRIGHTNESS_RATIO_RANGE
            report.look(peak > COMPARISON_MIN_SNR * noise and low <= ratio <= high, f"{label} brightness",
                        f"{peak:+.0f} ADU above background, {ratio:.1f}x the target")
        else:
            report.look(peak > COMPARISON_MIN_SNR * noise, f"{label} brightness",
                        f"{peak:+.0f} ADU above background")
    return report
