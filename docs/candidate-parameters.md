# Optional NextAstro candidate parameter auto-fill

Interactive candidate runs offer **Auto-fill missing candidate parameters from
NextAstro? (y/n)** when values are missing. Accepting fills the available values
without asking for each estimate again. Other interactive input confirmations
continue normally.

For an initialization file, add this to its existing `optional_info` object:

```json
"use_nextastro_candidate_estimator": true
```

Unattended candidate runs enable the estimator automatically when this setting
is absent, including unattended override (`-ov`) runs. They never stop to ask.
An explicit `false` disables the estimator without prompting. Interactive runs
offer the choice when the setting is absent; `true` selects it directly,
including interactive override runs. Normal confirmed-planet catalogue runs
retain their existing behaviour.

EXOTIC sends the candidate name, coordinates and current planetary parameters to
`https://archive.nextastro.org/api/candidate_parameters`. It preserves supplied
values and fills missing ones from Archive's confirmed-planet/TOI cache, Gaia DR3
stellar values, target-specific ratio calculations, and Archive's existing
EXOTIC setup starting approximations.
Gaia is matched with Archive's existing nearest-source lookup; the source ID and
separation are reported for review. Images and light curves are not uploaded.

Every filled value is labelled `catalogue`, `gaia`, `derived`, or `approximation` in the
console and the local `working_artifacts/NextAstroCandidateParameters.json`
report. Starting approximations and their widths are not measurements. For
example, without better information the starting radius ratio is 0.15 ± 0.125,
a/Rs is 20 ± 15, inclination is 85 ± 20 degrees, and effective temperature is
5500 ± 200 K. These are the existing Archive setup values, not a stellar fit.
Period and transit epoch have no invented defaults; unresolved values still
need input. Missing catalogue timing standards remain unknown.

Before using a fixed a/Rs guess, Archive calculates
`a/Rs = (G M_star P^2 / (4 pi^2))^(1/3) / R_star` from the supplied/catalogue
period and Gaia FLAME mass/radius. GSP-Phot radius is used if FLAME radius is
absent. If FLAME mass is absent, Gaia radius and `logg_gspphot` supply an inferred
mass through `M_star = g R_star^2/G`, where `g = 10^logg` in cm/s². The inferred
mass and its input columns/method are saved in the derivation report. The generic
setup gravity is never used to manufacture a stellar mass.
This neglects planet mass. Available input uncertainties are propagated
to first order assuming independence; missing uncertainty terms and omitted
covariance/systematics are recorded in the report. A missing error does not
prevent calculating the ratio.

For the gravity-derived mass, propagation uses the original period, gravity and
radius directly. Since `a/R = (g P^2/(4 pi^2 R))^(1/3)`, the fractional radius
error contributes with a coefficient of 1/3. This accounts for the shared radius
in the inferred mass. Remaining Gaia parameter covariance/systematics and
missing errors remain explicit limitations.

Before using the fixed Rp/Rs guess, Archive uses `sqrt(fractional transit depth)`
where a supplied or catalogue depth is available. Existing user and published
catalogue ratios remain preferred. Depth uncertainty, when present, is
propagated to first order. The depth estimate does not correct for dilution,
limb darkening or grazing geometry. A ratio can be calculated even when its
uncertainty must still use a labelled starting width.

To supply your own depth, put the following in `optional_info` alongside the
estimator opt-in (this example is a 1% transit with a 100 ppm depth uncertainty):

```json
"transit_depth": 10000,
"transit_depth_unc": 100,
"transit_depth_unit": "ppm"
```

Accepted units are `fraction` (default), `percent`/`%`, `ppt` and `ppm`.
These are estimator inputs; they do not add required planetary parameters.

Before using the fixed inclination guess, Archive uses total first-to-fourth
contact duration T14, period, a/Rs and Rp/Rs to estimate impact parameter and
inclination. When both eccentricity and argument of periastron are known, the
Winn duration approximation applies both the separation factor
`(1-e^2)/(1+e sin(omega))` and speed factor `sqrt(1-e^2)/(1+e sin(omega))`.
The report records `eccentric_T14` and the actual e/omega used. If either is
unknown, it records a `circular_T14` estimate and the missing orbital information.
Known e=0 remains circular without needing omega; omega=0 is a valid known angle.
Supplied/published inclinations and supplied e/omega remain preferred/preserved.
The existing EXOTIC omega convention is used (planet argument, degrees); an
ambiguous stellar-versus-planet argument is not automatically shifted by 180°.
This is a 0–90 degree estimate and does not infer the prograde/retrograde branch.

Durations are obtained from NASA/TOI records or supplied in `optional_info`:

```json
"transit_duration": 3.0,
"transit_duration_unc": 0.1,
"transit_duration_unit": "hours"
```

Available e/omega uncertainties can also be supplied in `optional_info` as
`"ecc_unc": 0.02` and `"omega_unc": 3.0` (degrees). Catalogue uncertainties
are used when available and appropriate to the nominal e/omega. Missing errors
never prevent a nominal estimate; unquantified terms are reported.

Accepted duration units are `hours` (default), `days`, `minutes`, and `seconds`.
Supply total T14, not the full-transit T23 or a light-curve width at half depth.
The calculation uses available geometry before fixed a/Rs/Rp/Rs guesses are
inserted. Missing input errors do not prevent the estimate. Approximate errors
are propagated in cos²(i) and transformed to angles, remaining finite at central
transits. Shared period dependence in an a/Rs calculated by this request is
accounted for. Available eccentricity/omega errors are propagated through both
correction factors. Missing errors, other covariance, approximation error and
systematics are recorded as limitations. The eccentric relation remains an
approximation; see [Winn's duration equations](https://arxiv.org/abs/1001.2010).

If T14 is too long for the chosen model, the endpoint returns the nearest
edge-on (90°) estimate and records the duration mismatch. The mismatch is not
part of the quoted uncertainty. Inputs that do not identify detached-transit
geometry (a/Rs<=1+Rp/Rs, or the corresponding eccentric transit separation), or
missing duration/period/ratios, leave the labelled fixed inclination available
and report any mathematical failure. An invalid known e outside 0<=e<1 is
reported; it is not silently replaced with a circular model.

Gaia GSP-Phot estimates temperature, log(g) and [M/H]. Its 16th/84th percentile
bounds supply asymmetric uncertainty widths. [M/H] is used as an approximation
to EXOTIC's [Fe/H]; those intervals omit systematic uncertainty. See the
[ESA GSP-Phot documentation](https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu8par/sec_cu8par_apsis/ssec_cu8par_apsis_gspphot.html).

Lookup failures are printed. A failed endpoint request leaves inputs unchanged;
partial responses label any starting approximations and report unavailable
catalogue/Gaia information. Users can fill remaining parameters through the
existing input workflow. The integration uses HTTP only and does not require an
Archive checkout or account.
