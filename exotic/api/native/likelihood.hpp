// Live likelihood, distinct from conditional-orbit replay: no replay prior
// correction or additional geometry rejection is applied to sampler points.
#pragma once

EXOTIC_EXPORT int exotic_loglikelihood(
    std::size_t draws, std::size_t dimension, const double* points,
    std::size_t count, const double* times, const double* exposures, std::size_t exposure_count,
    const double* base, const int* mapping, const double* data, const double* inverse_error,
    const double* airmass, const double* fixed_airmass, const double* weights,
    const unsigned char* baseline_mask, int law, int geometry, int baseline_mode,
    int airmass_mode, int observed_valid, int duration_prior,
    double expected_duration, double sigma_log_duration, std::size_t samples,
    double tolerance, double minimum_exposure, double padding,
    double* scores, char* error, std::size_t error_size
) noexcept {
    try {
        std::vector<double> model(count), residual(count), denominator, numerator, ratios;
        denominator.reserve(count); numerator.reserve(count); ratios.reserve(count);
        const double rad = pi / 180.;
        for (std::size_t d = 0; d < draws; ++d) {
            scores[d] = -1.e100; // Existing BAD_LOG_LIKELIHOOD, including invalid observations.
            if (!observed_valid) continue;
            double row[14];
            for (int j = 0; j < 14; ++j)
                row[j] = mapping[j] < 0 ? base[j] : points[d * dimension + mapping[j]];
            if (geometry) {
                double denom = 1. + row[3] * std::sin(row[5] * rad);
                if (std::abs(denom) <= 1.e-8) denom = std::numeric_limits<double>::epsilon();
                double scale = row[2] * (1. - row[3] * row[3]) / denom;
                if (std::abs(scale) <= 1.e-8) scale = std::numeric_limits<double>::epsilon();
                double cosine = row[9] / scale;
                if (cosine < -1.) cosine = -1.;
                if (cosine > 1.) cosine = 1.;
                row[4] = std::acos(cosine) * (180. / pi);
            }
            const Orbit orbit(row, law);
            double duration_loglike = 0.;
            if (duration_prior) {
                const double duration = orbit.duration();
                if (!std::isfinite(duration) || duration <= 0.) continue;
                const double value = std::log(duration / expected_duration);
                const double scaled = value / sigma_log_duration;
                duration_loglike = -.5 * (scaled * scaled);
            }
            // The reference scorer penalizes a failing trial model, while
            // library-loading/ABI errors are raised before sampling in Python.
            try {
                integrated_flux(count, times, exposures, exposure_count, orbit, row + 10,
                                samples, tolerance, minimum_exposure, padding, model.data());
                if (airmass_mode) for (std::size_t i = 0; i < count; ++i)
                    model[i] *= airmass_mode == 1 ? fixed_airmass[i] : std::exp(row[7] * airmass[i]);
                double baseline = row[8];
                if (!baseline_mode) {
                    denominator.clear(); numerator.clear(); ratios.clear();
                    for (std::size_t i = 0; i < count; ++i)
                        if (baseline_mask[i] && std::isfinite(model[i]) && model[i] != 0.) {
                            denominator.push_back(weights[i] * (model[i] * model[i]));
                            numerator.push_back((weights[i] * data[i]) * model[i]);
                            const double ratio = data[i] / model[i];
                            if (std::isfinite(ratio)) ratios.push_back(ratio);
                        }
                    baseline = 1.;
                    if (!denominator.empty()) {
                        const double denom = array_sum(denominator.data(), denominator.size());
                        if (!std::isfinite(denom) || denom <= 0.) {
                            if (!ratios.empty()) {
                                std::sort(ratios.begin(), ratios.end());
                                const std::size_t middle = ratios.size() / 2;
                                baseline = ratios.size() % 2 ? ratios[middle]
                                    : (ratios[middle - 1] + ratios[middle]) / 2.;
                            }
                        } else baseline = array_sum(numerator.data(), numerator.size()) / denom;
                    }
                    if (!std::isfinite(baseline)) baseline = 1.;
                }
                bool finite = true;
                for (std::size_t i = 0; i < count; ++i) {
                    model[i] *= baseline;
                    finite = finite && std::isfinite(model[i]);
                    const double value = (data[i] - model[i]) * inverse_error[i];
                    residual[i] = value * value;
                }
                if (!finite) continue;
                const double score = -.5 * array_sum(residual.data(), count) + duration_loglike;
                if (std::isfinite(score)) scores[d] = score;
            } catch (const std::exception&) {
                // Same per-point invalid-model penalty as single_loglike.
            }
        }
        return 0;
    } catch (const std::exception& exc) {
        if (error_size) std::snprintf(error, error_size, "%s", exc.what());
        return 1;
    } catch (...) {
        if (error_size) std::snprintf(error, error_size, "%s", "Unknown native likelihood error");
        return 2;
    }
}
