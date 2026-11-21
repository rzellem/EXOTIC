// Exact detector low-tail histogram and simultaneous local defect repair.
// Compact dense counts use at most one entry per input pixel; wide ADU spans
// use ordered sparse counts. Neither representation sorts/hashes the image.
#pragma once
#include <map>
#include <cstdint>

namespace exotic_image {
double low_threshold(const double* image, std::size_t n, std::int64_t gap, double population) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double low = std::numeric_limits<double>::infinity(), high = -low;
    for (std::size_t i = 0; i < n; ++i) if (std::isfinite(image[i])) {
        if (image[i] < low) low = image[i];
        if (image[i] > high) high = image[i];
    }
    if (!std::isfinite(low)) return nan;
    // Truncation is monotonic. Only the two extrema need scalar libm calls;
    // ordinary ADU bins use the exact truncating integer instruction below.
    low = std::trunc(low); high = std::trunc(high);
    const double span = high - low;
    std::vector<std::size_t> dense;
    std::map<double, std::size_t> sparse;
    const bool compact = std::isfinite(span) && span < static_cast<double>(n);
    if (compact) {
        dense.resize(static_cast<std::size_t>(span) + 1);
        if (low >= -9223372036854775808. && high < 9223372036854775808.) {
            const std::int64_t first = static_cast<std::int64_t>(low);
            for (std::size_t i = 0; i < n; ++i) if (std::isfinite(image[i]))
                ++dense[static_cast<std::size_t>(static_cast<std::int64_t>(image[i]) - first)];
        } else {
            for (std::size_t i = 0; i < n; ++i) if (std::isfinite(image[i]))
                ++dense[static_cast<std::size_t>(std::trunc(image[i]) - low)];
        }
    } else {
        for (std::size_t i = 0; i < n; ++i) if (std::isfinite(image[i]))
            ++sparse[std::trunc(image[i])];
    }
    double mode = low;
    std::size_t largest = 0;
    if (compact) {
        for (std::size_t j = 0; j < dense.size(); ++j) if (dense[j] > largest) {
            largest = dense[j]; mode = low + static_cast<double>(j);
        }
    } else {
        for (const auto& bin : sparse) if (bin.second > largest) {
            largest = bin.second; mode = bin.first;
        }
    }
    const auto count = [&](double bin) -> std::size_t {
        if (compact) {
            const double index = bin - low;
            if (index < 0. || index > span) return 0;
            return dense[static_cast<std::size_t>(index)];
        }
        const auto it = sparse.find(bin);
        return it == sparse.end() ? 0 : it->second;
    };
    // Keep the existing search extent, strict population comparison, and
    // double-precision threshold arithmetic (including very large ADU bins).
    for (int offset = 1; offset < 10000; ++offset) {
        const double threshold = mode - offset;
        bool quiet = true;
        for (std::int64_t step = 0; step < gap; ++step)
            if (!(static_cast<double>(count(threshold - static_cast<double>(step))) < population)) {
                quiet = false; break;
            }
        if (quiet) return threshold;
    }
    return nan;
}
}

EXOTIC_EXPORT int exotic_low_pixel_mask(
    std::size_t n, const double* image, std::int64_t gap, double population,
    unsigned char* mask, double* threshold, char* error, std::size_t capacity
) noexcept {
    try {
        *threshold = exotic_image::low_threshold(image, n, gap, population);
        for (std::size_t i = 0; i < n; ++i)
            mask[i] = std::isfinite(image[i]) && image[i] < *threshold;
        return 0;
    } catch (const std::exception& exc) {
        if (capacity) std::snprintf(error, capacity, "%s", exc.what());
        return 1;
    }
}

EXOTIC_EXPORT int exotic_repair_pixels(
    std::size_t ny, std::size_t nx, double* image, const unsigned char* reference,
    int detect_low, unsigned char* mask, std::size_t* unresolved,
    char* error, std::size_t capacity
) noexcept {
    try {
        const std::size_t n = ny * nx;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double threshold = detect_low ? exotic_image::low_threshold(image, n, 17, 25.) : nan;
        std::vector<std::size_t> bad;
        *unresolved = 0;
        for (std::size_t i = 0; i < n; ++i) {
            mask[i] = reference[i] || !std::isfinite(image[i]) || image[i] < threshold;
            if (mask[i]) bad.push_back(i);
        }
        if (bad.empty()) return 0;
        if (bad.size() == n) throw std::invalid_argument("Cannot infill a frame with no valid pixels");
        for (std::size_t i : bad) {
            const std::size_t y = i / nx, x = i % nx;
            double neighbours[8];
            std::size_t count = 0;
            for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
                if ((!dy && !dx) || (dy < 0 && !y) || (dx < 0 && !x)
                    || (dy > 0 && y + 1 == ny) || (dx > 0 && x + 1 == nx)) continue;
                const std::size_t j = (y + dy) * nx + x + dx;
                // The original combined mask excludes ALL defects, including
                // earlier repairs; only original good values are neighbours.
                if (!mask[j]) neighbours[count++] = image[j];
            }
            if (count) {
                std::sort(neighbours, neighbours + count);
                image[i] = count % 2 ? neighbours[count / 2]
                    : (neighbours[count / 2 - 1] + neighbours[count / 2]) / 2.;
            } else {
                image[i] = nan;
                mask[i] = 2;
                ++*unresolved;
            }
        }
        return 0;
    } catch (const std::exception& exc) {
        if (capacity) std::snprintf(error, capacity, "%s", exc.what());
        return 1;
    }
}
