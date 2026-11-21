// Batch C interface to the vendored MIT-licensed limbdarkfast C++ kernel.
// Orbital convention and midpoint exposure integration match EXOTIC's existing
// Python/pylightcurve path. Double precision, without fast-math.
#include "limbdarkfast.hpp"
#include "claret.hpp"
#include "photometry.hpp"
#include <cstddef>
#include <cstdio>
#include <exception>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#ifdef _WIN32
#define EXOTIC_EXPORT extern "C" __declspec(dllexport)
#else
#define EXOTIC_EXPORT extern "C" __attribute__((visibility("default")))
#endif

EXOTIC_EXPORT int exotic_hybrid4_abi_version() { return 6; }

#include "image_preparation.hpp"

EXOTIC_EXPORT int exotic_psf_fit(size_t n,const double* x,const double* y,const double* data,
    const double* initial,const double* lo,const double* hi,int method,exotic_photo::Svd svd,
    double* output,double* info,char* error,size_t capacity) {
    try { exotic_photo::fit(x,y,data,n,initial,lo,hi,method,svd,output,info); return 0; }
    catch(const std::exception& e){std::snprintf(error,capacity,"%s",e.what());return 1;}
}
EXOTIC_EXPORT int exotic_aperture_grid(const double* image,int nx,int ny,int ox,int oy,
    double xc,double yc,size_t na,const double* radii,size_t nb,const double* inner,
    const double* outer,int center,double ptol,double sigma,int iters,double* ap,double* bg,
    char* error,size_t capacity) {
    try {
        for(size_t i=0;i<na;++i) {
            if(std::isfinite(radii[i])&&radii[i]>0)
                exotic_photo::aperture(image,nx,ny,ox,oy,xc,yc,radii[i],center,ap+3*i);
            else for(int k=0;k<3;++k)ap[3*i+k]=exotic_photo::nan;
            for(size_t j=0;j<nb;++j){size_t k=i*nb+j;
                if(outer[k]==0){bg[3*k]=bg[3*k+1]=bg[3*k+2]=0;}
                else if(std::isfinite(inner[k])&&std::isfinite(outer[k])&&outer[k]>inner[k])
                    exotic_photo::background(image,nx,ny,ox,oy,xc,yc,inner[k],outer[k],center,ptol,sigma,iters,bg+3*k);
                else for(int c=0;c<3;++c)bg[3*k+c]=exotic_photo::nan;
            }
        }return 0;
    }catch(const std::exception& e){std::snprintf(error,capacity,"%s",e.what());return 1;}
}

namespace {
const double pi = 3.141592653589793238462643383279502884;

struct Orbit {
    int law;
    double radius, period, ars, ecc, inclination, omega, mid;
    double sini, cosi, reference_mid, eccentric_factor;

    explicit Orbit(const double* p,int model_law=0)
        : law(model_law), radius(p[0]), period(p[1]), ars(p[2]), ecc(p[3]),
          inclination(p[4] * pi / 180.0), omega(p[5] * pi / 180.0), mid(p[6]) {
        sini = std::sin(inclination);
        cosi = std::cos(inclination);
        reference_mid = mid;
        eccentric_factor = 1.0;
        if (ecc != 0.0) {
            const double at_mid = (omega < pi / 2.0 ? pi / 2.0 : 5.0 * pi / 2.0) - omega;
            double eccentric_mid = 2.0 * std::atan(std::sqrt((1.0-ecc)/(1.0+ecc)) * std::tan(at_mid/2.0));
            if (eccentric_mid < 0.0) eccentric_mid += 2.0 * pi;
            reference_mid = mid - (period / 2.0 / pi) * (eccentric_mid-ecc*std::sin(eccentric_mid));
            eccentric_factor = std::sqrt((1.0+ecc)/(1.0-ecc));
        }
    }

    double duration() const {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        if (!std::isfinite(period) || period <= 0.0 || !std::isfinite(radius) || radius < 0.0
            || !std::isfinite(ars) || ars <= 0.0 || !std::isfinite(inclination)
            || !std::isfinite(sini) || sini <= 0.0) return nan;
        const double denominator = std::max(std::numeric_limits<double>::epsilon(), 1.0+ecc*std::sin(omega));
        const double scale = ars * (1.0-ecc*ecc) / denominator;
        const double impact = scale * cosi;
        const double chord = (1.0+radius)*(1.0+radius) - impact*impact;
        if (!std::isfinite(chord) || chord <= 0.0 || !std::isfinite(scale) || scale <= 0.0) return nan;
        const double argument = std::max(-1.0, std::min(1.0, std::sqrt(chord)/(ars*sini)));
        return (period/pi) * std::asin(argument) * (std::sqrt(1.0-ecc*ecc)/denominator);
    }

    void flux(const double* times, std::size_t count, const double* weights, double* output) const {
        std::vector<double> anomaly;
        if (ecc != 0.0) {
            std::vector<double> mean(count);
            anomaly.resize(count);
            for (std::size_t i=0; i<count; ++i) {
                const double delta = times[i]-reference_mid;
                // Truncate toward zero, including before the reference epoch.
                mean[i] = (delta-std::trunc(delta/period)*period)*2.0*pi/period;
                anomaly[i] = mean[i];
            }
            bool converged = false;
            // Match the existing batch Newton solver's all-points stopping
            // rule, tolerance and iteration limit, rather than stopping each
            // point early and changing its last iterate.
            for (int iteration=0; iteration<10000; ++iteration) {
                converged = true;
                for (std::size_t i=0; i<count; ++i) {
                    const double old = anomaly[i];
                    anomaly[i] = old - (old-ecc*std::sin(old)-mean[i])/(1.0-ecc*std::cos(old));
                    if (!(std::abs(anomaly[i]-old) < 1.e-7)) converged = false;
                }
                if (converged) break;
            }
            if (!converged) throw std::runtime_error("Failed to find a solution in 10000 loops");
        }
        for (std::size_t i=0; i<count; ++i) {
            double x, y, z;
            if (ecc == 0.0) {
                const double phase = 2.0*pi*(times[i]-mid)/period;
                const double projected = ars*std::cos(phase);
                x = projected*sini;
                y = ars*std::sin(phase);
                z = -projected*cosi;
            } else {
                const double true_anomaly = 2.0*std::atan(eccentric_factor*std::tan(anomaly[i]/2.0));
                const double distance = ars*(1.0-ecc*ecc)/(1.0+ecc*std::cos(true_anomaly));
                const double a = std::cos(true_anomaly+omega);
                const double b = std::sin(true_anomaly+omega);
                x = distance*b*sini;
                y = distance*(-a);
                z = distance*(-b*cosi);
            }
            if (law == 1) {
                const double separation=x<0.?1.+5.*radius:std::sqrt(y*y+z*z);
                output[i]=exotic_claret::flux(separation,radius,weights);
            } else {
                const double value = limbdarkfast::hybrid_flux(std::hypot(y,z),radius,weights,limbdarkfast::Law::hybrid4);
                output[i] = x < 0.0 ? 1.0 : value;
            }
        }
    }
};

void integrated_flux(std::size_t count, const double* times, const double* exposures,
                     std::size_t exposure_count, const Orbit& orbit, const double* weights,
                     std::size_t samples, double tolerance, double minimum_exposure,
                     double padding, double* output) {
    orbit.flux(times,count,weights,output);
    if (!count || !exposure_count) return;
    auto exposure = [&](std::size_t i) { return exposures[exposure_count == 1 ? 0 : i]; };
    std::vector<std::size_t> candidates;
    double maximum_exposure = 0.0;
    for (std::size_t i=0; i<count; ++i) {
        const double value = exposure(i);
        if (std::isfinite(times[i]) && std::isfinite(value) && value > minimum_exposure) {
            candidates.push_back(i);
            maximum_exposure = std::max(maximum_exposure,value);
        }
    }
    const double duration = orbit.duration();
    if (std::isfinite(orbit.period) && orbit.period > 0.0 && std::isfinite(orbit.mid)
        && std::isfinite(duration) && duration > 0.0) {
        const double half_window = 0.5*duration + (0.5+padding)*maximum_exposure;
        candidates.erase(std::remove_if(candidates.begin(),candidates.end(),[&](std::size_t i) {
            double phase = std::fmod((times[i]-orbit.mid+0.25*orbit.period)/orbit.period,1.0);
            if (phase < 0.0) phase += 1.0;
            return !(std::abs((phase-0.25)*orbit.period) <= half_window);
        }),candidates.end());
    }
    if (candidates.empty()) return;
    std::vector<double> starts(candidates.size()), ends(candidates.size());
    std::vector<double> first(candidates.size()), last(candidates.size());
    for (std::size_t j=0; j<candidates.size(); ++j) {
        const std::size_t i = candidates[j];
        const double half = 0.5*exposure(i);
        starts[j] = times[i]-half;
        ends[j] = times[i]+half;
    }
    orbit.flux(starts.data(),starts.size(),weights,first.data());
    orbit.flux(ends.data(),ends.size(),weights,last.data());
    std::vector<std::size_t> active;
    for (std::size_t j=0; j<candidates.size(); ++j) {
        const double a = std::abs(first[j]-output[candidates[j]]);
        const double b = std::abs(last[j]-output[candidates[j]]);
        const double c = std::abs(last[j]-first[j]);
        if (!std::isnan(a) && !std::isnan(b) && !std::isnan(c) && std::max(a,std::max(b,c)) > tolerance)
            active.push_back(candidates[j]);
    }
    if (active.empty()) return;
    std::vector<double> sample_times, flux;
    if (samples > sample_times.max_size()/active.size())
        throw std::length_error("Exposure sample array exceeds addressable memory");
    sample_times.resize(active.size()*samples);
    flux.resize(sample_times.size());
    for (std::size_t j=0; j<active.size(); ++j)
        for (std::size_t k=0; k<samples; ++k) {
            const double offset = (static_cast<double>(k)+0.5)/static_cast<double>(samples)-0.5;
            sample_times[j*samples+k] = times[active[j]]+exposure(active[j])*offset;
        }
    orbit.flux(sample_times.data(),sample_times.size(),weights,flux.data());
    for (std::size_t j=0; j<active.size(); ++j) {
        double sum = 0.0;
        for (std::size_t k=0; k<samples; ++k) sum += flux[j*samples+k];
        output[active[j]] = sum/static_cast<double>(samples);
    }
}
} // namespace

EXOTIC_EXPORT int exotic_hybrid4_transit(
    std::size_t count, const double* times, const double* exposures, std::size_t exposure_count,
    const double* parameters, const double* weights, std::size_t samples,
    double tolerance, double minimum_exposure, double padding,
    double* output, char* error, std::size_t error_size
) noexcept {
    try {
        if (exposure_count != 0 && exposure_count != 1 && exposure_count != count)
            throw std::invalid_argument("Exposure array must be scalar or match times");
        if (!samples) throw std::invalid_argument("Exposure integration needs samples");
        integrated_flux(count,times,exposures,exposure_count,Orbit(parameters),weights,
                        samples,tolerance,minimum_exposure,padding,output);
        return 0;
    } catch (const std::exception& exc) {
        if (error_size) std::snprintf(error,error_size,"%s",exc.what());
        return 1;
    } catch (...) {
        if (error_size) std::snprintf(error,error_size,"%s","Unknown C++ orbit/exposure error");
        return 2;
    }
}

EXOTIC_EXPORT int exotic_hybrid4_flux(
    std::size_t count, const double* separation, const double* radius,
    const double* weights, double* output, char* error, std::size_t error_size
) noexcept {
    try {
        for (std::size_t i = 0; i < count; ++i) {
            output[i] = limbdarkfast::hybrid_flux(
                separation[i], radius[i], weights, limbdarkfast::Law::hybrid4);
        }
        return 0;
    } catch (const std::exception& exc) {
        if (error_size) std::snprintf(error, error_size, "%s", exc.what());
        return 1;
    } catch (...) {
        if (error_size) std::snprintf(error, error_size, "%s", "Unknown C++ hybrid4 error");
        return 2;
    }
}


EXOTIC_EXPORT int exotic_claret_transit(
    std::size_t count,const double* times,const double* exposures,std::size_t exposure_count,
    const double* parameters,const double* coefficients,std::size_t samples,
    double tolerance,double minimum_exposure,double padding,
    double* output,char* error,std::size_t error_size
) noexcept {
    try {
        if (exposure_count!=0 && exposure_count!=1 && exposure_count!=count)
            throw std::invalid_argument("Exposure array must be scalar or match times");
        if (!samples) throw std::invalid_argument("Exposure integration needs samples");
        integrated_flux(count,times,exposures,exposure_count,Orbit(parameters,1),coefficients,
                        samples,tolerance,minimum_exposure,padding,output);
        return 0;
    } catch (const std::exception& exc) {
        if (error_size) std::snprintf(error,error_size,"%s",exc.what());
        return 1;
    } catch (...) {
        if (error_size) std::snprintf(error,error_size,"%s","Unknown C++ Claret error");
        return 2;
    }
}

namespace {
// Pairwise float64 reduction, with eight interleaved accumulators. This keeps
// NumPy likelihood reductions and native reductions aligned without fast-math.
double array_sum(const double* a,std::size_t n) {
    if (n<8) {
        double result=-0.;
        for (std::size_t i=0;i<n;++i) result+=a[i];
        return result;
    }
    if (n<=128) {
        double r[8];
        for (int j=0;j<8;++j) r[j]=a[j];
        std::size_t i=8;
        for (;i+7<n;i+=8) for (int j=0;j<8;++j) r[j]+=a[i+j];
        double result=((r[0]+r[1])+(r[2]+r[3]))+((r[4]+r[5])+(r[6]+r[7]));
        for (;i<n;++i) result+=a[i];
        return result;
    }
    std::size_t half=n/2;half-=half%8;
    return array_sum(a,half)+array_sum(a+half,n-half);
}
}

#include "likelihood.hpp"

// Draw columns: rprs, per, ars, ecc, inc, omega, tmid, a2, baseline, b, u0..u3.
// Baseline mode 0 profiles it using the original mask/weights; mode 1 uses the
// supplied per-draw value (including fixed and explicitly sampled baselines).
EXOTIC_EXPORT int exotic_orbital_replay(
    std::size_t draws,std::size_t count,const double* parameters,
    const double* times,const double* exposures,std::size_t exposure_count,
    const double* data,const double* inverse_error,const double* airmass,
    const double* baseline_weights,const unsigned char* baseline_mask,
    int law,int geometry,int baseline_mode,double ecc,double omega,
    int duration_prior,double expected_duration,double sigma_log_duration,
    std::size_t samples,double tolerance,double minimum_exposure,double padding,
    double* scores,double* inclinations,double* impacts,double* expansions,
    char* error,std::size_t error_size
) noexcept {
    try {
        std::vector<double> model(count), residual(count), denominator, numerator, ratios;
        denominator.reserve(count);numerator.reserve(count);ratios.reserve(count);
        const double rad=pi/180.;
        const double new_factor=(1.-ecc*ecc)/(1.+ecc*std::sin(omega*rad));
        for (std::size_t d=0;d<draws;++d) {
            const double* row=parameters+14*d;
            double p[7];std::copy(row,row+7,p);
            p[3]=ecc;p[5]=omega;
            scores[d]=-std::numeric_limits<double>::infinity();
            inclinations[d]=row[4];impacts[d]=row[9];expansions[d]=0.;
            double prior=0.;
            if (geometry==1) {
                const double old_factor=(1.-row[3]*row[3])/(1.+row[3]*std::sin(row[5]*rad));
                const double old_upper=std::min(1.+row[0],row[2]*old_factor);
                const double new_upper=std::min(1.+row[0],row[2]*new_factor);
                expansions[d]=std::max(0.,1.-old_upper/new_upper);
                const double cosine=row[9]/(row[2]*new_factor);
                if (!(cosine>=0. && cosine<=1.) || row[9]>new_upper) continue;
                p[4]=std::acos(cosine)*(180./pi);
                inclinations[d]=p[4];
                prior=std::log(old_upper/new_upper);
            } else if (geometry==2) {
                impacts[d]=row[2]*new_factor*std::cos(row[4]*rad);
            }
            const Orbit orbit(p,law);
            integrated_flux(count,times,exposures,exposure_count,orbit,row+10,samples,tolerance,
                            minimum_exposure,padding,model.data());
            for (std::size_t i=0;i<count;++i) model[i]*=std::exp(row[7]*airmass[i]);
            double baseline=row[8];
            if (!baseline_mode) {
                denominator.clear();numerator.clear();ratios.clear();
                for (std::size_t i=0;i<count;++i) if (baseline_mask[i] && std::isfinite(model[i]) && model[i]!=0.) {
                    denominator.push_back(baseline_weights[i]*(model[i]*model[i]));
                    numerator.push_back((baseline_weights[i]*data[i])*model[i]);
                    const double ratio=data[i]/model[i];
                    if (std::isfinite(ratio)) ratios.push_back(ratio);
                }
                const double denom=array_sum(denominator.data(),denominator.size());
                baseline=1.;
                if (!denominator.empty()) {
                    if (!std::isfinite(denom) || denom<=0.) {
                        if (!ratios.empty()) {
                            std::sort(ratios.begin(),ratios.end());
                            const std::size_t middle=ratios.size()/2;
                            baseline=ratios.size()%2?ratios[middle]:(ratios[middle-1]+ratios[middle])/2.;
                        }
                    } else baseline=array_sum(numerator.data(),numerator.size())/denom;
                }
                if (!std::isfinite(baseline)) baseline=1.;
            }
            for (std::size_t i=0;i<count;++i) {
                model[i]*=baseline;
                const double value=(data[i]-model[i])*inverse_error[i];
                residual[i]=value*value;
            }
            double score=-.5*array_sum(residual.data(),residual.size());
            if (duration_prior) {
                const double duration=orbit.duration();
                if (!std::isfinite(duration) || duration<=0.) continue;
                const double value=std::log(duration/expected_duration)/sigma_log_duration;
                score-=.5*(value*value);
            }
            if (std::isfinite(score)) scores[d]=score+prior;
        }
        return 0;
    } catch (const std::exception& exc) {
        if (error_size) std::snprintf(error,error_size,"%s",exc.what());
        return 1;
    } catch (...) {
        if (error_size) std::snprintf(error,error_size,"%s","Unknown native orbital replay error");
        return 2;
    }
}
