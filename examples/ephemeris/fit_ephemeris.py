import numpy as np
import statsmodels.api as sm
from exotic.api.ephemeris import ephemeris_fitter
import matplotlib.pyplot as plt

if __name__ == "__main__":
    Tc = np.array([  # measured mid-transit times
       2461656.170979  , 2460683.06352087, 2461312.22680483,
       2461840.72721957, 2461404.50457126, 2459614.88352437,
       2459967.2136158 , 2461250.70825625, 2460196.5097846 ,
       2459444.30884179, 2460297.17833986, 2460842.44956614,
       2461270.28536414, 2459447.10519731, 2459497.43952496,
       2460023.14013537, 2460445.37816892, 2461605.83844971,
       2459161.88372907, 2460878.80452072, 2460766.95256289,
       2460926.34200058, 2461712.09519907, 2460355.90036206,
       2459489.04986809, 2459379.99494715, 2459187.04825632,
       2460224.47242434, 2460979.4719767 , 2460406.23304265,
       2461700.91035602, 2460448.17416082, 2461611.43282057,
       2460132.19662872, 2460876.00842332, 2461446.45102159,
       2460395.04580418, 2460920.74932793, 2459463.87955256,
       2461756.83564536
    ])

    Tc_error = np.array([
       0.00086083, 0.00078861, 0.00086634, 0.00093534, 0.00075317,
       0.0008555 , 0.0007527 , 0.00078389, 0.00075229, 0.00076776,
       0.00042222, 0.00098135, 0.00075493, 0.00022053, 0.00038568,
       0.00070884, 0.00044426, 0.00043817, 0.00079945, 0.00084854,
       0.00053842, 0.00078082, 0.00079287, 0.00018737, 0.00018376,
       0.00062508, 0.00015116, 0.00034341, 0.00061096, 0.00083818,
       0.00076459, 0.00076027, 0.00081125, 0.00030851, 0.00015512,
       0.0007488 , 0.00027584, 0.00027871, 0.00080871, 0.00086118
    ])

    # labels for a legend
    labels = np.array([
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar',
        'TESS', 'TESS', 'EPW', 'ExoClock', 'Unistellar'
    ])

    P = 2.7962868  # orbital period for your target

    Tc_norm = Tc - Tc.min()  # normalize the data to the first observation
    # print(Tc_norm)
    orbit = np.rint(Tc_norm / P)  # number of orbits since first observation (rounded to nearest integer)
    # print(orbit)

    # make a n x 2 matrix with 1's in the first column and values of orbit in the second
    A = np.vstack([np.ones(len(Tc)), orbit]).T

    # perform the weighted least squares regression
    res = sm.WLS(Tc, A, weights=1.0 / Tc_error ** 2).fit()
    # use sm.WLS for weighted LS, sm.OLS for ordinary LS, or sm.GLS for general LS

    params = res.params  # retrieve the slope and intercept of the fit from res
    std_dev = np.sqrt(np.diagonal(res.normalized_cov_params))

    slope = params[1]
    slope_std_dev = std_dev[1]
    intercept = params[0]
    intercept_std_dev = std_dev[0]

    # 3 sigma clip based on residuals
    calculated = orbit * slope + intercept
    residuals = (Tc - calculated) / Tc_error
    mask = np.abs(residuals) < 3
    Tc = Tc[mask]
    Tc_error = Tc_error[mask]
    labels = labels[mask]

    # print(res.summary())
    # print("Params =",params)
    # print("Error matrix =",res.normalized_cov_params)
    # print("Standard Deviations =",std_dev)

    print("Weighted Linear Least Squares Solution")
    print("T0 =", intercept, "+-", intercept_std_dev)
    print("P =", slope, "+-", slope_std_dev)

    # min and max values to search between for fitting
    bounds = {
        'P': [P - 0.1, P + 0.1],  # orbital period
        'T0': [intercept - 0.1, intercept + 0.1]  # mid-transit time
    }

    # used to plot red overlay in O-C figure
    prior = {
        'P': [slope, slope_std_dev],  # value from WLS (replace with literature value)
        'T0': [intercept, intercept_std_dev]  # value from WLS (replace with literature value)
    }

    lf = ephemeris_fitter(Tc, Tc_error, bounds, prior=prior, labels=labels)

    lf.plot_triangle()
    plt.subplots_adjust(top=0.9, hspace=0.2, wspace=0.2)
    plt.savefig("posterior.png")
    plt.close()

    fig, ax = lf.plot_oc()
    plt.tight_layout()
    plt.savefig("oc.png")
    plt.show()
    plt.close()

    fig, ax = lf.plot_periodogram()
    plt.tight_layout()
    plt.savefig("periodogram.png")
    plt.show()
    plt.close()

