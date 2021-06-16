import numpy as np
from dataclasses import dataclass

from scipy.optimize import least_squares
from scipy.ndimage import median_filter
from scipy.interpolate import RectBivariateSpline

from photutils import aperture_photometry
from photutils import CircularAperture

from astropy.io import fits

from skimage.registration import phase_cross_correlation

try:
    from api.elca import lc_fitter, binner, transit
except ImportError:
    from .api.elca import lc_fitter, binner, transit


def mesh_box_prev(pos, box):
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0] - box, pos[0] + box + 1)
    y = np.arange(pos[1] - box, pos[1] + box + 1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


# Method fits a 2D gaussian function that matches the star_psf to the star image and returns its pixel coordinates
def fit_centroid_prev(data, pos, init=None, lossfn='linear', box=10):
    if init is None:  # if init is none, then set the values
        init = [-1, 5, 5, 0]
    else:  # remove rotation
        init = np.delete(init, 3)

    # estimate the amplitude and centroid
    if init[0] == -1:
        # subarray of data around star
        xv, yv = mesh_box_prev(pos, box)

        # amplitude guess
        init[0] = np.max(data[yv, xv])

        # weighted sum to estimate center
        wx = np.sum(np.unique(xv) * data[yv, xv].sum(0)) / np.sum(data[yv, xv].sum(0))
        wy = np.sum(np.unique(yv) * data[yv, xv].sum(1)) / np.sum(data[yv, xv].sum(1))
        pos = [wx, wy]
        # estimate std by calculation of FWHM
        x, y = data[yv, xv].sum(0), data[yv, xv].sum(1)
        init[1] = estimate_sigma_prev(x)
        init[2] = estimate_sigma_prev(y)

        # Background Estimate
        # compute the average from 1/4 of the lowest values in the background
        init[3] = np.mean(np.sort(data[yv, xv].flatten())[:int(data[yv, xv].flatten().shape[0] * 0.25)])
    # print('init priors for centroid:',init)
    # print('init2:',init)

    # recenter data on weighted average of light (peak amplitude)
    xv, yv = mesh_box_prev(pos, box)

    # pars = x,y, a,sigx,sigy, rotate
    def fcn2min(pars):
        model = star_psf_prev(xv, yv, *pars)
        return (data[yv, xv] - model).flatten()  # method for LS
        # return np.sum( (data[yv,xv]-model)**2 ) # method for minimize

    lo = [pos[0] - box, pos[1] - box, 0, 1, 1, 0]
    up = [pos[0] + box, pos[1] + box, 100000, 40, 40, np.max(data[yv, xv])]
    res = least_squares(fcn2min, x0=[*pos, *init], bounds=[lo, up], loss=lossfn, jac='3-point')

    return np.insert(res.x, 5, 0)


# defines the star point spread function as a 2D Gaussian
def star_psf_prev(x, y, x0, y0, a, sigx, sigy, b):
    gaus = a * np.exp(-(x - x0) ** 2 / (2 * sigx ** 2)) * np.exp(-(y - y0) ** 2 / (2 * sigy ** 2)) + b
    return gaus


# Method uses the Full Width Half Max to estimate the standard deviation of the star's psf
def estimate_sigma_prev(x, maxidx=-1):
    if maxidx == -1:
        maxidx = np.argmax(x)
    lower = np.abs(x - 0.5 * np.max(x))[:maxidx].argmin()
    upper = np.abs(x - 0.5 * np.max(x))[maxidx:].argmin() + maxidx
    FWHM = upper - lower
    return FWHM / (2 * np.sqrt(2 * np.log(2)))


def sigma_clip_prev(ogdata, sigma=3, dt=20):
    mdata = median_filter(ogdata, dt)
    res = ogdata - mdata
    std = np.nanmedian([np.nanstd(np.random.choice(res,50)) for i in range(100)])
    #std = np.nanstd(res) # biased from large outliers
    return np.abs(res) > sigma*std


# Method calculates the average flux of the background
def skybg_phot_prev(x0, y0, data, r=25, dr=5, samp=3, debug=False):
    # determine img indexes for aperture region
    xv, yv = mesh_box_prev([x0, y0], int(np.round(r + dr)))

    # derive indexs on a higher resolution grid and create aperture mask
    px, py, mask = sky_annulus_prev(x0, y0, r=r, samp=xv.shape[0] * samp)

    # interpolate original data onto higher resolution grid
    subdata = data[yv, xv]
    model = RectBivariateSpline(np.unique(xv), np.unique(yv), subdata)

    # evaluate data on highres grid
    pz = model.ev(px, py)

    # zero out pixels larger than radius
    pz[~mask] = 0
    pz[pz < 0] = 0

    quarterMask = pz < np.percentile(pz[mask], 50)
    pz[~quarterMask] = 0

    # scale area back to original grid, total flux in sky annulus
    parea = pz.sum() * np.diff(px).mean() * np.diff(py[:, 0]).mean()

    if debug:
        print('mask area=', mask.sum() * np.diff(px).mean() * np.diff(py[:, 0]).mean())
        print('true area=', 2 * np.pi * r * dr)
        print('subdata flux=', subdata.sum())
        print('bg phot flux=', parea)
        import pdb
        pdb.set_trace()

    # return bg value per pixel
    bgmask = mask & quarterMask
    avgBackground = pz.sum() / bgmask.sum()
    return avgBackground


# Method defines the annulus used to do a background subtraction
def sky_annulus_prev(x0, y0, r=25, dr=5, samp=10):
    xv, yv = mesh_box_prev([x0, y0], r + dr + 1, npts=samp)
    rv = ((xv - x0) ** 2 + (yv - y0) ** 2) ** 0.5
    mask = (rv > r) & (rv < (r + dr))  # sky annulus mask
    return xv, yv, mask


# Method that computes and returns the total flux of the given star
# calls the phot function for flux calculation which includes the background subtraction
def getFlux_prev(photoData, xPix, yPix, apertureRad, annulusRad):
    bgSub, totalFlux = phot_prev(xPix, yPix, photoData, r=apertureRad, dr=annulusRad, debug=False, bgsub=True)
    return bgSub, totalFlux  # return the total flux for the given star in the one image


# Method calculates the flux of the star (uses the skybg_phot method to do backgorund sub)
def phot_prev(x0, y0, data, r=25, dr=5, samp=5, debug=False, bgsub=True):
    if bgsub:
        # get the bg flux per pixel
        bgflux = skybg_phot_prev(x0, y0, data, r, dr, samp)
    else:
        bgflux = 0

    positions = [(x0, y0)]
    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(data - bgflux, apertures)
    # print(phot_table[0][3])
    bgSubbed = phot_table[0][3]

    rawPhot_Table = aperture_photometry(data, apertures)
    raw = rawPhot_Table[0][3]
    return bgSubbed, raw


# Method that gets and returns the current phase of the target
def getPhase_prev(curTime, pPeriod, tMid):
    phase = (curTime - tMid - 0.5 * pPeriod) / pPeriod % 1
    return phase - 0.5


@dataclass
class PhotometryOutput:
    finXTargCent: list
    finYTargCent: list
    finXRefCent: list
    finYRefCent: list
    goodFluxes: list
    nonBJDTimes: list
    nonBJDPhases: list
    goodAirmasses: list
    goodTargets: list
    goodReferences: list
    goodTUnc: list
    goodRUnc: list
    goodNormUnc: list
    goodResids: np.array
    bestlmfit: None
    bestCompStar: int
    comp_coords: list


def photometry(input_files, compStarList, targx, targy, pDict, ld0, ld1, ld2, ld3):
    minAperture = max(1, int(2 * max(targsigX, targsigY)))
    maxAperture = int(5 * max(targsigX, targsigY) + 1)
    minAnnulus = 2
    maxAnnulus = 5
    distFC = 10
    firstImageData = input_files[0]

    # fit centroids for first image to determine priors to be used later
    for compCounter in range(0, len(compStarList)):
        print('\n\n***************************************************************')
        print('Determining Optimal Aperture and Annulus Size for Comp Star #' + str(compCounter + 1))
        print('***************************************************************')

        # #just in case comp star drifted off and timeSortedNames had to be altered, reset it for the new comp star
        # timeSortedNames = tsnCopy

        UIprevRPX, UIprevRPY = compStarList[compCounter]
        print('Target X: ' + str(round(targx)) + ' Target Y: ' + str(round(targy)))
        refx, refy, refamplitude, refsigX, refsigY, retrot, refoff = fit_centroid_prev(firstImageData, [UIprevRPX, UIprevRPY])
        print('Comparison X: ' + str(round(refx)) + ' Comparison Y: ' + str(round(refy)) + '\n')

        # determines the aperture and annulus combinations to iterate through based on the sigmas of the LM fit
        aperture_min = int(3 * np.nanmax([targsigX, targsigY]))
        aperture_max = int(5 * np.nanmax([targsigX, targsigY]))

        # Run through only 5 different aperture sizes, all interger pixel values
        aperture_step = np.nanmax([1, (aperture_max + 1 - aperture_min) // 5])  # forces step size to be at least 1
        aperture_sizes = np.arange(aperture_min, aperture_max + 1, aperture_step)
        if aperture_min <= 1:
            aperture_sizes = np.arange(1, 10, 2)

        # single annulus size
        annulus_sizes = [5]

        target_fits = {}
        ref_fits = {}
        reg_trans = {}

        for apertureR in aperture_sizes:  # aperture loop
            for annulusR in annulus_sizes:  # annulus loop # no need
                # fileNumber = 1
                print('Testing Comparison Star #' + str(compCounter + 1) + ' with a ' + str(
                    apertureR) + ' pixel aperture and a ' + str(annulusR) + ' pixel annulus.')
                for fileNumber, fileName in enumerate(input_files):
                    hdul = fits.open(name=fileName, memmap=False, cache=False, lazy_load_hdus=False,
                                     ignore_missing_end=True)

                    extension = 0
                    image_header = hdul[extension].header
                    while image_header["NAXIS"] == 0:
                        extension += 1
                        image_header = hdul[extension].header

                    # IMAGES
                    imageData = hdul[extension].data

                    # Find the target star in the image and get its pixel coordinates if it is the first file
                    if fileNumber == 0:
                        # Initializing the star location guess as the user inputted pixel coordinates
                        prevTPX, prevTPY, prevRPX, prevRPY = UIprevTPX, UIprevTPY, UIprevRPX, UIprevRPY  # 398, 275, 419, 203
                        prevTSigX, prevTSigY, prevRSigX, prevRSigY = targsigX, targsigY, refsigX, refsigY
                        prevImageData = imageData  # no shift should be registered

                    # ------ CENTROID FITTING ----------------------------------------

                    # corrects for any image shifts that result from a tracking slip
                    # shift, error, diffphase = phase_cross_correlation(ImageData, imageData)
                    if fileNumber in reg_trans.keys():
                        shift, error, diffphase = reg_trans[fileNumber]
                    else:
                        shift, error, diffphase = phase_cross_correlation(prevImageData, imageData)
                        reg_trans[fileNumber] = [shift, error, diffphase]

                    xShift = shift[1]
                    yShift = shift[0]

                    prevTPX = prevTPX - xShift
                    prevTPY = prevTPY - yShift
                    prevRPX = prevRPX - xShift
                    prevRPY = prevRPY - yShift

                    # set target search area
                    txmin = int(prevTPX) - distFC  # left
                    txmax = int(prevTPX) + distFC  # right
                    tymin = int(prevTPY) - distFC  # top
                    tymax = int(prevTPY) + distFC  # bottom

                    # boolean that represents if either the target or comp star gets too close to the detector
                    driftBool = False

                    # check if your target star is too close to the edge of the detector
                    if txmin <= 0 or tymin <= 0 or txmax >= len(imageData) or tymax >= len(imageData[0]):
                        print('*************************************************************************************')
                        print('WARNING: In image ' + str(
                            fileNumber) + ', your target star has drifted too close to the edge of the detector.')
                        # tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                        print('All the remaining images after image #' + str(fileNumber - 1) + ' will be ignored')
                        driftBool = True

                        # mask off the rest of timeSortedNames and then ignore the rest of the procedure until

                    # Set reference search area
                    rxmin = int(prevRPX) - distFC  # left
                    rxmax = int(prevRPX) + distFC  # right
                    rymin = int(prevRPY) - distFC  # top
                    rymax = int(prevRPY) + distFC  # bottom

                    # check if the reference is too close to the edge of the detector
                    if rxmin <= 0 or rymin <= 0 or rxmax >= len(imageData[0]) or rymax >= len(imageData):
                        print('*************************************************************************************')
                        print('WARNING: In image ' + str(
                            fileNumber) + ', your reference star has drifted too close to the edge of the detector.')
                        # tooClose = int(input('Enter "1" to pick a new comparison star or enter "2" to continue using the same comp star, with the images with all the remaining images ignored \n'))
                        print('All the remaining images after image #' + str(
                            fileNumber - 1) + ' will be ignored for this comparison star')
                        print('*************************************************************************************')
                        driftBool = True

                    # if the star isn't too close, then proceed as normal
                    if not driftBool:
                        targSearchA = imageData[tymin:tymax, txmin:txmax]
                        refSearchA = imageData[rymin:rymax, rxmin:rxmax]

                        targPos = [prevTPX, prevTPY]

                        # get minimum background value bigger than 0
                        targImFlat = np.sort(np.array(targSearchA).ravel())

                        # Initialize the variable
                        tGuessBkg = 0
                        for el in targImFlat:
                            if el > 0:
                                tGuessBkg = el
                                break

                        refImFlat = np.sort(np.array(refSearchA).ravel())
                        for rel in refImFlat:
                            if rel > 0:
                                rGuessBkg = rel
                                break

                        # Guess at Gaussian Parameters and feed them in to help gaussian fitter

                        tGuessAmp = targSearchA.max() - tGuessBkg
                        if tGuessAmp < 0:
                            print('Error: the Darks have a higher pixel counts than the image itself')
                        myPriors = [tGuessAmp, prevTSigX, prevTSigY, 0, tGuessBkg]

                        # tx, ty, tamplitude, tsigX, tsigY, toff = fit_centroid(imageData, targPos,
                        #                                                     init=myPriors, box=distFC)
                        if fileNumber in target_fits.keys():
                            tx, ty, tamplitude, tsigX, tsigY, toff = target_fits[fileNumber]
                        else:
                            tx, ty, tamplitude, tsigX, tsigY, trot, toff = fit_centroid_prev(imageData, targPos,
                                                                                             init=myPriors)
                            target_fits[fileNumber] = [tx, ty, tamplitude, tsigX, tsigY, toff]

                        currTPX = tx
                        currTPY = ty

                        # append to list of target centroid positions for later plotting
                        xTargCent.append(currTPX)
                        yTargCent.append(currTPY)

                        rGuessAmp = refSearchA.max() - rGuessBkg
                        myRefPriors = [rGuessAmp, prevRSigX, prevRSigY, 0, rGuessBkg]
                        # rx, ry, ramplitude, rsigX, rsigY, roff = fit_centroid(imageData, [prevRPX, prevRPY],
                        # init=myRefPriors, box=distFC)
                        if fileNumber in ref_fits.keys():
                            rx, ry, ramplitude, rsigX, rsigY, roff = ref_fits[fileNumber]
                        else:
                            rx, ry, ramplitude, rsigX, rsigY, rrot, roff = fit_centroid_prev(imageData, [prevRPX, prevRPY],
                                                                                             init=myRefPriors)
                            ref_fits[fileNumber] = [rx, ry, ramplitude, rsigX, rsigY, roff]
                        currRPX = rx
                        currRPY = ry

                        # append to list of reference centroid positions for later plotting
                        xRefCent.append(currRPX)
                        yRefCent.append(currRPY)

                        if tamplitude < 0 or tsigX < 0 or tsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                            print('Could not fit 2D gaussian to Target for File Number' + str(fileNumber))

                        elif ramplitude < 0 or rsigX < 0 or rsigY < 0:  # gets rid of negative amplitude values that indicate it couldn't fit gaussian
                            print('Could not fit 2D gaussian to Comparison Star for File Number' + str(fileNumber))

                        else:
                            # ------FLUX CALCULATION WITH BACKGROUND SUBTRACTION----------------------------------

                            # gets the flux value of the target star and subtracts the background light
                            tFluxVal, tTotCts = getFlux_prev(imageData, currTPX, currTPY, apertureR, annulusR)
                            # FIXME centroid position is way off from user input for star

                            targetFluxVals.append(tFluxVal)  # adds tFluxVal to the total list of flux values of target star
                            targUncertanties.append(
                                np.sqrt(tFluxVal))  # uncertanty on each point is the sqrt of the total counts

                            # gets the flux value of the reference star and subracts the background light
                            rFluxVal, rTotCts = getFlux_prev(imageData, currRPX, currRPY, apertureR, annulusR)

                            referenceFluxVals.append(
                                rFluxVal)  # adds rFluxVal to the total list of flux values of reference star
                            refUncertanties.append(np.sqrt(rFluxVal))

                            # # TIME
                            # currTime = getJulianTime(hDul)
                            # timesListed.append(currTime)

                            # ORBITAL PHASE
                            currentPhase = getPhase_prev(currTime, pDict['pPer'], pDict['midT'])
                            phasesList.append(currentPhase)  # adds to list of phases

                            # # AIRMASS
                            # airMass = getAirMass(hDul)  # gets the airmass at the time the image was taken
                            # airMassList.append(airMass)  # adds that airmass value to the list of airmasses

                            # UPDATE PIXEL COORDINATES and SIGMAS
                            # target
                            prevTPX = currTPX
                            prevTPY = currTPY
                            prevTSigX = tsigX
                            prevTSigY = tsigY
                            # reference
                            prevRPX = currRPX
                            prevRPY = currRPY
                            prevRSigX = rsigX
                            prevTSigY = rsigY

                        # UPDATE FILE COUNT
                        prevImageData = imageData
                        # fileNumber = fileNumber + 1
                        # hDul.close()  # close the stream

                    # otherwise, mask off the rest of the files from time sorted names including the current one
                    else:
                        print("\nFiltering data to account for drifting target.")
                        # timeSortedNames = timeSortedNames[:fileNumber]

                        # TIME
                        timesListed = timesListed[:fileNumber]

                        # AIRMASS
                        airMassList = airMassList[:fileNumber]

                        # ALL IMAGES
                        sortedallImageData = sortedallImageData[:fileNumber]

                        boollist = boollist[:fileNumber]

                        break

                # EXIT THE FILE LOOP

                # NORMALIZE BY REF STAR
                # Convert the raw flux values to arrays and then divide them to get the normalized flux data
                rawFinalFluxData = np.array(targetFluxVals)[boollist]

                # Convert Everything to numpy Arrays
                arrayFinalFlux = np.array(rawFinalFluxData)  # finalFluxData
                arrayTargets = np.array(targetFluxVals)[boollist]  # finalFluxData
                arrayTimes = np.array(timesListed)[boollist]
                arrayPhases = np.array(phasesList)[boollist]
                arrayTargets = np.array(targetFluxVals)[boollist]
                arrayReferences = np.array(referenceFluxVals)[boollist]
                arrayAirmass = np.array(airMassList)[boollist]
                arrayTUnc = np.array(targUncertanties)[boollist]
                arrayRUnc = np.array(refUncertanties)[boollist]

                arrayNormUnc = rawFinalFluxData ** 0.5

                # Execute sigma_clip
                try:
                    filtered_data = sigma_clip_prev(arrayFinalFlux, sigma=3)
                except TypeError:
                    filtered_data = sigma_clip_prev(arrayFinalFlux, sigma=3)

                # -----LM LIGHTCURVE FIT--------------------------------------

                prior = {
                    'rprs': pDict['rprs'],  # Rp/Rs
                    'ars': pDict['aRs'],  # a/Rs
                    'per': pDict['pPer'],  # Period [day]
                    'inc': pDict['inc'],  # Inclination [deg]
                    'u0': ld0[0], 'u1': ld1[0], 'u2': ld2[0], 'u3': ld3[0],  # limb darkening (nonlinear)
                    'ecc': pDict['ecc'],  # Eccentricity
                    'omega': 0,  # Arg of periastron
                    'tmid': pDict['midT'],  # time of mid transit [day]
                    'a1': arrayFinalFlux.mean(),  # max() - arrayFinalFlux.min(), #mid Flux
                    'a2': 0,  # Flux lower bound
                }

                phase = (arrayTimes[~filtered_data] - prior['tmid']) / prior['per']
                prior['tmid'] = pDict['midT'] + np.floor(phase).max() * prior['per']
                upper = pDict['midT'] + 25 * pDict['midTUnc'] + np.floor(phase).max() * (
                            pDict['pPer'] + 25 * pDict['pPerUnc'])
                lower = pDict['midT'] - 25 * pDict['midTUnc'] + np.floor(phase).max() * (
                            pDict['pPer'] - 25 * pDict['pPerUnc'])

                if np.floor(phase).max() - np.floor(phase).min() == 0:
                    print('Estimated mid-transit not in observation range (check priors or observation time)')
                    print('start:', arrayTimes[~filtered_data].min())
                    print('  end:', arrayTimes[~filtered_data].max())
                    print('prior:', prior['tmid'])

                mybounds = {
                    'rprs': [pDict['rprs'] - 3 * pDict['rprsUnc'], pDict['rprs'] + 3 * pDict['rprsUnc']],
                    'tmid': [max(lower, arrayTimes[~filtered_data].min()), min(arrayTimes[~filtered_data].max(), upper)],
                    'ars': [pDict['aRs'] - 5 * pDict['aRsUnc'], pDict['aRs'] + 5 * pDict['aRsUnc']],

                    'a1': [0, 3 * max(arrayFinalFlux[~filtered_data])],
                    'a2': [-3, 3],
                    # 'a3':[0, max(arrayFinalFlux[~filtered_data])]
                }

                myfit = lc_fitter(
                    arrayTimes[~filtered_data],
                    arrayFinalFlux[~filtered_data],
                    arrayNormUnc[~filtered_data],
                    arrayAirmass[~filtered_data],
                    prior,
                    mybounds,
                    mode='lm'
                )

                for k in myfit.bounds.keys():
                    print("  {}: {:.6f}".format(k, myfit.parameters[k]))  # , myfit.errors[k]))

                print('The Residual Standard Deviation is: %' + str(
                    round(100 * myfit.residuals.std() / np.median(myfit.data), 6)))
                print('The Mean Squared Error is: ' + str(round(np.sum(myfit.residuals ** 2), 6)) + '\n')
                if minSTD > myfit.residuals.std():  # If the standard deviation is less than the previous min
                    bestCompStar = compCounter + 1
                    comp_cords = compStarList[compCounter]
                    minSTD = myfit.residuals.std()  # set the minimum standard deviation to that

                    arrayNormUnc = arrayNormUnc * np.sqrt(
                        myfit.chi2 / myfit.data.shape[0])  # scale errorbars by sqrt(rchi2)
                    minAnnulus = annulusR  # then set min aperature and annulus to those values
                    minAperture = apertureR
                    # gets the centroid trace plots to ensure tracking is working
                    finXTargCentArray = np.array(xTargCent)[boollist]
                    finYTargCentArray = np.array(yTargCent)[boollist]
                    finXRefCentArray = np.array(xRefCent)[boollist]
                    finYRefCentArray = np.array(yRefCent)[boollist]

                    # APPLY DATA FILTER
                    # apply data filter sets the lists we want to print to correspond to the optimal aperature
                    finXTargCent = finXTargCentArray[~filtered_data]
                    finYTargCent = finYTargCentArray[~filtered_data]
                    finXRefCent = finXRefCentArray[~filtered_data]
                    finYRefCent = finYRefCentArray[~filtered_data]
                    # sets the lists we want to print to correspond to the optimal aperature
                    goodFluxes = arrayFinalFlux[~filtered_data]
                    nonBJDTimes = arrayTimes[~filtered_data]
                    nonBJDPhases = arrayPhases[~filtered_data]
                    goodAirmasses = arrayAirmass[~filtered_data]
                    goodTargets = arrayTargets[~filtered_data]
                    goodReferences = arrayReferences[~filtered_data]
                    goodTUnc = arrayTUnc[~filtered_data]
                    goodRUnc = arrayRUnc[~filtered_data]
                    goodNormUnc = arrayNormUnc[~filtered_data]
                    goodResids = myfit.residuals
                    bestlmfit = myfit

                # Reinitialize the the arrays to be empty
                # airMassList = []
                phasesList = []
                # timesListed = []
                targetFluxVals = []
                referenceFluxVals = []
                targUncertanties = []
                refUncertanties = []
                normUncertainties = []
                xTargCent = []
                yTargCent = []
                xRefCent = []
                yRefCent = []

            # Exit aperture loop
        # Exit annulus loop
    # Exit the Comp Stars Loop

    print('\n*********************************************')
    print('Best Comparison Star: #' + str(bestCompStar))
    print('Minimum Residual Scatter: ' + str(round(minSTD / np.median(goodFluxes) * 100, 4)) + '%')
    print('Optimal Aperture: ' + str(minAperture))
    print('Optimal Annulus: ' + str(minAnnulus))
    print('********************************************\n')

    return PhotometryOutput(finXTargCent, finYTargCent, finXRefCent, finYRefCent, goodFluxes, nonBJDTimes, nonBJDPhases,
                            goodAirmasses, goodTargets, goodReferences, goodTUnc, goodRUnc, goodNormUnc, goodResids,
                            bestlmfit, bestCompStar, comp_cords)


if __name__ == '__main__':
    targsigX = 0
    targsigY = 0
    UIprevTPX = 0
    UIprevTPY = 0
    currTime = 0
