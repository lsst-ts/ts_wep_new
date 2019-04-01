import numpy as np
import matplotlib.pyplot as plt

from lsst.ts.wep.cwfs.lib import cyMath

plt.switch_backend('Agg')


def ZernikeAnnularEval(Z, x, y, e, nMax=28):
    """Calculate the wavefront surface in the basis of annular Zernike
    polynomial.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of annular Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    e : float
        Obscuration value. It is 0.61 in LSST.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 28.)

    Returns
    -------
    numpy.ndarray
        Wavefront surface.
    """

    # Check the preconditions
    Z = _checkPrecondition(Z, x, y, nMax=nMax)

    # Calculate the wavefront
    return cyMath.ZernikeAnnularEval(Z, x.flatten(), y.flatten(), e).reshape(x.shape)


def ZernikeAnnularGrad(Z, x, y, e, axis, nMax=22):
    """Evaluate the gradident of annular Zernike polynomials in a certain
    direction.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of annular Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    e : float
        Obscuration value. It is 0.61 in LSST.
    axis : str
        It can be "dx", "dy", "dx2", "dy2", or "dxy".
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 22.)

    Returns
    -------
    numpy.ndarray
        Integration elements of gradient part in pupul x and y directions.
    """

    # Check the preconditions
    Z = _checkPrecondition(Z, x, y, nMax=nMax)

    # Calculate the integration elements
    return cyMath.ZernikeAnnularGrad(Z, x.flatten(), y.flatten(), e,
                                     axis).reshape(x.shape)


def ZernikeAnnularJacobian(Z, x, y, e, order, nMax=22):
    """Evaluate the Jacobian of annular Zernike polynomials in a certain order.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of annular Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    e : float
        Obscuration value. It is 0.61 in LSST.
    order : str
        Order of Jocobian Matrix. It can be "1st" or "2nd".
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 22.)

    Returns
    -------
    numpy.ndarray
        Jacobian elements in pupul x and y directions in a certain order.
    """

    # Check the preconditions
    Z = _checkPrecondition(Z, x, y, nMax=nMax)

    # Calculate the Jacobian
    return cyMath.ZernikeAnnularJacobian(Z, x.flatten(), y.flatten(), e,
                                         order).reshape(x.shape)


def ZernikeAnnularFit(S, x, y, numTerms, e, nMax=28):
    """Get the coefficients of annular Zernike polynomials by fitting the
    wavefront surface.

    Parameters
    ----------
    S : numpy.ndarray
        Wavefront surface to be fitted.
    x : numpy.ndarray
        Normalized x coordinate between -1 and 1 (pupil coordinate).
    y : numpy.ndarray
        Normalized y coordinate between -1 and 1 (pupil coordinate).
    numTerms : int
        Number of annular Zernike terms used in the fit.
    e : float
        Obscuration ratio of annular Zernikes.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 28.)

    Returns
    -------
    numpy.ndarray
        Coefficients of annular Zernike polynomials by the fitting.
    """

    # Check the dimensions of x and y are the same or not
    if (x.shape != y.shape):
        print("x & y are not the same size")

    # Get the value that is finite
    SFinite = S[:].copy()
    xFinite = x[:].copy()
    yFinite = y[:].copy()

    finiteIndex = np.isfinite(SFinite+xFinite+yFinite)

    SFinite = SFinite[finiteIndex]
    xFinite = xFinite[finiteIndex]
    yFinite = yFinite[finiteIndex]

    # Do the fitting
    H = np.zeros([len(SFinite), numTerms])

    for ii in range(numTerms):
        Z = np.zeros(numTerms)
        Z[ii] = 1
        H[:, ii] = ZernikeAnnularEval(Z, xFinite, yFinite, e, nMax=nMax)

    # Solve the equation: H*Z = S => Z = H^(-1)S
    Z = np.linalg.lstsq(H, S, rcond=None)[0]

    return Z


def ZernikeGrad(Z, x, y, axis, nMax=22):
    """Evaluate the gradident of Zernike polynomials in a certain axis.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    axis : str
        Integration direction. It can be "dx" or "dy".
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 22.)

    Returns
    -------
    numpy.ndarray
        Integration elements of gradient part in pupul x and y directions.
    """

    # Calculate the integration elements
    # Use obscuration (e) = 0 for standard Zernike polynomials
    return ZernikeAnnularGrad(Z, x, y, 0, axis, nMax=nMax)


def ZernikeJacobian(Z, x, y, order, nMax=22):
    """Evaluate the Jacobian of Zernike polynomials in a certain order.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    order : str
        Order of Jocobian Matrix. It can be "1st" or "2nd".
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 22.)

    Returns
    -------
    numpy.ndarray
        Jacobian elements in pupul x and y directions in a certain order.
    """

    # Calculate the Jacobian elements
    # Use obscuration (e) = 0 for standard Zernike polynomials
    return ZernikeAnnularJacobian(Z, x, y, 0, order, nMax=nMax)


def ZernikeEval(Z, x, y, nMax=28):
    """Calculate the wavefront surface in the basis of Zernike polynomial.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 28.)

    Returns
    -------
    numpy.ndarray
        Wavefront surface.
    """

    # Calculate the wavefront surface
    # Use obscuration (e) = 0 for standard Zernike polynomials
    return ZernikeAnnularEval(Z, x, y, 0, nMax=nMax)


def ZernikeFit(S, x, y, numTerms, nMax=28):
    """Get the coefficients of Zernike polynomials by fitting the wavefront
    surface.

    Parameters
    ----------
    S : numpy.ndarray
        Wavefront surface to be fitted.
    x : numpy.ndarray
        Normalized x coordinate between -1 and 1 (pupil coordinate).
    y : numpy.ndarray
        Normalized y coordinate between -1 and 1 (pupil coordinate).
    numTerms : int
        Number of Zernike terms used in the fit.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 28.)

    Returns
    -------
    numpy.ndarray
        Coefficients of Zernike polynomials by the fitting.
    """

    # Do the fitting to get coefficients of Zernike polynomials
    # Use obscuration (e) = 0 for standard Zernike polynomials
    return ZernikeAnnularFit(S, x, y, numTerms, 0, nMax=nMax)


def ZernikeMaskedFit(S, x, y, numTerms, mask, e, nMax=28):
    """Fit the wavefront surface on pupil (e.g. under the mask) to a linear
    combination of normal/ annular Zernike polynomials.

    Parameters
    ----------
    S : numpy.ndarray
        Wavefront surface to be fitted.
    x : numpy.ndarray
        Normalized x coordinate between -1 and 1 (pupil coordinate).
    y : numpy.ndarray
        Normalized y coordinate between -1 and 1 (pupil coordinate).
    numTerms : int
        Number of normal/ annular Zernike terms used in the fit.
    mask : numpy.ndarray[int]
        Mask used.
    e : float
         Obscuration ratio of annular Zernikes.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 28.)

    Returns
    -------
    numpy.ndarray
        Coefficients of normal/ annular Zernike polynomials by the fitting.
    """

    # Get S, x, y elements in mask
    j, i = np.nonzero(mask[:])
    S = S[i, j]
    x = x[i, j]
    y = y[i, j]

    # Calculate coefficients of normal/ spherical Zernike polynomials
    return ZernikeAnnularFit(S, x, y, numTerms, e, nMax=nMax)


def _checkPrecondition(Z, x, y, nMax=22):
    """Check the preconditions before the evaluation related to Zernike
    polynomials.

    Parameters
    ----------
    Z : numpy.ndarray
        Coefficient of Zernike polynomials.
    x : numpy.ndarray
        X coordinate on pupil plane.
    y : numpy.ndarray
        Y coordinate on pupil plane.
    nMax : int, optional
        Maximum number of Zernike terms. (the default is 22.)

    Returns
    -------
    numpy.ndarray
        Coefficient of Zernike polynomials.
    """

    # Check the dimensions of x and y are the same or not
    if (x.shape != y.shape):
        print("x & y are not the same size")
        exit()

    # Check the number of terms of annular Zernike polynomials
    if (len(Z) > nMax):
        print("Some Zernike related functions are not implemented with >%d terms." % nMax)
        return
    elif (len(Z) < nMax):
        # Put the higher order terms as zero to make sure nMax terms of
        # polynomials
        Z = np.hstack((Z, np.zeros(nMax-len(Z))))

    return Z


def padArray(inArray, dim):
    """Extend the boundary of image.

    For example, the input image is 120x120 matrix. This function will create
    an image such as 140x140 (dim x dim) matrix and put the input image in the
    center of new image.

    Parameters
    ----------
    inArray : numpy.ndarray
        Input central image.
    dim : int
        Dimension of new extended image.

    Returns
    -------
    numpy.ndarray
        Extended image from the dimension of inArray to dim x dim.

    Raises
    ------
    Exception
        Check the dimension of inArray is n by n or not.
    Exception
        Check the extending dimension is bigger than the dimension of inArray
        or not.
    """

    # Check the conditions
    m, n = inArray.shape
    if (m != n):
        raise Exception("padArray: array is not square.")

    if (m > dim):
        raise Exception("padArray: array is larger than dimension.")

    # Extend the boundary of image by creating a bigger matrix and putting the
    # input image in the center
    out = np.zeros([dim, dim])
    ii = int(np.floor((dim-m)/2))

    # Put the original image in the center of extended image
    out[ii:ii+m, ii:ii+m] = inArray

    return out


def extractArray(inArray, dim):
    """Extract the central image.

    For example, the input image is a 140x140 matrix. This function will
    extract the central matrix with the dimension of 120x120 (dim x dim).

    Parameters
    ----------
    inArray : numpy.ndarray
        Input image.
    dim : int
        Dimension of extracted image.

    Returns
    -------
    numpy.ndarray
        Extracted central image from the dimension of inArray to dim x dim.

    Raises
    ------
    Exception
        Check the dimension of inArray is n by n or not.
    Exception
        Check the extracted dimension is smaller than the dimension of inArray
        or not.
    """

    # Check the conditions
    m, n = inArray.shape
    if (m != n):
        raise Exception("extractArray: array is not square.")

    if (m < dim):
        raise Exception("extractArray: array is smaller than dimension")

    # Calculate the begining index to extract the central image
    ii = int(np.floor((m-dim)/2))

    # Extract the cetral image
    out = inArray[ii:ii+dim, ii:ii+dim]

    return out


def plotImage(imageDonut, title=None, mask=None, show=True, fitParameters=[],
              saveFilePath=None):
    """Show the wavefront (donut) image.

    Parameters
    ----------
    imageDonut : numpy.ndarray[float]
        Input donut Image.
    title : str, optional
        Title of image. (the default is None.)
    mask : numpy.ndarray[int], optional
        Mask. (the default is None.)
    show : bool, optional
        Show the figure or not. (the default is True.)
    fitParameters : list, optional
        Fitting parameter of circle (center position and radius) (the default
        is [].)
    saveFilePath : str, optional
        File path to save the image. (the default is None.)
    """

    if mask is not None:
        image = np.where(mask == 0, np.nan, imageDonut)
    else:
        image = imageDonut

    plt.figure()

    # Plot the fitted circle
    if (fitParameters):
        if (len(fitParameters) == 3):
            theta = np.linspace(0, 2*np.pi, 101)
            x = fitParameters[0] + fitParameters[2]*np.cos(theta)
            y = fitParameters[1] + fitParameters[2]*np.sin(theta)
            plt.plot(x, y, "b")
        else:
            print("fitParameters should have 3 elements.")

    # Plot the wavefront image
    plt.imshow(image, origin="lower")
    plt.colorbar()
    if title:
        plt.title(title)

    if (saveFilePath is not None):
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()

    if show:
        plt.show()


def getConfigValue(configFilePath, varName, index=1):
    """Get the value of certain variable defined in the configuration file.

    Parameters
    ----------
    configFilePath : str
        Configuration file path.
    varName : str
        Name of variable.
    index : int, optional
        Index of value. (the default is 1.)

    Returns
    -------
    float, int, or str
        Variable value.
    """

    # Read the file
    fid = open(configFilePath)

    # Search for the value of certain variable
    value = None
    for line in fid:

        # Strip the line
        line = line.strip()

        # Get the element of line
        lineArray = line.split()

        # Get the value
        if (len(lineArray) > 0) and (lineArray[0] == varName):
            value = lineArray[index]
            break

    # Close the file
    fid.close()

    # Change the value type if necessary
    try:
        value = float(value)
        if (value == int(value)):
            value = int(value)
    except ValueError:
        pass

    return value


if __name__ == "__main__":
    pass
