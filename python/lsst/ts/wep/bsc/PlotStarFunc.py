import numpy as np
import matplotlib.pyplot as plt


def plotRaDecl(wavefrontSensors, starMap, neighborStarMap, stddevSplit,
               saveFilePath=None):
    """Plot stars in (Ra, Dec) and label the candidate stars and neighboring
    stars.

    Parameters
    ----------
    wavefrontSensors : dict
        Collection of the wavefront sensor with corner coordinates. The
        dictionary key is the detector name. The item is a list contains the
        detector's corners information.
    starMap : dict
        Collection of star information. The dictionary key is the detector
        name. The item is the stars with the type of StarData.
    neighborStarMap : dict
        Collection of information of neighboring stars on sensors. The
        dictionary key is the detector name. The item is the neighboring stars
        with the type of NbrStar.
    stddevSplit : float
        Value to decide the condition if the sensor crosses the RA=0.
    saveFilePath : str, optional
        Save image to file path.
    """

    # Handle the condition if the sensor across the RA=0
    acrossRA0 = False

    # Collect all sensorX postions of waverfront sensors
    allSensorX = np.array([])
    for detector, corners in wavefrontSensors.items():
        # Get the ra position on four corners
        sensorX = np.array([corners[0][0], corners[1][0], corners[2][0],
                            corners[3][0]])
        allSensorX = np.append(allSensorX, sensorX)

    # Decide the sensors will cross RA=0 or not based on the standard deviation
    if np.std(allSensorX) >= stddevSplit:
        acrossRA0 = True

    # Plot the figure
    plt.figure()
    for detector in wavefrontSensors:
        wavefrontSensor = wavefrontSensors[detector]
        stars = starMap[detector]
        neighboringStar = neighborStarMap[detector]
        _plotSingleRaDecl(wavefrontSensor, stars, neighboringStar, acrossRA0)

    plt.xlabel("RA (degree)")
    plt.ylabel("Decl (degree)")

    # Save the figure and close or just show the iamge
    if (saveFilePath is not None):
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def _plotSingleRaDecl(wavefrontSensor, stars, neighboringStar, acrossRA0):
    """Plot stars in (Ra, Dec) and label the candidate stars and neighboring
    stars for singe CCD.

    Parameters
    ----------
    wavefrontSensor : list
        List of wavefront sensor corners.
    stars : StarData
        Star information.
    neighboringStar : NbrStar
        Information of neighboring stars on sensors.
    acrossRA0 : bool
        Sensors across RA=0 or not.
    """

    # Sensor corners in Ra, Decl
    sensorX = np.array([wavefrontSensor[0][0], wavefrontSensor[1][0],
                        wavefrontSensor[2][0], wavefrontSensor[3][0]])
    sensorY = np.array([wavefrontSensor[0][1], wavefrontSensor[1][1],
                        wavefrontSensor[2][1], wavefrontSensor[3][1]])

    # Star positions in Ra, Decl
    starX = stars.getRA()
    starY = stars.getDecl()

    # Map the candidate stars and neighboring stars
    candidateX = np.array([])
    candidateY = np.array([])

    neighborStarMapX = np.array([])
    neighborStarMapY = np.array([])

    # Get the candidate stars and neighboring stars information in star map
    for candidateStar, neighboringStars in neighboringStar.getId().items():

        # Get the candidate stars
        raDecl = neighboringStar.getRaDecl()
        candidateX = np.append(candidateX, raDecl[candidateStar][0])
        candidateY = np.append(candidateY, raDecl[candidateStar][1])

        # Get the neighboring stars
        for star in neighboringStars:
            neighborStarMapX = np.append(neighborStarMapX, raDecl[star][0])
            neighborStarMapY = np.append(neighborStarMapY, raDecl[star][1])

    # Shift the coordinates if sensors cross RA=0
    if (acrossRA0):
        # Shift Ra position for the plotting
        sensorX[np.where(sensorX > 180)] = sensorX[np.where(sensorX > 180)] - 360
        if (len(starX)):
            starX[np.where(starX > 180)] = starX[np.where(starX > 180)] - 360
        if (len(neighborStarMapX)):
            neighborStarMapX[np.where(neighborStarMapX > 180)] = \
                neighborStarMapX[np.where(neighborStarMapX > 180)] - 360
        if (len(candidateX)):
            candidateX[np.where(candidateX > 180)] = \
                candidateX[np.where(candidateX > 180)] - 360

    # Rearrange points to plot the quadrilateral
    sensorX, sensorY = _getQuadrilateral(sensorX, sensorY)

    # Plot the figure
    plt.plot(sensorX, sensorY, "b")
    plt.plot(starX, starY, "bx")
    plt.plot(neighborStarMapX, neighborStarMapY, "go")
    plt.plot(candidateX, candidateY, "ro")


def _getQuadrilateral(Xvalues, Yvalues):
    """Rearrange X and Y values and append the first point for the plotting.

    Parameters
    ----------
    Xvalues : numpy.ndarray
        X values.
    Yvalues : numpy.ndarray
        Y values.

    Returns
    -------
    numpy.ndarray
        Rearranged and appended X values
    numpy.ndarray
        Rearranged and appended Y values

    Raises
    -------
    ValueError
        The length of X and Y should be 4.
    """

    # Check the length of X and Y should be 4
    if (len(Xvalues) != 4) or (len(Yvalues) != 4):
        raise ValueError("The length of X and Y should be 4.")

    # Find the diagonal
    pair = [[2, 3], [1, 3], [1, 2]]
    indexDiag = -1
    for ii in range(1, 4):
        # Get line equation by "y = m * x + b"
        m = (Yvalues[ii] - Yvalues[0]) / (Xvalues[ii] - Xvalues[0])
        if np.isfinite(m):

            b = Yvalues[0] - m * Xvalues[0]

            dis1 = m * Xvalues[pair[ii - 1][0]] + b - \
                Yvalues[pair[ii - 1][0]]
            dis2 = m * Xvalues[pair[ii - 1][1]] + b - \
                Yvalues[pair[ii - 1][1]]

            # Find the diagonal
            if (dis1 * dis2 < 0):
                indexDiag = ii

    if (indexDiag != -1):
        Xtemp = Xvalues[2]
        Ytemp = Yvalues[2]

        # Put point to the diagonal by swaping
        Xvalues[2] = Xvalues[indexDiag]
        Yvalues[2] = Yvalues[indexDiag]

        Xvalues[indexDiag] = Xtemp
        Yvalues[indexDiag] = Ytemp

    Xvalues = np.append(Xvalues, Xvalues[0])
    Yvalues = np.append(Yvalues, Yvalues[0])

    return Xvalues, Yvalues


def plotStarInPixelOnDetector(stars, neighboringStar, xyDim=None,
                              saveFilePath=None):
    """Plot stars in pixel and label the candidate stars and neighboring stars.

    Parameters
    ----------
    stars : StarData
        Star information.
    neighboringStar : NbrStar
        Information of neighboring stars.
    xyDim : tuple, optional
        Dimensions of XY (dim X, dim Y). (the default is None.)
    saveFilePath : str, optional
        Save image to file path.
    """

    # Star positions in pixel
    starX = stars.getRaInPixel()
    starY = stars.getDeclInPixel()

    # Map the neighboring stars and candidate stars
    candidateX = np.array([])
    candidateY = np.array([])

    neighborStarMapX = np.array([])
    neighborStarMapY = np.array([])

    # Get the candidate stars and neighboring stars information in star map
    raDeclInPixel = neighboringStar.getRaDeclInPixel()
    for candidateStar, neighboringStars in neighboringStar.getId().items():

        # Get the candidate stars
        candidateX = np.append(candidateX, raDeclInPixel[candidateStar][0])
        candidateY = np.append(candidateY, raDeclInPixel[candidateStar][1])

        # Get the neighboring stars
        for star in neighboringStars:
            neighborStarMapX = np.append(neighborStarMapX,
                                         raDeclInPixel[star][0])
            neighborStarMapY = np.append(neighborStarMapY,
                                         raDeclInPixel[star][1])

    # Plot the figure
    plt.figure()
    plt.plot(starX, starY, "bx")
    plt.plot(neighborStarMapX, neighborStarMapY, "go")
    plt.plot(candidateX, candidateY, "ro")

    if (xyDim is not None):
        plt.xlim(0, xyDim[0])
        plt.ylim(0, xyDim[1])

    plt.xlabel("x-RA (pixel)")
    plt.ylabel("y-Decl (pixel)")

    # Save the figure and close or just show the iamge
    if (saveFilePath is not None):
        plt.savefig(saveFilePath, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


if __name__ == "__main__":
    pass
