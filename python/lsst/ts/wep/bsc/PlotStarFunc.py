import numpy as np
import matplotlib.pyplot as plt

def plotRaDecl(wavefrontSensors, starMap, neighborStarMap, stddevSplit):
    """
    
    Plot stars in (Ra, Dec) and label the candidate stars and neighboring stars.
    
    Arguments:
        wavefrontSensors {[dict]} -- Collection of the wavefront sensor with corner coordinates.
        starMap {[dict]} -- Collection of star information.
        neighborStarMap {[dict]} -- Collection of information of neighboring stars on sensors.
        stddevSplit {[float} -- Value to decide the condition if the sensor crosses the RA=0.
    """

    # Handle the condition if the sensor across the RA=0
    acrossRA0 = False

    # Collect all sensorX postions of waverfront sensors
    allSensorX = np.array([])
    for detector, corners in wavefrontSensors.items():
        # Get the ra position on four corners
        sensorX = np.array([corners[0][0], corners[1][0], corners[2][0], corners[3][0]])
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
        __plotSingleRaDecl(wavefrontSensor, stars, neighboringStar, stddevSplit, acrossRA0)

    plt.xlabel("RA (degree)")
    plt.ylabel("Decl (degree)")
    plt.show()

def __plotSingleRaDecl(wavefrontSensor, stars, neighboringStar, stddevSplit, acrossRA0):
    """
    
    Plot stars in (Ra, Dec) and label the candidate stars and neighboring stars for singe CCD.
    
    Arguments:
        wavefrontSensor {[list]} -- List of wavefront sensor.
        stars {[StarData]} -- Star information.
        neighboringStar {[list]} -- Information of neighboring stars on sensors.
        stddevSplit {[float} -- Value to decide the condition if the sensor crosses the RA=0.
        acrossRA0 {[bool]} -- Sensors across RA=0 or not.
    """

    # Sensor corners in Ra, Decl
    sensorX = np.array([wavefrontSensor[0][0], wavefrontSensor[1][0], wavefrontSensor[2][0], 
                        wavefrontSensor[3][0]])
    sensorY = np.array([wavefrontSensor[0][1], wavefrontSensor[1][1], wavefrontSensor[2][1], 
                        wavefrontSensor[3][1]])

    # Star positions in Ra, Decl
    starX = np.array([stars.RA])
    starY = np.array([stars.Decl])

    # Map the neighboring stars and candidate stars
    neighborStarMapX = []
    neighborStarMapY = []

    candidateX = []
    candidateY = []

    # Get the candidate stars and neighboring stars information in star map
    for candidateStar, neighboringStars in neighboringStar.SimobjID.items():

        # Get the candidate stars
        candidateX = np.append(candidateX, neighboringStar.RaDecl[candidateStar][0]) 
        candidateY = np.append(candidateY, neighboringStar.RaDecl[candidateStar][1]) 

        # Get the neighboring stars
        for star in neighboringStars:
            neighborStarMapX = np.append(neighborStarMapX, neighboringStar.RaDecl[star][0])
            neighborStarMapY = np.append(neighborStarMapY, neighboringStar.RaDecl[star][1])

    # Shift the coordinates if sensors cross RA=0
    if (acrossRA0):
        # Shift Ra position for the plotting
        sensorX[np.where(sensorX>180)] = sensorX[np.where(sensorX>180)] - 360
        if (len(starX)):
            starX[np.where(starX>180)] = starX[np.where(starX>180)] - 360
        if (len(neighborStarMapX)):
            neighborStarMapX[np.where(neighborStarMapX>180)] = \
                                    neighborStarMapX[np.where(neighborStarMapX>180)] - 360
        if (len(candidateX)):
            candidateX[np.where(candidateX>180)] = candidateX[np.where(candidateX>180)] - 360

    # Rearrange points to plot the quadrilateral
    sensorX, sensorY = __getQuadrilateral(sensorX,sensorY)

    # Plot the figure
    plt.plot(sensorX, sensorY, "b")
    plt.plot(starX, starY, "bx")
    plt.plot(neighborStarMapX, neighborStarMapY, "go")
    plt.plot(candidateX, candidateY, "ro")

def plotPixel(stars, neighboringStar):
    """
    
    Plot stars in pixel and label the candidate stars and neighboring stars. 
    
    Arguments:
        stars {[StarData]} -- Star Information.
        neighboringStar {[list]} -- Information of neighboring stars on sensors.
    """

    # Star positions in pixel
    starX = np.array(stars.RAInPixel)     
    starY = np.array(stars.DeclInPixel)

    # Map the neighboring stars and candidate stars
    neighborStarMapX = []
    neighborStarMapY = []

    candidateX = []
    candidateY = []

    # Get the candidate stars and neighboring stars information in star map
    for candidateStar, neighboringStars in neighboringStar.SimobjID.items():

        # Get the candidate stars
        candidateX = np.append(candidateX, neighboringStar.RaDeclInPixel[candidateStar][0]) 
        candidateY = np.append(candidateY, neighboringStar.RaDeclInPixel[candidateStar][1]) 

        # Get the neighboring stars
        for star in neighboringStars:
            neighborStarMapX = np.append(neighborStarMapX, neighboringStar.RaDeclInPixel[star][0])
            neighborStarMapY = np.append(neighborStarMapY, neighboringStar.RaDeclInPixel[star][1])

    # Plot the figure
    plt.figure()
    plt.plot(starX, starY, "bx")
    plt.plot(neighborStarMapX, neighborStarMapY, "go")
    plt.plot(candidateX, candidateY, "ro")

    plt.xlabel("x-RA (pixel)")
    plt.ylabel("y-Decl (pixel)")

    plt.show()

def __getQuadrilateral(Xvalues, Yvalues):
    """
    
    Rearrange X and Y values and append the first point for the plotting.
    
    Arguments:
        Xvalues {[float]} -- X values.
        Yvalues {[float]} -- Y values.
    
    Returns:
        [float] -- Rearranged and appended X, Y values
    """

    # Check the length of X and Y should be 4
    if (len(Xvalues)!=4) or (len(Yvalues)!=4):
        print("The length of X and Y should be 4.")
    else:
        # Find the diagonal
        pair = [[2,3], [1,3], [1,2]]
        indexDiag = []
        for ii in range(1,4):
            # Get line equation by "y = m * x + b"
            m = (Yvalues[ii]-Yvalues[0])/(Xvalues[ii]-Xvalues[0])
            if np.isfinite(m):

                b = Yvalues[0] - m*Xvalues[0]
                
                dis1 = m*Xvalues[pair[ii-1][0]] + b - Yvalues[pair[ii-1][0]]
                dis2 = m*Xvalues[pair[ii-1][1]] + b - Yvalues[pair[ii-1][1]]
                
                # Find the diagonal
                if (dis1*dis2<0):
                    indexDiag = ii 

    if (indexDiag):
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