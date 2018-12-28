import os
import subprocess
from enum import Enum

import lsst.ts.wep


class FilterType(Enum):
    U = 1
    G = 2
    R = 3
    I = 4
    Z = 5
    Y = 6


class BscDbType(Enum):
    LocalDb = 1


def getModulePath(module=lsst.ts.wep, startIdx=1, endIdx=-4):
    """Get the path of module.
    Parameters
    ----------
    module : str, optional
        Module name. (the default is lsst.ts.ofc.)
    startIdx : int, optional
        Start index. (the default is 1.)
    endIdx : int, optional
        End index. (the default is -4.)
    Returns
    -------
    str
        Directory path of module based on the start and end indexes.
    """

    # Get the path of module
    modulePathList = os.path.dirname(module.__file__).split(
                                os.sep)[int(startIdx):int(endIdx)]
    modulePath = os.path.join(os.sep, *modulePathList)

    return modulePath


def runProgram(command, binDir=None, argstring=None):
    """Run the program w/o arguments.

    Parameters
    ----------
    command : str
        Command of application.
    binDir : str, optional
        Directory of binary application. (the default is None.)
    argstring : str, optional
        Arguments of program. (the default is None.)

    Raises
    ------
    RuntimeError
        Error running of command.
    """

    # Directory of binary application
    if (binDir is not None):
        command = os.path.join(binDir, command)

    # Arguments for the program
    if (argstring is not None):
        command += (" " + argstring)

    # Call the program w/o arguments
    if (subprocess.call(command, shell=True) != 0):
        raise RuntimeError("Error running: %s" % command)


def writeFile(filePath, content):

    with open(filePath, "w") as file:
        file.write(content)


def readPhoSimSettingData(folderPath, fileName, atype):
    """Read the PhoSim setting data (segmentation or focal plane layout).

    Parameters
    ----------
    folderPath : str
        Path to folder.
    fileName : str
        File name ("segmentation.txt", "focalplanelayout.txt").
    atype : str
        Type of data to read ("readOutDim", "darkCurrent", "fieldCenter",
        "eulerRot").

    Returns
    -------
    dict
        Needed CCD data.
    
    Raises
    ------
    ValueError
        File can not be read.
    ValueError
        Type is not correct.
    """

    # Check the file name
    if fileName not in ("segmentation.txt", "focalplanelayout.txt"):
        raise ValueError("'%s' can not be read." % fileName)

    # Check the type
    if atype not in ("readOutDim", "darkCurrent", "fieldCenter", "eulerRot"):
        raise ValueError("'%s' can not be read." % atype)

    # Get the file path
    pathToFile = os.path.join(folderPath, fileName)

    # Amplifier list (only list the scientific ccd here)
    ampList = ["C00", "C01", "C02", "C03", "C04", "C05", "C06", "C07",
               "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17"]

    # Open the file to read
    ccdData = {}
    fid = open(pathToFile)
    for line in fid:
        line = line.strip()

        # Get each element
        lineElement = line.split()

        data = []
        # Analyze the sensor name to find the amplifier
        if (fileName == "segmentation.txt"):
            
            sensorNameStr = lineElement[0].split("_")
            if (len(sensorNameStr)==3):
                if sensorNameStr[2] in ampList:
                    name = lineElement[0]
                    # Get the segmentation in txt file
                    if (atype == "readOutDim"):
                        # parallel prescan, serial overscan, serial prescan,
                        # parallel overscan (pixel)
                        data = lineElement[15:19]
                    elif (atype == "darkCurrent"):
                        data = lineElement[13:15]

        elif (fileName == "focalplanelayout.txt"):

                # Analyze the sensor name to make sure this line of data is needed
                sensorNameStr = lineElement[0].split("_")
                if (len(sensorNameStr) == 2 or len(sensorNameStr) == 3):
                    if (atype == "fieldCenter"):
                        # Collect the field center:
                        # x position (microns), y position (microns), pixel size (microns)
                        # number of x pixels, number of y pixels
                        data = lineElement[1:6]
                    elif (atype == "eulerRot"):
                        # Collect the euler Rotation (degrees)
                        data = lineElement[10:13]

        # Collect the data
        if (data):
            ccdData.update({lineElement[0]:data})

    # Close the file
    fid.close()

    return ccdData


if __name__ == "__main__":
    pass
