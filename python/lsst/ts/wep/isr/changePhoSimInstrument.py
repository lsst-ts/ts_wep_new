import os, shutil, unittest
from lsst.ts.wep.Utility import getModulePath

def updateData(folderPath, fileName, newData, atype, sameCcdList=None):
    """

    Update the setting such as the readout dimension (parallel prescan, serial overscan,
    serial prescan, parallel overscan (pixel)) or dark current.

    Arguments:
        folderPath {[string]} -- Path to folder.
        fileName {[string]} -- File name.
        newData {[string]} -- New data to replace the original one.
        atype {[string]} -- Type of data to change ("readOutDim" or "darkCurrent").

    Keyword Arguments:
        sameCcdList {[list]} -- List of ccd that keeps the same value (default: {None}).
    """

    # Get he path to file
    pathToFile = os.path.join(folderPath, fileName)

    # Replace the segmentation
    with open(pathToFile, "r") as input_file, open(fileName, "w") as output_file:
        for line in input_file:

            # Update the details in file (segmentation.txt)
            if (not line.startswith("#") and len(line)>100):
                # Get each element
                lineElement = line.split()

                # Update the values
                if sameCcdList is not None:
                    # Get the ccd name
                    sensorNameStr = lineElement[0].split("_")
                    ccdName = "_".join(sensorNameStr[0:2])

                    # Only update the ccd that is not in the list
                    if ccdName not in sameCcdList:
                        lineElement = __changeValue(lineElement, newData, atype)

                elif sameCcdList is None:
                    lineElement = __changeValue(lineElement, newData, atype)

                # If the type is wrong, return and finish this funtion.
                if lineElement is None:
                    return

                # Write the line with new data in the same format of the original one
                if (len(lineElement[1]) == 1):
                    lineFormat = "{}     {}"
                else:
                    lineFormat = "{}  {}"

                if (len(lineElement[3]) == 1):
                    lineFormat += "  {}     {}"
                elif (len(lineElement[3]) == 3):
                    lineFormat += "  {}   {}"
                else:
                    lineFormat += "  {}  {}"

                if (len(lineElement[4]) == 3):
                    lineFormat += "   {}"
                else:
                    lineFormat += "  {}"

                lineFormat += "  {}"

                if (len(lineElement[6]) == 2):
                    lineFormat += " {}"
                else:
                    lineFormat += "  {}"

                if (len(lineElement) == 38):
                    lineFormat += " {} {} {} {} {} {} {} {}  {}  {}  {}" + \
                                    "  {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}"
                elif (len(lineElement) == 30):
                    lineFormat += " {} {} {} {} {} {} {} {}  {}  {}  {}" + \
                                    "  {} {} {} {} {} {} {} {} {} {} {} {}"

                newLine = lineFormat.format(*lineElement)

                output_file.write(newLine)
                output_file.write("\n")

            # Keep the comments in file
            else:
                newLine = line
                output_file.write(newLine)

    # Move the file back
    shutil.copy2(fileName, pathToFile)

    # Delete the file in the folder
    os.remove(fileName)

def __changeValue(lineElement, newData, atype):

    if (atype == "readOutDim"):
        lineElement[15:19] = newData
    elif (atype == "darkCurrent"):
        lineElement[13:15] = newData
    else:
        print("No type of %s. The type should be 'readOutDim' or 'darkCurrent'." % atype)
        lineElement = None

    return lineElement

def readData(folderPath, fileName, atype):
    """

    Read the data in segmentation or focal plane layout.

    Arguments:
        folderPath {[string]} -- Path to folder.
        fileName {[string]} -- File name ("segmentation.txt", "focalplanelayout.txt").
        atype {[string]} -- Type of data to read ("readOutDim", "darkCurrent", "fieldCenter", "eulerRot").

    Returns:
        [dict] -- Needed CCD data.

    Raises:
        ValueError -- File can not be read.
        ValueError -- Type is not correct.
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

class changePhoSimInstrumentTest(unittest.TestCase):
    
    """ 
    Test the function to update instument segmentation file
    """

    def setUp(self):

        # Get the path of module
        modulePath = getModulePath()

        # Path of data folder
        self.folderPath = os.path.join(modulePath, "test")
        self.fileName = "segmentation.txt"

    def testReadData(self):

        # Modigy the segmentaion: parallel prescan, serial overscan, serial prescan, parallel overscan (pixel)
        newData = ["1", "2", "3", "4"]
        sameCcdList = ["R44_S10", "R44_S01"]
        updateData(self.folderPath, self.fileName, newData, "readOutDim", sameCcdList=sameCcdList)

        # Read the file and print the new segmentation
        ccdData = readData(self.folderPath, self.fileName, "readOutDim")
        self.assertEqual(ccdData["R44_S10_C17"], ["4", "0", "1", "0"])
        self.assertEqual(ccdData["R43_S10_C17"], newData)

        # Update back to original one
        newData = ["4", "0", "1", "0"]
        updateData(self.folderPath, self.fileName, newData, "readOutDim")

if __name__ == '__main__':

    # Do the unit test
    unittest.main() 
