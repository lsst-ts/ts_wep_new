import os, shutil
from astropy.io import fits


class PhoSimImgAdaptor(object):

    def __init__(self, dataFolderPath, pathDestination):

        # Define the path
        self.dataFolderPath = dataFolderPath
        self.pathDestination = pathDestination
        self.pathImage = None
        self.pathBias = None
        self.pathDC = None
        self.pathFlatDome = None

    def config(self, dataFolderPath=None, pathDestination=None, rawDir=None, biasDir=None, darkDir=None, flatDir=None):
        """
        
        Configure the directory.
        
        Keyword Arguments:
            dataFolderPath {str]} -- Path to simulated data. (default: {None})
            pathDestination {[str]} -- Path of destination. (default: {None})
            rawDir {[str]} -- Directory of raw image. (default: {None})
            biasDir {[str]} -- Directory of bias image. (default: {None})
            darkDir {[str]} -- Directory of dark current image. (default: {None})
            flatDir {[str]} -- Directory of flat dome image. (default: {None})
        """

        # Path to input/ output directory
        if (dataFolderPath is not None):
            self.dataFolderPath = dataFolderPath

        if (pathDestination is not None):
            self.pathDestination = pathDestination

        # Path to calibration products folder
        if (rawDir is not None):
            self.pathImage = os.path.join(self.dataFolderPath, rawDir)

        if (biasDir is not None):
            self.pathBias = os.path.join(self.dataFolderPath, biasDir)

        if (darkDir is not None):
            self.pathDC = os.path.join(self.dataFolderPath, darkDir)

        if (flatDir is not None):
            self.pathFlatDome = os.path.join(self.dataFolderPath, flatDir)

    def rearrangeFileForButler(self, aVisit=None, eVisit=None, aFilter=None, atype=None, overwrite=False):
        """
        
        Rearrange the PhoSim images to the folder path needed by data butler to use.
        
        Keyword Arguments:
            aVisit {[int]} -- Visit time of amplifier images. (default: {None})
            eVisit {[int]} -- Visit time of electronic images. (default: {None})
            aFilter {[string]} -- Filter name (u, g, r, i, z, y). (default: {None})
            atype {[string]} -- Type of arrangement: raw, bias, flat, dark. (default: {None})
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})
        
        Returns:
            [string] -- Renamed name of amplifier/ electronic images.
        """

        # Arrange files for Butler to use       
        newAmpImgName = None
        newElecImgName = None

        if (atype =="raw"):
            folderPath = self.pathImage
        elif (atype == "bias"):
            folderPath = self.pathBias
        elif (atype == "flat"):
            folderPath = self.pathFlatDome
        elif (atype == "dark"):
            folderPath = self.pathDC

        # Analyze the file names
        elecImgName, ampImgName = self.__analyzeFileName(folderPath)

        if (aVisit is not None):
            newAmpImgName = self.__arrangeFile(folderPath, self.pathDestination, ampImgName, int(aVisit), 
                                               aFilter, atype, overwrite=overwrite)
        
        if (eVisit is not None):
            if (atype == "raw"):
                newElecImgName = self.__arrangeFile(folderPath, self.pathDestination, elecImgName, int(eVisit), 
                                                    aFilter, atype, overwrite=overwrite)

        return newAmpImgName, newElecImgName

    def __arrangeFile(self, folderPath, pathImage, fileNames, visit, aFilter, atype, overwrite=False):
        """

        Arrange the files to specific folders for following the rule to import
        the image files into data butler.

        Arguments:
            folderPath {[string]} -- Original path to image files.
            pathImage {[string]} -- Destination path to image files.
            fileNames {[string]} -- Names of image files.
            visit {int} -- Visit time to map to LsstSimMapper. This input should be
                           removed in the future.
            aFilter {[string]} -- Filter name (u, g, r, i, z, y).
            atype {[string]} -- Type of arrangement: raw, bias, flat, dark.

        Keyword Arguments:
            overwrite {[boolean]} -- Overwrite the existed files or not. (default: {False})

        Returns:
            [string] -- New image file names to follow the rule of data butler.
        """

        # Arrange the files for data bulter to use

        newFileNames = []
        for file in fileNames:
            name = file.replace(".fits.gz", "")
            numStr = name.split("_")

            # Generate the image holder and new image name
            if (numStr[1] == "e"):
                pathToNewFolder = self.__generateFolder(pathImage, "eimage")

                # Generate the folder name. "v" means the visit. "f" means the filter.
                visitFolderName = "v" + str(visit) + "-f" + aFilter

                # Generate the new image name
                # Consider the WFS 
                if (numStr[4] == "R00" and numStr[5] == "S22") or (numStr[4] == "R04" and numStr[5] == "S20") or (
                    numStr[4] == "R40" and numStr[5] == "S02") or (numStr[4] == "R44" and numStr[5] == "S00"):
                    newName = "_".join(["eimage", str(visit), numStr[4], numStr[5], numStr[6], numStr[7]]) + ".fits.gz"
                else:
                    newName = "_".join(["eimage", str(visit), numStr[4], numStr[5], numStr[6]]) + ".fits.gz"

            elif (numStr[1] == "a"):
                pathToNewFolder = self.__generateFolder(pathImage, atype)

                # Generate the folder name. "v" means the visit. "f" means the filter.
                if atype in ("bias", "dark"):
                    visitFolderName = "v" + str(visit)
                elif atype in ("raw", "flat"):
                    visitFolderName = "v"+ str(visit) +"-f" + aFilter

                # Generate the new image name
                if atype in ("bias", "flat", "dark"):
                    newName = "_".join(["imsim", str(visit), numStr[4], numStr[5], numStr[6]]) + ".fits.gz"
                elif (atype == "raw"):
                    newName = "_".join(["imsim", str(visit), numStr[4], numStr[5], numStr[6], numStr[7]]) + ".fits.gz"

            # Collect the new file names
            newFileNames.append(newName)

            # Check the existence of folder
            pathToNewFolder = self.__generateFolder(pathToNewFolder, visitFolderName)

            # Generate the exposureID (e.g. "E000") folder
            if (atype == "raw"):
                pathToNewFolder = self.__generateFolder(pathToNewFolder, numStr[-1])

            # Generate the chip folder
            pathToNewFolder = self.__generateFolder(pathToNewFolder, numStr[4])

            # Generate the region of chip folder for amplifier images
            if (numStr[1] == "a"):
                pathToNewFolder = self.__generateFolder(pathToNewFolder, numStr[5])

            # Copy the file to the new folder with new name
            newPosition = os.path.join(pathToNewFolder, newName)

            # Check the file exists or not and check the overwrite is true or false
            if (not os.path.exists(newPosition) or overwrite is True):
                shutil.copy2(os.path.join(folderPath, file), newPosition)
            else:
                print("File (%s) exists already." % newPosition)

        return newFileNames

    def collectHeaderInfo(self):
        """
        
        Collect the headers from PhoSim images.
        
        Returns:
            [dictionary] -- Dictionary of header in the key of image name.
        """

        # Collect the header file from PhoSim

        # Define the dictionary of header files
        headerFile = {}

        # Analyze the file names
        elecImgName, ampImgName = self.__analyzeFileName(self.pathImage)
        headerFile.update(self.__getHeader(self.pathImage, elecImgName))
        headerFile.update(self.__getHeader(self.pathImage, ampImgName))

        return headerFile

    def importToButler(self, fileFolder, fileNames, atype, aFilter=None):
        """
        
        Get the exposures as a disctionary of image name and data butler. This function will also
        do the check of header format.
        
        Arguments:
            fileFolder {[string]} -- Path to image folder.
            fileNames {[string]} -- Image name.
            atype {[string]} -- Type of arrangement: eimage, raw, bias, flat, dark.

        Keyword Arguments:
            aFilter {[string]} -- Filter name (u, g, r, i, z, y). (default: {"None"})
        """

        # Check the header file fulfill the rule
        # (This should be removed in the future.)
        for file in fileNames:
            self.__checkButlerFile(fileFolder, file, atype, aFilter=aFilter)

    def __getHeader(self, fileFolder, fileNames):
        """

        Collect the header information of PhoSim output image files. This is for future
        analysis. This function might not be needed in the final.

        Arguments:
            fileFolder {[string]} -- Path to image folder.
            fileNames {[string]} -- Image name.

        Returns:
            [dictionary] -- Dictionary of header with the key contains "e/a" + "raft" + "sensor"
        """

        # Get the information in header of FITS file
        headerFile = {}
        for fileName in fileNames:

            # Generate the file path
            imagePath = os.path.join(fileFolder, fileName)

            # Get tag name
            name = fileName.replace(".fits.gz", "")
            numStr = name.split("_")

            # Get the header files
            # Electronic image
            if (numStr[1] == "e"):
                keyName = "_".join(["e", numStr[4], numStr[5]])
            # Amplifier image
            elif (numStr[1] == "a"):
                keyName = "_".join(["a", numStr[4], numStr[5], numStr[6]])

            # Read the file to get the header file and put into the dictionary
            headerFile[keyName] = fits.getheader(imagePath, 0)

        return headerFile

    def __checkButlerFile(self, fileFolder, fileName, atype, aFilter=None):
        """

        Check and update the header of FITS files if it is needed.

        Arguments:
            fileFolder {[string]} -- Path to file folder.
            fileName {[string]} -- Name of image file.
            atype {[string]} -- Type of arrangement: eimage, raw, bias, flat, dark.

        Keyword Arguments:
            aFilter {[string]} -- Filter name (u, g, r, i, z, y) (default: {"None"}).
        """
        
        # Modify the header of FITS file for butler to import the data

        # Generate the file path
        numStr = fileName.split("_")
        if (atype == "eimage"):
            exposureId = numStr[-1].split(".")[0]
            imagePath = os.path.join(fileFolder, numStr[0], "v"+numStr[1]+"-f"+aFilter,
                                     exposureId, numStr[2], fileName)
        elif (atype == "raw"):
            exposureId = numStr[-1].split(".")[0]
            imagePath = os.path.join(fileFolder, atype, "v"+numStr[1]+"-f"+aFilter,
                                     exposureId, numStr[2], numStr[3], fileName)
        elif atype in ("bias", "dark"):
            imagePath = os.path.join(fileFolder, atype, "v"+numStr[1], numStr[2], 
                                     numStr[3], fileName)
        elif (atype == "flat"):
            imagePath = os.path.join(fileFolder, atype, "v"+numStr[1]+"-f"+aFilter, 
                                     numStr[2], numStr[3], fileName)

        # Get the header
        data, header = fits.getdata(imagePath, header=True)

        # Update the header file
        # Need to modify this in the future once the system is done in PhoSim.
        if (header["RADESYS"]=="J2000" or header["CTYPE1"]!="RA---TAN"):

            # Modify the values
            header["RADESYS"] = "ICRS"
            header["CTYPE1"] = "RA---TAN"

            if (numStr[0] == "imsim"):
                header["VERSION"] = 13357

            if atype in ("bias", "dark", "flat"):
                header.append(card=("EXTTYPE", "IMAGE"))

            # Write into the FITS file
            fits.writeto(imagePath, data, header, overwrite=True)
            print ("Modify the header to import to bulter.")
        else:
            print ("No modification in the header.")

    def __generateFolder(self, pathToFolder, newFolder):
        """

        Generate the folder if the folder does not exist.

        Arguments:
            pathToFolder {[string]} -- Path to folder.
            newFolder {[string]} -- Name of new folder.

        Returns:
            [string] -- Path to generated folder.
        """

        # Generate the path of new folder
        pathToNewFolder = os.path.join(pathToFolder, newFolder)

        # Check the existence of folder
        if (not os.path.exists(pathToNewFolder)):
            # Create the folder if it does not exist
            os.makedirs(pathToNewFolder)

        return pathToNewFolder

    def __analyzeFileName(self, pathImage):
        """

        Analyze the file names of PhoSim output to get the names of electronic and
        amplifier images individually.

        Arguments:
            pathImage {[string]} -- Path to PhoSim output images.

        Returns:
            [string] -- Names of electronic and amplifier images.
        """

        # Get all files
        fileList = os.listdir(pathImage)

        # Analysis file name
        elecImgName = []
        ampImgName = []
        for fileName in fileList:

            # Find ".fits.gz"
            if fileName.endswith(".fits.gz"):

                # Analyze the main file name
                numStr = fileName.split("_")

                # Find the electronic image
                if (numStr[1]=="e"):
                    elecImgName.append(fileName)
                # Find the amplifier image
                elif (numStr[1]=="a"):
                    ampImgName.append(fileName)

        return elecImgName, ampImgName


if __name__ == "__main__":
    pass
