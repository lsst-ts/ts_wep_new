import os

from lsst.ts.wep.SciWFDataCollector import SciWFDataCollector
from lsst.ts.wep.FitsFormat import FitsFormat

if __name__ == "__main__":

    # Instantiate the FitsFormat
    fitsDir = "/home/ttsai/Document/phosimObsData/raw"
    fitsFormat = FitsFormat()
    fitsFormat.config(fitsDir=fitsDir)

    # Raw data path
    daqDataPath = "/home/ttsai/Document/daqData/ncsaSciImg/comcam_a_99999999_f0_R10_S00_C00"
    
    # Get the daq data in numpy array
    data = fitsFormat.getData(daqDataPath)

    # Get the metadata
    dataDict = fitsFormat.getMetaDataFromFileName(daqDataPath)
    dataDict = fitsFormat.addDefaultFakeData(dataDict)

    # Write to the new fits file
    fitsFileName = os.path.basename(daqDataPath)
    fitsFileNameStr = fitsFileName.split("_")
    fitsFileName = "lsst_" + "_".join(fitsFileNameStr[1:]) + "_E000.fits"
    fitsFilePath = fitsFormat.writeNewFits(data, fitsFileName)

    # Update the header
    fitsFormat.config(fitsFilePath=fitsFilePath)
    fitsFormat.updateHeader(dataDict)

    # Gzip the file
    fitsFormat.gzipFits()    

    # Ingest the fits image
    pathOfRawData = "/home/ttsai/Document/phosimObsData/raw"
    destinationPath = "/home/ttsai/Document/phosimObsData/input"
    sciWFData = SciWFDataCollector()
    sciWFData.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

    fitsFileArg = "lsst_a_99999999*.fits.gz"
    sciWFData.ingestSimImages(fitsFileArg=fitsFileArg)
