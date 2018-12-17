import os

from lsst.ts.wep.SciWFDataCollector import SciWFDataCollector

if __name__ == "__main__":

    # This script is to demonstrate the ingestion of PhoSim images, make the fake 
    # flat calibratioin products, and import them to the data butler.
    
    # Directory
    homeDir = os.path.expanduser("~")
    pathOfRawData = os.path.join(homeDir, "Document", "phosimObsData", "raw")
    calibDestDir = os.path.join(homeDir, "Document", "phosimObsData", "calibs")
    destinationPath = os.path.join(homeDir, "Document", "phosimObsData", "input")

    # Instantiate
    sciWfDataCollector = SciWFDataCollector()
    sciWfDataCollector.config(pathOfRawData=pathOfRawData, destinationPath=destinationPath)

    # Get the mapper
    sciWfDataCollector.setMapper()

    # Import the raw data
    sciWfDataCollector.ingestSimImages()

    # Get the calibration products
    sciWfDataCollector.getCalibsData(calibDestDir)

    # Ingest the calibration products
    sciWfDataCollector.ingestCalibs(calibDestDir)
