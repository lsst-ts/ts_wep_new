import os
from lsst.ts.wep.Utility import runProgram, writeFile, \
    getObsLsstCmdTaskConfigDir


class CamDataCollector(object):

    def __init__(self, destDir):
        """Initialize the camera data collector class.

        Parameters
        ----------
        destDir : str
            Destination directory.
        """

        self.destDir = destDir

    def genPhoSimMapper(self):
        """Generate the PhoSim mapper."""

        fileName = "_mapper"
        filePath = os.path.join(self.destDir, fileName)

        content = "lsst.obs.lsst.phosim.PhosimMapper"

        writeFile(filePath, content)

    def ingestCalibs(self, calibFiles):
        """Ingest the calibration files.

        Parameters
        ----------
        calibFiles : str
            Calibration files.
        """

        command = "ingestCalibs.py"
        argstring = "%s %s --validity 99999 --output %s" % (
            self.destDir, calibFiles, self.destDir)
        runProgram(command, argstring=argstring)

    def ingestImages(self, imgFiles):
        """Ingest the amplifier image files.

        Parameters
        ----------
        imgFiles : str
            Image files.
        """

        self._ingestImagesByButler(imgFiles)

    def _ingestImagesByButler(self, imgFiles, configFile=None):
        """Ingest the images by butler.

        Parameters
        ----------
        imgFiles : str
            Image files.
        configFile : str, optional
            Config override file(s). (the default is None.)
        """

        command = "ingestImages.py"

        argstring = "%s %s" % (self.destDir, imgFiles)
        if (configFile is not None):
            argstring += " --configfile %s" % configFile

        runProgram(command, argstring=argstring)

    def ingestEimages(self, imgFiles):
        """Ingest the PhoSim eimage files.

        Parameters
        ----------
        imgFiles : str
            Image files.
        """

        eimgConfigFile = os.path.join(getObsLsstCmdTaskConfigDir(), "phosim",
                                      "ingestEimg.py")
        self._ingestImagesByButler(imgFiles, configFile=eimgConfigFile)


if __name__ == "__main__":
    pass
