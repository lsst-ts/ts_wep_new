import os
import shutil
import warnings

from lsst.ts.wep.Utility import runProgram, writeFile


class CamIsrWrapper(object):

    def __init__(self, destDir):
        """Initialize the camera ISR wrapper class.

        ISR: Instrument signature removal.

        Parameters
        ----------
        destDir : str
            Destination directory.
        """

        self.destDir = destDir

        self.doBias = False
        self.doDark = False
        self.doFlat = False
        self.doFringe = False
        self.doDefect = False

        self.isrConfigFilePath = None

    def config(self, doBias=False, doDark=False, doFlat=False,
               doFringe=False, doDefect=False, fileName="isr_config.py"):
        """Do the ISR configuration.

        ISR: Instrument signature removal.

        Parameters
        ----------
        doBias : bool, optional
            Do the bias correction. (the default is False.)
        doDark : bool, optional
            Do the dark correction. (the default is False.)
        doFlat : bool, optional
            Do the flat correction. (the default is False.)
        doFringe : bool, optional
            Do the fringe correction. (the default is False.)
        doDefect : bool, optional
            Do the defect correction. (the default is False.)
        fileName : str, optional
            Config override file name. (the default is "isr_config.py".)
        """

        self.doBias = doBias
        self.doDark = doDark
        self.doFlat = doFlat
        self.doFringe = doFringe
        self.doDefect = doDefect

        self._setIsrConfigfile(fileName)

    def _setIsrConfigfile(self, fileName):
        """Set the ISR configuration file.

        ISR: Instrument signature removal.

        Parameters
        ----------
        fileName : str
            ISR configuration file name.
        """

        filePath = os.path.join(self.destDir, fileName)

        content = "config.isr.doBias=%s\n" % self.doBias
        content += "config.isr.doDark=%s\n" % self.doDark
        content += "config.isr.doFlat=%s\n" % self.doFlat
        content += "config.isr.doFringe=%s\n" % self.doFringe
        content += "config.isr.doDefect=%s\n" % self.doDefect

        try:
            writeFile(filePath, content)
            self.isrConfigFilePath = filePath
        except Exception:
            raise

    def doISR(self, inputDir, rerunName="run1"):
        """Do the ISR.

        ISR: Instrument signature removal.

        Parameters
        ----------
        inputDir : str
            Input data directory.
        rerunName : str, optional
            Rerun name. (the default is "run1".)
        """

        # Work around the bug of obs_lsst in w_2019_20 that can not do the ISR
        # continuously. This part should be removed in the final.
        rerunDirPath = os.path.join(inputDir, "rerun")
        if os.path.exists(rerunDirPath):
            shutil.rmtree(rerunDirPath)
            warnings.warn("Rerun dir exists. Remove it to work around.",
                          category=UserWarning)

        # Do the ISR
        command = "runIsr.py"

        argstring = "%s --id --rerun=%s" % (inputDir, rerunName)
        if (self.isrConfigFilePath is not None):
            argstring += " --configfile %s" % self.isrConfigFilePath

        runProgram(command, argstring=argstring)


if __name__ == "__main__":
    pass
