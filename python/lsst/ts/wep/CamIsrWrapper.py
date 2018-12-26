import os

from lsst.ts.wep.Utility import runProgram, writeFile


class CamIsrWrapper(object):

    def __init__(self, destDir):
        
        self.destDir = destDir

        self.doBias = False
        self.doDark = False
        self.doFlat = False
        self.doFringe = False
        self.doDefect = False

        self.isrConfigFilePath = None

    def config(self, doBias=False, doDark=False, doFlat=False,
               doFringe=False, doDefect=False, fileName="isr_config.py"):

        self.doBias = doBias
        self.doDark = doDark
        self.doFlat = doFlat
        self.doFringe = doFringe
        self.doDefect = doDefect

        self._setIsrConfigfile(fileName)

    def _setIsrConfigfile(self, fileName):

        filePath = os.path.join(self.destDir, fileName)

        content = "config.isr.doBias=%s\n" % self.doBias
        content += "config.isr.doDark=%s\n" % self.doDark
        content += "config.isr.doFlat=%s\n" % self.doFlat
        content += "config.isr.doFringe=%s\n" % self.doFringe
        content += "config.isr.doDefect=%s\n" % self.doDefect

        try:
            writeFile(filePath, content)
            self.isrConfigFilePath = filePath
        except Exception as e:
            raise

    def doISR(self, inputDir, rerunName="run1"):
        
        command = "runIsr.py"

        argstring = "%s --id --rerun=%s" % (inputDir, rerunName)
        if (self.isrConfigFilePath is not None):
            argstring += " --configfile %s" % self.isrConfigFilePath

        runProgram(command, argstring=argstring)


if __name__ == "__main__":
    pass
