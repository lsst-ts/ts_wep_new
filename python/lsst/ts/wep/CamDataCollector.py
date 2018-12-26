import os
from lsst.ts.wep.Utility import runProgram, writeFile


class CamDataCollector(object):

    def __init__(self, destDir):

        self.destDir = destDir

    def genPhoSimMapper(self):

        fileName = "_mapper"
        filePath = os.path.join(self.destDir, fileName)

        content = "lsst.obs.lsst.phosim.PhosimMapper"

        writeFile(filePath, content)

    def ingestCalibs(self, calibFiles):
        
        command = "ingestCalibs.py"
        argstring = "%s %s --validity 99999 --output %s" % (
                        self.destDir, calibFiles, self.destDir)
        runProgram(command, argstring=argstring)

    def ingestImages(self, imgFiles):
        
        command = "ingestImages.py"
        argstring = "%s %s" % (self.destDir, imgFiles)
        runProgram(command, argstring=argstring)


if __name__ == "__main__":
    pass
