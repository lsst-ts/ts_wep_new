import sys

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.cwfs.CompensationImageDecorator import CompensationImageDecorator


class WfEstimator(object):

    def __init__(self, instruFolderPath, algoFolderPath):
        """Initialize the wavefront estimator class.

        Parameters
        ----------
        instruFolderPath : str
            Path to instrument directory.
        algoFolderPath : str
            Path to algorithm directory.
        """

        self.algo = Algorithm(algoFolderPath)
        self.inst = Instrument(instruFolderPath)
        self.ImgIntra = CompensationImageDecorator()
        self.ImgExtra = CompensationImageDecorator()
        self.opticalModel = ""

        self.sizeInPix = 0

    def getAlgo(self):
        """Get the algorithm object.

        Returns
        -------
        Algorithm
            Algorithm object.
        """

        return self.algo

    def getInst(self):
        """Get the instrument object.

        Returns
        -------
        Instrument
            Instrument object.
        """

        return self.inst

    def getIntraImg(self):
        """Get the intra-focal donut image.

        Returns
        -------
        CompensationImageDecorator
            Intra-focal donut image.
        """

        return self.ImgIntra

    def getExtraImg(self):
        """Get the extra-focal donut image.

        Returns
        -------
        CompensationImageDecorator
            Extra-focal donut image.
        """

        return self.ImgExtra

    def getOptModel(self):
        """Get the optical model.

        Returns
        -------
        str
            Optical model.
        """

        return self.opticalModel

    def getSizeInPix(self):
        """Get the donut image size in pixel defined by the config() function.

        Returns
        -------
        int
            Donut image size in pixel
        """

        return self.sizeInPix

    def reset(self):
        """

        Reset the calculation for the new input images with the same algorithm
        settings.
        """

        self.algo.reset()

    def config(self, solver="exp", instName="lsst", opticalModel="offAxis",
               defocalDisInMm=None, sizeInPix=120, debugLevel=0):
        """Configure the TIE solver.

        Parameters
        ----------
        solver : str, optional
            Algorithm to solve the Poisson's equation in the transport of
            intensity equation (TIE). It can be "fft" or "exp" here. (the
            default is "exp".)
        instName : str, optional
            Instrument name. It is "lsst" in the baseline. (the default is
            "lsst".)
        opticalModel : str, optional
            Optical model. It can be "paraxial", "onAxis", or "offAxis". (the
            default is "offAxis".)
        defocalDisInMm : float, optional
            Defocal distance in mm. (the default is None.)
        sizeInPix : int, optional
            Wavefront image pixel size. (the default is 120.)
        debugLevel : int, optional
            Show the information under the running. If the value is higher,
            the information shows more. It can be 0, 1, 2, or 3. (the default
            is 0.)

        Raises
        ------
        ValueError
            Wrong instrument name.
        ValueError
            No intra-focal image.
        ValueError
            Wrong Poisson solver name.
        ValueError
            Wrong optical model.
        """

        # Check the inputs and assign the parameters used in the TIE
        # Need to change the way to hold the instance of Instrument and
        # Algorithm

        # Update the isnstrument name
        if (defocalDisInMm is not None):
            instName = instName + str(int(10*defocalDisInMm))

        if instName not in ("lsst", "lsst05", "lsst10", "lsst15", "lsst20",
                            "lsst25", "comcam10", "comcam15", "comcam20"):
            raise ValueError("Instrument can not be '%s'." % instName)

        # Set the available wavefront image size (n x n)
        self.sizeInPix = int(sizeInPix)

        # Configurate the instrument
        self.inst.config(instName, self.sizeInPix)

        if solver not in ("exp", "fft"):
            raise ValueError("Poisson solver can not be '%s'." % solver)
        else:
            self.algo.config(solver, self.inst, debugLevel=debugLevel)

        if opticalModel not in ("paraxial", "onAxis", "offAxis"):
            raise ValueError("Optical model can not be '%s'." % opticalModel)
        else:
            self.opticalModel = opticalModel

    def setImg(self, fieldXY, image=None, imageFile=None, defocalType=None):
        """Set the wavefront image.

        Parameters
        ----------
        fieldXY : tuple or list
            Position of donut on the focal plane in degree for intra- and
            extra-focal images.
        image : numpy.ndarray, optional
            Array of image. (the default is None.)
        imageFile : str, optional
            Path of image file. (the default is None.)
        defocalType : str, optional
            Type of image. It should be "intra" or "extra". (the default is
            None.)

        Raises
        ------
        ValueError
            Wrong defocal type.
        """

        # Check the defocal type
        if defocalType not in (self.ImgIntra.INTRA, self.ImgIntra.EXTRA):
            raise ValueError("Defocal type can not be '%s'." % defocalType)

        # Read the image and assign the type
        if (defocalType == self.ImgIntra.INTRA):
            self.ImgIntra.setImg(fieldXY, image=image, imageFile=imageFile,
                                 atype=defocalType)
        elif (defocalType == self.ImgIntra.EXTRA):
            self.ImgExtra.setImg(fieldXY, image=image, imageFile=imageFile,
                                 atype=defocalType)

    def calWfsErr(self, tol=1e-3, showZer=False, showPlot=False):
        """Calculate the wavefront error.

        Parameters
        ----------
        tol : float, optional
            [description] (the default is 1e-3.)
        showZer : bool, optional
            Decide to show the annular Zernike polynomails or not. (the default
            is False.)
        showPlot : bool, optional
            Decide to show the plot or not. (the default is False.)

        Returns
        -------
        numpy.ndarray
            Coefficients of Zernike polynomials (z4 - z22).

        Raises
        ------
        RuntimeError
            Input image shape is wrong.
        """

        # Check the image size
        for img in (self.ImgIntra, self.ImgExtra):
            d1, d2 = img.image.shape
            if (d1 != self.sizeInPix) or (d2 != self.sizeInPix):
                raise RuntimeError("Input image shape is (%d, %d), not required (%d, %d)" % (
                    d1, d2, self.sizeInPix, self.sizeInPix))

        # Calculate the wavefront error.
        # Run cwfs
        self.algo.runIt(self.inst, self.ImgIntra, self.ImgExtra,
                        self.opticalModel, tol=tol)

        # Show the Zernikes Zn (n>=4)
        if (showZer):
            self.algo.outZer4Up(showPlot=showPlot)

        return self.algo.zer4UpNm

    def outParam(self, filename=None):
        """Put the information of images, instrument, and algorithm on terminal
        or file.

        Parameters
        ----------
        filename : str, optional
            Name of output file. (the default is None.)
        """

        # Write the parameters into a file if needed.
        if (filename is not None):
            fout = open(filename, "w")
        else:
            fout = sys.stdout

        # Write the information of image and optical model
        if (self.ImgIntra.name is not None):
            fout.write("Intra image: \t %s\n" % self.ImgIntra.name)

        if (self.ImgIntra.fieldX is not None):
            fout.write("Intra image field in deg =(%6.3f, %6.3f)\n" % (
                self.ImgIntra.fieldX, self.ImgIntra.fieldY))

        if (self.ImgExtra.name is not None):
            fout.write("Extra image: \t %s\n" % self.ImgExtra.name)

        if (self.ImgExtra.fieldX is not None):
            fout.write("Extra image field in deg =(%6.3f, %6.3f)\n" % (
                self.ImgExtra.fieldX, self.ImgExtra.fieldY))

        if (self.opticalModel is not None):
            fout.write("Using optical model:\t %s\n" % self.opticalModel)

        # Read the instrument file
        if (self.inst.filename is not None):
            self._readConfigFile(fout, self.inst, "instrument")

        # Read the algorithm file
        if (self.algo.filename is not None):
            self._readConfigFile(fout, self.algo, "algorithm")

        # Close the file
        if (filename is not None):
            fout.close()

    def _readConfigFile(self, fout, config, configName):
        """Read the configuration file

        Parameters
        ----------
        fout : file
            File instance.
        config : Instrument or Algorithm
            Instance of configuration. It is Instrument or Algorithm here.
        configName : str
            Name of configuration.
        """

        # Create a new line
        fout.write("\n")

        # Open the file
        fconfig = open(config.filename)
        fout.write("---" + configName + " file: --- %s ----------\n" % config.filename)

        # Read the file information
        iscomment = False
        for line in fconfig:
            line = line.strip()
            if (line.startswith("###")):
                iscomment = ~iscomment
            if (not(line.startswith("#")) and (not iscomment) and len(line) > 0):
                fout.write(line + "\n")

        # Close the file
        fconfig.close()


if __name__ == "__main__":
    pass
