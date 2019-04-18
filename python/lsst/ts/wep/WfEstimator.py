import sys

from lsst.ts.wep.cwfs.Instrument import Instrument
from lsst.ts.wep.cwfs.Algorithm import Algorithm
from lsst.ts.wep.cwfs.CompensationImageDecorator import CompensationImageDecorator
from lsst.ts.wep.Utility import DefocalType


class WfEstimator(object):

    def __init__(self, instDir, algoDir):
        """Initialize the wavefront estimator class.

        Parameters
        ----------
        instDir : str
            Path to instrument directory.
        algoDir : str
            Path to algorithm directory.
        """

        self.inst = Instrument(instDir)
        self.algo = Algorithm(algoDir)

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

    def config(self, solver="exp", instName="lsst15", opticalModel="offAxis",
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
            "lsst15".)
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

        if instName not in ("lsst05", "lsst10", "lsst15", "lsst20", "lsst25",
                            "comcam10", "comcam15", "comcam20"):
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

    def setImg(self, fieldXY, defocalType, image=None, imageFile=None):
        """Set the wavefront image.

        Parameters
        ----------
        fieldXY : tuple or list
            Position of donut on the focal plane in degree for intra- and
            extra-focal images.
        defocalType : enum 'DefocalType'
            Defocal type of image.
        image : numpy.ndarray, optional
            Array of image. (the default is None.)
        imageFile : str, optional
            Path of image file. (the default is None.)

        Raises
        ------
        ValueError
            Wrong defocal type.
        """

        if (defocalType == DefocalType.Intra):
            img = self.ImgIntra
        elif (defocalType == DefocalType.Extra):
            img = self.ImgExtra

        img.setImg(fieldXY, defocalType, image=image, imageFile=imageFile)

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
            d1, d2 = img.getImg().shape
            if (d1 != self.sizeInPix) or (d2 != self.sizeInPix):
                raise RuntimeError("Input image shape is (%d, %d), not required (%d, %d)" % (
                    d1, d2, self.sizeInPix, self.sizeInPix))

        # Calculate the wavefront error.
        # Run cwfs
        self.algo.runIt(self.ImgIntra, self.ImgExtra, self.opticalModel, tol=tol)

        # Show the Zernikes Zn (n>=4)
        if (showZer):
            self.algo.outZer4Up(showPlot=showPlot)

        return self.algo.getZer4UpInNm()


if __name__ == "__main__":
    pass
