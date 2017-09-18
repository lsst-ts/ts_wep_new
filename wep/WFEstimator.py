import os, sys

from cwfs.Instrument import Instrument
from cwfs.Algorithm import Algorithm
from cwfs.Image import Image
from cwfs.Tool import plotImage

class WFEstimator(object):

	# Constant
	INTRA = "intra"
	EXTRA = "extra"

	def __init__(self, instruFolder, algoFolderPath):

		self.algo = Algorithm(algoFolderPath)
		self.inst = Instrument(instruFolder)
		self.ImgIntra = Image()
		self.ImgExtra = Image()
		self.opticalModel = None

	def config(self, solver="exp", instName="lsst", opticalModel="offAxis", debugLevel=0):
		"""
		
		Configure the TIE solver.
		
		Keyword Arguments:
			solver {str} -- Algorithm to solve the Poisson's equation in the transport of 
							intensity equation (TIE). It can be "fft" or "exp" here. 
							(default: {"exp"})
			instName {str} -- Instrument name. It is "lsst" in the baseline. (default: {"lsst"})
			opticalModel {str} -- Optical model. It can be "paraxial", "onAxis", or "offAxis". 
								  (default: {"offAxis"})
			debugLevel {number} -- Show the information under the running. If the value is higher, 
								   the information shows more. It can be 0, 1, 2, or 3. 
								   (default: {0})
		
		Raises:
			ValueError -- Wrong instrument name.
			ValueError -- No intra-focal image.
			ValueError -- Wrong Poisson solver name.
			ValueError -- Wrong optical model.
		"""

		# Check the inputs and assign the parameters used in the TIE
		# Need to change the way to hole the classes of Instrument and Algorithm

		if (instName != "lsst"):
			raise ValueError("Instrument can not be '%s'." % instName)
		else:
			if (not self.ImgIntra.sizeinPix):
				raise ValueError("The intra-focal image does not be set yet.")
			else:
				self.inst.config(instName, self.ImgIntra.sizeinPix)

		if solver not in ("exp", "fft"):
			raise ValueError("Poisson solver can not be '%s'." % solver)
		else:
			self.algo.config(solver, self.inst, debugLevel=debugLevel)
		
		if opticalModel not in ("paraxial", "onAxis", "offAxis"):
			raise ValueError("Optical model can not be '%s'." % opticalModel)
		else:
			self.opticalModel = opticalModel

	def setImg(self, fieldXY, image=None, imageFile=None, defocalType=None):
		"""
		
		Set the wavefront image.
		
		Arguments:
			fieldXY {[float]} -- Position of donut on the focal plane in degree for intra- and extra-focal
                             	 images.
		
		Keyword Arguments:
			image {[float]} -- Array of image. (default: {None})
			imageFile {[str]} -- Path of image file. (default: {None})
			defocalType {[str]} -- Type of image. It should be "intra" or "extra". (default: {None})
		
		Raises:
			ValueError -- Wrong defocal type.
		"""
		
		# Check the defocal type
		if defocalType not in (self.INTRA, self.EXTRA):
			raise ValueError("Defocal type can not be '%s'." % defocalType)

		# Read the image and assign the type
		if (defocalType == self.INTRA):
			self.ImgIntra.setImg(fieldXY, imageFile=imageFile, atype=defocalType)
		elif (defocalType == self.EXTRA):
			self.ImgExtra.setImg(fieldXY, imageFile=imageFile, atype=defocalType)

	def calWfsErr(self, tol=1e-3, showZer=False, showPlot=False):
		"""
		
		Calculate the wavefront error.
		
		Keyword Arguments:
			tol {number} -- Tolerance of difference of coefficients of Zk polynomials compared with
                            the previours iteration. (default: {1e-3})
			showZer {bool} -- Decide to show the annular Zernike polynomails or not. (default: {False})
			showPlot {bool} -- Decide to show the plot or not. (default: {False})
		
		Returns:
			[float] -- Coefficients of Zernike polynomials (z4 - z22).
		"""

		# Calculate the wavefront error.
		# Run cwfs
		self.algo.runIt(self.inst, self.ImgIntra, self.ImgExtra, self.opticalModel, tol=tol)

		# Show the Zernikes Zn (n>=4)
		if (showZer):
			self.algo.outZer4Up(showPlot=showPlot)

		return self.algo.zer4UpNm

	def outParam(self, filename=None):
	    """
	    
	    Put the information of images, instrument, and algorithm on terminal or file.
	    	    
	    Keyword Arguments:
	        filename {[str]} -- Name of output file. (default: {None})
	    """

	    # Write the parameters into a file if needed.
	    if (filename is not None):
	        fout = open(filename, "w")
	    else:
	        fout = sys.stdout

	    # Write the information of image and optical model
	    if self.ImgIntra is not None:
		    fout.write("intra image: \t %s \t field in deg =(%6.3f, %6.3f)\n" %
		               (self.ImgIntra.name, self.ImgIntra.fieldX, self.ImgIntra.fieldY))

	    if self.ImgExtra is not None:
		    fout.write("extra image: \t %s \t field in deg =(%6.3f, %6.3f)\n" %
		               (self.ImgExtra.name, self.ImgExtra.fieldX, self.ImgExtra.fieldY))

	    if self.opticalModel is not None:
		    fout.write("Using optical model:\t %s\n" % self.opticalModel)
	    
	    # Read the instrument file
	    if self.algo is not None:
		    self.__readConfigFile(fout, self.inst, "instrument")

	    # Read the algorithm file
	    if self.algo is not None:
		    self.__readConfigFile(fout, self.algo, "algorithm")

	    # Close the file
	    if (filename is not None):
	        fout.close()

	def __readConfigFile(self, fout, config, configName):
	    """
	    
	    Read the configuration file
	    
	    Arguments:
	        fout {[file]} -- File instance.
	        config {[metadata]} -- Instance of configuration. It is Instrument or Algorithm here.
	        configName {[str]} -- Name of configuration.
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

    # Define the instrument folder
    instruFolder = "/Users/Wolf/Documents/stash/ts_lsst_wep_27/instruData"

    # Define the algorithm folder 
    algoFolderPath = "/Users/Wolf/Documents/stash/ts_lsst_wep_27/algo"

    # Define the image folder and image names
    # Image data -- Don't know the final image format. 
    # It is noted that image.readFile inuts is based on the txt file  
    imageFolderPath = "/Users/Wolf/Documents/stash/ts_lsst_wep_27/tests/testImages/LSST_NE_SN25"
    intra_image_name = "z11_0.25_intra.txt"
    extra_image_name = "z11_0.25_extra.txt"

    # Path to image files
    intraImgFile = os.path.join(imageFolderPath, intra_image_name)
    extraImgFile = os.path.join(imageFolderPath, extra_image_name)

    # Field XY position
    fieldXY = [1.185, 1.185]

    # Declare the WFEsitmator
    wfsEst = WFEstimator(instruFolder, algoFolderPath)

    # Setup the images
    wfsEst.setImg(fieldXY, imageFile=intraImgFile, defocalType="intra")
    wfsEst.setImg(fieldXY, imageFile=extraImgFile, defocalType="extra")

    # Setup the configuration
    # If the configuration is reset, the images are needed to be set again. This bug should
    # be solved.
    wfsEst.config(solver="exp", debugLevel=0)

    # # Evaluate the wavefront error
    # plotImage(wfsEst.ImgIntra.image, title="intra image")
    # plotImage(wfsEst.ImgExtra.image, title="extra image")

    # # Evalute the wavefront error
    zer4UpNm = wfsEst.calWfsErr(showZer=True)

    # # Show the compensated images
    # plotImage(wfsEst.ImgIntra.image, title="Compensated intra image")
    # plotImage(wfsEst.ImgExtra.image, title="Compensated extra image")

    # Show the parameters
    # wfsEst.outParam()






