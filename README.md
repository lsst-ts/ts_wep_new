# Wavefront Estimation Pipeline (WEP)

*This module calculates the wavefront error in annular Zernike polynomials up to 22 terms based on the intra- and extra-focal donut images in the large synoptic survey telescope (LSST). The image quality is decided by z4-z22. The wavefront error is determined by solving the transport of intensity equation (TIE) that approximates the change of intensity mainly comes from the wavefront error.*

## 1. Version History

The version history is [here](./doc/VersionHistory.md).

*Author: Te-Wei Tsai*
<br/>
*Date: 2-15-2019*

## 2. Platform

- *CentOS 7*
- *python: 3.6.6*
- *scientific pipeline (newinstall.sh from master branch)*

## 3. Needed Package

- *lsst_sims (tag: sims_w_2019_02)*
- *lsst_distrib (tag: w_2019_02)*
- *obs_lsst - master branch (commit: 69b4a98)*
- *phosim_utils - master branch (commit: b8d87d9)*
- *scikit-image*

## 4. Compile cwfs
*To compile the code, execute the following command at the directory of WEP:*
```
python builder/setup.py build_ext --build-lib python/lsst/ts/wep/cwfs/lib
```

## 5. Install the LSST Packages, obs_lsst, and phosim_utils

*1. Setup the LSST environment by `source $LSST_DIR/loadLSST.bash`. LSST_DIR is the directory of scientific pipeline.*
<br/>
*2. Install the lsst_sims by `eups distrib install lsst_sims -t sims_w_2019_02`.*
<br/>
*3. Install the lsst_distrib by `eups distrib install lsst_distrib -t w_2019_02`.*
<br/>
*4. Fix the path by `curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python`. The [shebangtron repo](https://github.com/lsst/shebangtron) has the further discussion of this.*
<br/>
*5. Clone the repository of [obs_lsst](https://github.com/lsst/obs_lsst) to some other directory. Under the obs_lsst directory, use `setup -k -r .` to setup the package in eups and use `scons` to build the module. It is noted that the build process is only needed for the first time.*
<br/>
*6. Do the step 5 for the repository of [phosim_utils](https://github.com/lsst-dm/phosim_utils.git).*

## 6. Pull the Built Image from Docker Hub

*Pull the built docker image by `docker pull lsstts/aos:w_2019_02`. The scientific pipeline and lsst packages are installed already. For the details of docker image, please follow the [docker aos image](https://hub.docker.com/r/lsstts/aos).*

## 7. DM Command Line Task (obs_lsst and phosim_utils)

*1. Make the faked flat images. Flats only need to be made once. They can then be shared between repos. The flats can be faked with (1) all sensors, (2) corner wavefront sensors, or (3) user-specified sensor list.*
```
cd $work_dir
mkdir fake_flats
cd fake_flats/
makeGainImages.py
cd ..
```
*2. Repackage the PhoSim output amplifiers. The data needs to be put in single 16 extension MEFs (Multi-Extension FITS) for processing.*
```
phosim_repackager.py $phosim_amp_dir --out_dir=repackaged_files
```
*3. Make the repository for butler to use, ingest the images, and ingest the calibration products.*
```
mkdir input
echo lsst.obs.lsst.phosim.PhosimMapper > input/_mapper
ingestImages.py input repackaged_files/*.fits
ingestCalibs.py input fake_flats/* --validity 99999 --output input`
```
*4. Make the config override file to turn only flat field on.*
```
echo "config.isr.doBias=False
config.isr.doDark=False
config.isr.doFlat=True
config.isr.doFringe=False
config.isr.doDefect=False" >isr_config.py
```
*5. Run the instrument signature removal (ISR).*
```
runIsr.py input --id --rerun=run1 --configfile isr_config.py
```

## 8. Use of Module

*1. Setup the DM environment.*
```
source $path_of_lsst_scientific_pipeline/loadLSST.bash
setup sims_catUtils -t sims_w_2019_02
```
*2. Setup the WEP environment.*
```
export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep/python
```

## 9. Integrate with SAL (Need to be refined with the SAL module)

*Some environment paths defined in ts_sal/setup.env need to be modified to use lsst stack with SAL.*

*Need to setup the following path variables: LSST_SDK_INSTALL, OSPL_HOME, PYTHON_BUILD_VERSION, and PYTHON_BUILD_LOCATION.*

*1. `PYTHON_BUILD_LOCATION=$lsst_stack_python_directory`. e.g. `PYTHON_BUILD_LOCATION=/home/ttsai/Document/lsst16/python/miniconda3-4.5.4`.*

*2. In ts_sal/setup.env, use `LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SAL_HOME}/lib` instead of `LD_LIBRARY_PATH=${SAL_HOME}/lib`.*

## 10. SAL XML Model (Need to be refined with the SAL module)

*The SAL xml model are in the 'sal_interfaces' of ts_xml repository. The branch used is the develop branch. The CSC keywords are 'tcsOfc' and 'tcsWEP' for active optics to use. The xml files can be found in the related directory. The way to generate the SAL py libraries can follow the ts_sal manual.*

## 11. Content

*This module contains the following classes ([class diagram](./doc/wepClassDiag.png)):*

- **ButlerWrapper**: Wrapper of DM butler class to get the raw and post-ISR CCD image.
- **CamDataCollector**: Ingest the amplifier images and calibration products based on the DM command line task.
- **CamIsrWrapper**: Do the ISR and assemble the CCD images based on the DM command line task.
- **SourceSelector**: Query the bright star catalog (BSC) to select the available target to calculate the wavefront error.
- **SourceProcessor**: Process the post-ISR images to get the clean star images with measured optical field position. The deblending algorithm is used to get the single target bright star image if the neighboring star exists.
- **WfEstimator**: Calculate the wavefront error in annular Zernike polynomials up to 22 terms based on the defocal donut images.
- **DefocalImage**: Defocal image class that provides the accessor methods.
- **DonutImage**: Donut image class that provides the accessor methods.
- **WepController**: High level class to use the WEP package.
- **Utility**: Utility functions used in WEP.
- **PlotUtil**: Plot utility functions used in WEP.

*There are the following modules in WEP:*

- **Bright Star Catalog (bsc)**: This module queries the bright star catalog and gets the scientific target. The class diagram is [here](./doc/bscClassDiag.png).
    - **CamFactory**: Camera factory to create the concrete camera object.
    - **CameraData**: Camera data class as the parent of specific camera child class.
    - **ComCam**: Commissioning camera class. The parent class is the CameraData class.
    - **LsstCam**: LSST camera class to use the corner wavefront sensor. The parent class is the CameraData class.
    - **LsstFamCam**: Lsst camera class to use the full-array mode (FAM). The wavefront sensor is the scientific sensor. The parent class is the CameraData class.
    - **DatabaseFactory**: Database factory to create the concrete database object.
    - **DefaultDatabase**: Default database class as the parent of specific database child class.
    - **LocalDatabase**: Local database class. The parent class is the DefaultDatabase class.
    - **LocalDatabaseForStarFile**: Local database class to read the star file. The parent class is the LocalDatabase class.
    - **StarData**: Star data class for the scientific target star.
    - **NbrStar**: Neighboring star class to have the bright star and the related neighboring stars.
    - **Filter**: Filter class to provide the scientific target star magnitude boundary.
    - **WcsSol**: Wavefront coordinate system (WCS) solution class to map the sky position to camera pixel position and vice versa.
    - **PlotStarFunc**: Plot funtions used in this bsc module.

- **Stars Deblending (deblend)**: This module does the image deblending. The class diagram is [here](./doc/deblendClassDiag.png).
    - **AdapThresImage**: Adapted threshold image class to get the donut centor according to the binary image by the adapted threshold method.
    - **BlendedImageDecorator**: Blended image decorator class to do the donut deblending.
    - **nelderMeadModify**: Do the numerical optimation according to the Nelder-Mead algorithm.

- **Curvature Wavefront Sensor (cwfs)**: This module calculates the wavefront error by solving the TIE. The class diagram is [here](./doc/cwfsClassDiag.png).
    - **Algorithm**: Algorithm class to solve the TIE to get the wavefront error.
    - **CompensationImageDecorator**: Compensation image decorator class to project the donut image from the image plane to the pupil plane.
    - **Image**: Image class to have the function to get the donut center.
    - **Instrument**: Instrument class to have the instrument information used in the Algorithm class to solve the TIE.
    - **Tool**: Annular Zernike polynomials related functions.

- **Control Interface (ctrlIntf)**: This module provides the interface classes to the main telescope active optics system (MTAOS). The factory pattern is applied to support the multiple instruments. The class diagram is [here](./doc/ctrlIntfClassDiag.png).
    - **WEPCalculationFactory**: Factory for creating the correct WEP calculation based off the camera type currently being used.
    - **WEPCalculation**: Base class for converting the wavefront images into wavefront errors.
    - **WEPCalculationOfPiston**: The child class of WEPCalculation that gets the defocal images by the camera piston.
    - **WEPCalculationOfLsstCam**: The concrete child class of WEPCalculation of the LSST camera (corner wavefront sensor).
    - **WEPCalculationOfComCam**: The concrete child class of WEPCalculationOfPiston of the commionning camera (ComCam).
    - **WEPCalculationOfLsstFamCam**: The concrete child class of WEPCalculationOfPiston of the LSST full-array mode (FAM) camera.
    - **SensorWavefrontData**: Sensor wavefront data class that has the information of sensor Id, list of donut, master donut, and wavefront error.
    - **WcsData**: Contains the world coordinate system (WCS) data of a camera.
    - **AstWcsSol**: AST world coordinate system (WCS) solution provided by DM team.
    - **RawExpData**: Raw exposure data class to populate the information of visit, snap, and data directory.

## 12. Example Script

- **mapSensorAndFieldIdx.py**: Map the sensor name to the field point index based on the sensor's position on the ideal focal plane.
