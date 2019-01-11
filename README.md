# Wavefront Estimation Pipeline (WEP)

*This module is used to calculate the wavefront error in annular Zernike polynomials up to 22 terms (z4-z22) based on the intra- and extra-focal images in the large synoptic survey telescope (LSST). The main idea is to use the transport of intensity (TIE) and assume the change of intensity only comes from the wavefront error.*

## 1. Version History

*Version 1.0*
<br/>
*Finish the WEP in totally ideal condition.*
<br/>
<br/>
*Version 1.0.1*
<br/>
*Integrate the DM cmd task and implement the high-level WepController class.*
<br/>
<br/>
*Version 1.1.0*
<br/>
*Update the WEP to use the obs_lsst and scientific pipeline of sims_w_2018_47.*
<br/>

*Author: Te-Wei Tsai*
<br/>
*Date: 1-11-2019*

## 2. Platform

- *python: 3.6.6*
- *scientific pipeline (newinstall.sh from master branch)*

## 3. Needed Package

- *lsst_sims (-t sims_w_2018_47)*
- *lsst_distrib (-t w_2018_47)*
- *obs_lsst - master branch*
- *phosim_utils - master branch*
- *scikit-image*

## 4. Compile cwfs
*To compile the code, at the directory of WEP, execute:*
```
python builder/setup.py build_ext --build-lib python/lsst/ts/wep/cwfs/lib
```

## 5. Install the LSST Packages, obs_lsst, and phosim_utils

*1. Setup the LSST environment by `source $LSST_DIR/loadLSST.bash`. LSST_DIR is the directory of scientific pipeline.*
<br/>
*2. Install the lsst_sims by `eups distrib install lsst_sims -t sims_w_2018_47`.*
<br/>
*3. Install the lsst_distrib by `eups distrib install lsst_distrib -t w_2018_47`.*
<br/>
*4. Fix the path by `curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python`. The [shebangtron repo](https://github.com/lsst/shebangtron) has the further discussion for this.*
<br/>
*5. Clone the repository of [obs_lsst](https://github.com/lsst/obs_lsst) to some other directory. Under the obs_lsst directory, use `setup -k -r .` to setup the package in eups and use `scons` to build the module. It is noted that the build process is only needed for the first time.*
<br/>
*6. Do the step 5 for the repository of [phosim_utils](https://github.com/lsst-dm/phosim_utils.git).*

## 6. Pull the Built Image from Docker Hub

*Pull the built docker image by `docker pull lsstts/aos:w_2018_47`. The scientific pipeline and lsst packages are installed already. For the details of docker image, please follow the [docker aos image](https://hub.docker.com/r/lsstts/aos).*

## 7. DM Command Line Task (obs_lsst and phosim_utils)

*1. Make the faked flat images. Flats only need to be made once. They can then be shared between repos.*
```
cd $work_dir
mkdir fake_flats
cd fake_flats/
makeGainImages.py
cd ..
```
*2. Repackage the PhoSim output amplifers. The data needs to be put in single 16 extension MEFs (Multi-Extension FITS) for processing.*
```
phosim_repackager.py $phosim_amp_dir --out_dir=repackaged_files
```
*3. Make the repository, ingest the images, and ingest the calibs.*
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
*5. Run ISR.*
```
runIsr.py input --id --rerun=run1 --configfile isr_config.py
```

## 8. Use of Module

*1. Setup the DM environment:*
```
source $path_of_lsst_scientific_pipeline/loadLSST.bash
setup sims_catUtils -t sims_w_2018_47
```
*2. Setup the WEP environment:*
```
export PYTHONPATH=$PYTHONPATH:$path_to_ts_tcs_wep/python
```

## 9. Integrate with SAL

*Some environment paths defined in ts_sal/setup.env need to be modified to use lsst stack with SAL.*

*Need to setup the following path variables: LSST_SDK_INSTALL, OSPL_HOME, PYTHON_BUILD_VERSION, and PYTHON_BUILD_LOCATION.*

*1. `PYTHON_BUILD_LOCATION=$lsst_stack_python_directory`. e.g. `PYTHON_BUILD_LOCATION=/home/ttsai/Document/lsst14/python/miniconda3-4.3.21`.*

*2. In ts_sal/setup.env, use `LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SAL_HOME}/lib` instead of `LD_LIBRARY_PATH=${SAL_HOME}/lib`.*

## 10. SAL XML Model

*The SAL xml model are in the 'sal_interfaces' of ts_xml repository. The branch used is the develop branch. The CSC keywords are 'tcsOfc' and 'tcsWEP' for active optics to use. The xml files can be found in the related directory. The way to generate the SAL py libraries can follow the ts_sal manual.*

## 11. Content

*This module contains the following classes:*

- **ButlerWrapper**: Wrapper of DM butler calss to get the raw and post-ISR CCD image.
- **CamDataCollector**: Ingest the repackaged PhoSim amplifier images and fake flat calibration products based on the DM cmd task.
- **CamIsrWrapper**: Do the ISR and assemble the CCD images based on the DM cmd task.
- **SourceSelector**: Query the bright star catalog (BSC) in University of Washington (UW) to select the available target to do the WEP. If the star data is in the local database already, the query is also available.
- **SourceProcessor**: Process the post-ISR images to get the clean star images with measured optical coordinate (field x, y). The deblending algorithm is used to get the single target star image if the neighboring stars exist.
- **WfEstimator**: Calculate the wavefront error in annular Zernike polynomials up to 22 terms based on the defocal star donut images.
- **DefocalImage**: Container for the defocal images.
- **DonutImage**: Container for the donut images.
- **WepController**: High level class to use the WEP package.
- **Utility**: Utility functions used in WEP.
- **PlotUtil**: Plot utility functions used in WEP.

## 12. Example Script

- **wfsCommu.py**: Use the WEPController to issue the event and publish the telemetry.
- **mapSensorAndFieldIdx.py**: Map the sensor name to the field point index based on the sensor's position on the ideal focal plane.

## 13. Target for Future Release

- *TIE is used as the main algorithm, which is based on the single source. However, for the LSST normal case, this is not true. The initial idea here is to normalize the intensities of multiple sources.*
- *No boundary consideration of TIE studied.*
- *The use of instrument signature removal (ISR) in WEP traces to data management (DM) ISR library, which needs to customize the details/ strategies in the future release.*
- *The deblending algorithm assumes the neighboring stars have the same optical condition as the bright star. This algorithm can only handle one neighboring star that has certain magnitude and distance compared with the bright star.*
- *The algorithm to calculate the centroid of star needs a clean background.*
- *No system error is considered.*
- *No image quality determination included.*
- *No robust signal-to-noise ratio (SNR) calculation included.*
- *No master donut images by migration included.*
- *No vignette correction included.*
- *World coordinate system (WCS) is based on the focal plane with the parallax model. However, the defocal images are used in TIE. The difference and compensation between real and calculated pixel positions are not considered yet.*
- *The local BSC database is not constructed. Need to use the Scheduler to give a reasonable survey route to minimize the calculation time. Another choice is to use SkyCoord() in Astropy. The ref is at: "http://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html".*
- *The mechanism to update the BSC is not included.*
- *No statistics/ strategy of selecting wavefront sensors on full-focal plane of LSST camera included.*
- *The calculation time is much longer than the spec (14 sec).*
- *The pipeline framework (e.g. luigi) and parallel calculation are not included.*
- *The commissioning camera (ComCam) mapper is mocked based on the central raft of LSST mapper.*
- *Update the annular Zernike polynomials to Z37.*
