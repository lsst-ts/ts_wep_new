# Wavefront Estimation Pipeline (WEP)

*This module calculates the wavefront error in annular Zernike polynomials up to 22 terms based on the intra- and extra-focal donut images in the large synoptic survey telescope (LSST).*

## 1. Version History

The version history is [here](./doc/devUse/VersionHistory.md).

Author: Te-Wei Tsai

## 2. Platform

- *CentOS 7*
- *python: 3.7.2*
- *scientific pipeline (newinstall.sh from master branch)*

## 3. Needed Package

- *lsst_sims (tag: sims_w_2019_15)*
- *lsst_distrib (tag: w_2019_15)*
- *phosim_utils - master branch (commit: 7b02084)*
- *scikit-image*
- *[documenteer](https://github.com/lsst-sqre/documenteer) (optional)*

## 4. Install the LSST Packages, phosim_utils, and ts_wep

*1. Setup the LSST environment by `source $LSST_DIR/loadLSST.bash`. LSST_DIR is the directory of scientific pipeline.* \
*2. Install the lsst_sims by `eups distrib install lsst_sims -t sims_w_2019_15`.* \
*3. Install the lsst_distrib by `eups distrib install lsst_distrib -t w_2019_15`.* \
*4. Fix the path by `curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python`. The [shebangtron repo](https://github.com/lsst/shebangtron) has the further discussion of this.* \
*5. Clone the repository of [phosim_utils](https://github.com/lsst-dm/phosim_utils.git) to some other directory. Under the phosim_utils directory, use `setup -k -r . -t sims_w_2019_15` to setup the package in eups and use `scons` to build the module. It is noted that the build process is only needed for the first time.* \
*6. Under the directory of ts_wep, do:*

```bash
setup -k -r .
scons
```

## 5. Pull the Built Image from Docker Hub

*Pull the built docker image by `docker pull lsstts/aos:w_2019_15`. The scientific pipeline and lsst packages are installed already. For the details of docker image, please follow the [docker aos image](https://hub.docker.com/r/lsstts/aos).*

## 6. DM Command Line Task (obs_lsst and phosim_utils)

*1. Make the faked flat images. Flats only need to be made once. They can then be shared between repos. The flats can be faked with (1) all sensors, (2) corner wavefront sensors, or (3) user-specified sensor list.*

```bash
cd $work_dir
mkdir fake_flats
cd fake_flats/
makeGainImages.py
cd ..
```

*2. Repackage the PhoSim output amplifiers. The data needs to be put in single 16 extension MEFs (Multi-Extension FITS) for processing.*

```bash
phosim_repackager.py $phosim_amp_dir --out_dir=repackaged_files
```

*3. Make the repository for butler to use, ingest the images, and ingest the calibration products.*

```bash
mkdir input
echo lsst.obs.lsst.phosim.PhosimMapper > input/_mapper
ingestImages.py input repackaged_files/*.fits
ingestCalibs.py input fake_flats/* --validity 99999 --output input`
```

*4. Make the config override file to turn only flat field on.*

```bash
echo "config.isr.doBias=False
config.isr.doDark=False
config.isr.doFlat=True
config.isr.doFringe=False
config.isr.doDefect=False" >isr_config.py
```

*5. Run the instrument signature removal (ISR).*

```bash
runIsr.py input --id --rerun=run1 --configfile isr_config.py
```

## 7. Use of Module

*1. Setup the DM environment.*

```bash
source $path_of_lsst_scientific_pipeline/loadLSST.bash
cd $path_of_phosim_utils
setup -k -r . -t sims_w_2019_15
```

*2. Setup the WEP environment.*

```bash
cd $path_of_ts_wep
setup -k -r .
```

## 8. Content

*The details of cotent is [here](./doc/devUse/Content.md).*

## 9. Example Script

- **mapSensorAndFieldIdx.py**: Map the sensor name to the field point index based on the sensor's position on the ideal focal plane.
- **deblendEimg.py**: Do the deblending of eimage. The deblending algorithm now can only work on the eimage.

## 10. Build the Document

*The user can use `package-docs build` to build the documentation. The documenteer is needed for this. To clean the built documents, use `package-docs clean`. See [Building single-package documentation locally](https://developer.lsst.io/stack/building-single-package-docs.html) for further details.*
