# ELECTRA Changelog
===================
Record of changes for ELECTRA project.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
ELECTRA uses a vMAJOR.MINOR.PATCH versioning scheme. The
* MAJOR version is increased when functionality that is not backwards compatible is added, the
* MINOR version is increased when backwards compatible functionality is added, and the
* PATCH version is increased when bug fixes are added.

## Version History
------------------
* v0.5.2
* v0.5.1
* v0.5.0
* v0.4.5
* v0.4.4
* v0.4.3
* v0.4.2
* v0.4.1
* v0.4.0
* v0.3.1
* v0.3.0
* v0.2.0
* v0.1.0


## 2022-02-20: v0.5.2
---------------------
### Fixed
* Corrected diffusion transition nodes of conduction system. Transition branches were 10 nodes long before. Now they are 10 nodes long except if a bifurcation is
  found. In this case the last transition node is the one located at the bifurcation.



## 2022-02-07: v0.5.1
---------------------
### Fixed
* Corrected an issue with the determination of fiber direction unit vectors. Now Electra normalizes the given fiber direction vectors internally


## 2022-01-02: v0.5.0
---------------------
### Added
* New applications: ElectraSim and ElectraPre

### Added
* Old application: ELECTRA-console

### Changed
* Json script interface. Check tutorial and examples folders for more info


## 2021-03-04: v0.4.5
---------------------
### Added
* Maleckar2009 electrophysiology model (Mathematical simulations of ligand-gated and cell-type specific effects on the action potential of human atrium, Maleckar et al. 2008, Progress in biophysics and molecular biology)

### Changed
* Included modifications on the IKCa current of the Courtemanche1998 model according to Engel and Peñaranda provided by Chiara Celotto

### Fixed
* Improved integration of gate variable equations in Courtemanche1998 model
* Improved integration of gate variable equations in Bueno2008 model
* Improved integration of gate variable equations in Stewart2009 model
* Adaptive diffusion timestep could exceed the maximum timestep in some cases. Now this is corrected.


## 2021-01-21: v0.4.4
---------------------
### Added
* Modified version of the Gong2020 electrophysiology model using the gate variable equations for INa from Tentusscher2006 model

### Fixed
* Improved integration of gate variable equations in Grandi2011a model


## 2021-01-20: v0.4.3
---------------------
### Added
* Tentusscher2006 electrophysiology model (Alternans and spiral breakup in a human ventricular tissue model , Tentusscher et al. 2006,  Am J Physiol Heart Circ Physiol)

### Fixed
* Improved integration of gate variable equations in Gong2020 model using the Rush Larsen method


## 2021-01-01: v0.4.2
---------------------
### Added
* Gong2020 electrophysiology model (Quantitative analysis of variability in an integrated model of human ventricular electrophysiology and b-adrenergic signaling, Gong et al. 2020, JMCC)

### Changed
* Objects *ap models* and *ap type* in input script (.json) are renamed to *ep models* and *ep type*
* electrophysiology model *Ohara* is renamed to *Ohara2011m*
* electrophysiology model *Bueno* is renamed to *Bueno2008*
* electrophysiology model *PaciVentri* is renamed to *Paci2013v*
* electrophysiology model *Courtemanche* is renamed to *Courtemanche1998*
* electrophysiology model *GrandiAtri* is renamed to *Grandi2011a*
* electrophysiology model *Maccannell* is renamed to *Maccannell2007*
* electrophysiology model *Stewart* is renamed to *Stewart2009*


## 2020-07-19: v0.4.1
---------------------
### Added
* Assignment of conduction system tree diffusivity on nodesets


## 2020-06-16: v0.4.0
---------------------
### Added
* Binary output of cells state at the end of the simulation
* Ability to start the simulation with the cells state initialized by loading a binary file
* Implementation of the ventricular conduction system

### Changed
* Significant changes on the Electra interface. Now there are two sections for tissue and conduction system setup
* Changed the conductivity attribute in the material section to diffusivity

### Fixed
* Corrections in meshfree implementation
* Corrected fibers assignment for homogeneous fiber orientation


## 2020-06-08: v0.3.1
---------------------
### Added
* Paci et al 2013 action potential model for stem cell derived human ventricular cell

### Changed
* The sign of the stimulus amplitude is now reversed. A positive amplitude value is required to trigger cell depolarization

### Fixed
* Corrections in meshfree implementation
* Corrected fibers assignment for homogeneous fiber orientation


## 2020-05-19: v0.3.0
---------------------
### Added
* Support for immersed grid Mixed Collocataion meshfree models
* Support for different dilatation coefficient for surface and interior nodes for meshfree support domains
* Calculation of APD in post processing

### Fixed
* Correction in cAF parameter application in Courtemanche action potential model


## 2020-03-19: v0.2.0
---------------------
### Added
* Time-varying action potential model parameters using a loading curve

### Changed
* Replacement of [Armadillo](http://arma.sourceforge.net/) library dependency with [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)


## 2020-02-16: v0.1.0
---------------------
* Experimental private release
