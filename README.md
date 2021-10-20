# CPSPEC manual
<!-- #### October 20, 2021 -->
<!-- #### Tenyo Kawamura (Kavli IPMU, University of Tokyo) -->

## Introduction
**CPSPEC** is an astrophysical data reduction software for timing.
Various timing properties, such as power spectra and cross spectra, are calculated from raw data.
The output data are currently written with plane text. 
However, they will be written in fits format which **XSPEC** can read, and the reasons of this change in format are listed as follows: 
- **XSPEC** are extremely useful as a fitting software. 
	Although it is designed for spectral fitting, it can be used for timing fitting as well.
- If users are aimed at simultaneous spectral-timing fitting, it is beneficial to take care of both spectral and timing data within **XSPEC**.

## Prerequisite
To use **HEASoft**, users are required to install following software:
- **HEASoft**
- **Python** with standard packages + **Astropy**

## Caveat
Obviously, the software still has room for improvement.
We strongly recommend that users check the results and would like to hear any issues found.

We summarize the consistensy of the software with others.
The software has been applied for NICER and NuSTAR data.
We have confirmed that power spectra calculated with the software are consistent with those obtained from **powspec**, and that lag-energy and lag-frequency spectra are in good agreement with those in previous papers.
