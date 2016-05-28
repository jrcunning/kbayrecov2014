[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.53231.svg)](http://dx.doi.org/10.5281/zenodo.53231)

This repository includes data and analysis scripts to accompany:

# Patterns of bleaching and recovery of _Montipora capitata_ in Kāneʻohe Bay, Hawaiʻi, USA
### Authors: Ross Cunning, Raphael Ritson-Williams, Ruth D. Gates
### Journal: _Marine Ecology Progress Series_
### Link: [doi:10.3354/meps11733](http://dx.doi.org/10.3354/meps11733)  

-----

### Description:
This work analyzes bleaching and recovery dynamics of *Symbiodinium* in *Montipora capitata* corals following a thermal stress event in Kāneʻohe Bay, Oʻahu (Hawaiʻi, USA) in late 2014. Corals were tagged and resampled 6 times between October 2014 and May 2015, and extracted DNA from each sample was used to measure symbiont to host cell ratios for both clade C and clade D *Symbiodinium*.

### Contents:
#### Scripts:
* **setup.R:** R script that imports and quality controls qPCR data in preparation for analysis.

* **analysis.R:** R script for all data analyses and figures presented in the manuscript.

* **suppmethods.R:** R script for analyses presented in the Supplement (fluorescence normalization, copy number estimation, ITS2 analysis, environmental data analysis).

* **ITS2_analysis.txt:** Shell script for bioinformatic analysis of _Symbiodinium_ ITS2 data.

#### Data:
* **data/coast_n83.shp/:** Hawaii coastline shapefile, originally downloaded from [here](http://files.hawaii.gov/dbedt/op/gis/data/coast_n83.shp.zip)

* **data/bleachedpair.png:** Photograph of a bleached and non-bleached pair of _M. capitata_ colonies used in Figure 1 (photo credit: R. Ritson-Williams).

* **data/qPCR/:** Directory containing .csv files of symbiont and host quantification data exported directly from Applied Biosystems StepOnePlus Software qPCR platform.

* **data/ITS2/:** Directory containing ITS2 sequence data .fastq files (raw data) and OTU tables (final result of bioinformatic analysis). Note that intermediate/temporary data files in bioinformatic analysis are not included in the repository.

* **data/supp/:** Directory containing data for fluorescence normalization and gene copy number estimation, presented in the Supplement.

* **data/temp_light/:** Directory containing temperature and light data recorded during the study, originally downloaded from Zenodo.
    * Links to Zenodo data repositories:
        + Temperature data: [doi:10.5281/zenodo.53226](http://dx.doi.org/10.5281/zenodo.53226)
        + Light data: [doi:10.5281/zenodo.53227](http://dx.doi.org/10.5281/zenodo.53227)

#### Output:
* **output/:** Directory containing the figures and tables produced by the analysis scripts and presented in the manuscript and the Supplement.
