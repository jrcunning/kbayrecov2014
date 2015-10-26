# Symbiosis recovery dynamics after bleaching in *Montipora capitata*

### Author: Ross Cunning
### Description:
This R project contains data and scripts to analyze recovery dynamics of *Symbiodinium* in *Montipora capitata* corals following a thermal stress event in Kaneohe Bay, Oahu (Hawaii, USA) in late summer 2014. Corals were tagged and resampled 6 times between October 2014 and May 2015, and extracted DNA from each sample was used to measure symbiont to host cell ratios for both clade C and clade D *Symbiodinium*.

### Contents:
**data/qPCR/:** Directory containing .csv files of symbiont and host quantification data exported directly from Applied Biosystems StepOnePlus Software qPCR platform.

**data/coast_n83.shp/:** Hawaii coastline shapefile, accessed from [here](http://files.hawaii.gov/dbedt/op/gis/data/coast_n83.shp.zip)

**data/bleachedpair.png:** Photograph of a bleached and non-bleached pair of _M. capitata_ colonies used in Figure 1 (photo: R. Ritson-Williams)

**data/supp/:** Directory containing data for Supplementary Methods

**setup.R:** R script that imports and quality controls qPCR data in preparation for analysis

**analysis.R:** R script for all data analysis and generation of figures

**suppmethods.R:** R script for all analyses presented in the Supplementary Methods (fluorescence normalization and copy number estimation)
