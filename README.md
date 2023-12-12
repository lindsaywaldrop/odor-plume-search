# odor-plume-search

This is the repository with code and data to support the manuscript ``Physical properties of odorants affect behavior of trained detection dogs during close-quarters searches" currently in revision. 

A compiled copy of the manuscript is located in `doc/` as `main-manuscript.pdf` and the supplementary info is also in the `doc/` folder as `supplementary-info.pdf`. 

Data from the processed kinematic tracks required to reproduce figures in the manuscript is located in `results/csv-files` as `kinematic-final-tracks_2023-08-25.csv`. The results of other analyses are also in this folder. 

Raw data including the unprocessed kinematic tracks and ethogram files are located in the `data/` folder in `kinematics/` and `ethogram/`, respectively. 

## Instructions to reproduce figures in manuscript

This repository is design to work with RStudio (version 2023.06.0+421) as an RStudio project. It was most recently knitted under R version 4.3.1. It uses the following packages available through CRAN: 

 - bookdown (v 0.35)
 - knitr (v 1.43)
 - tidyr (v 1.3.0)
 - kableExtra (v 1.3.4)
 - ggplot2 (v 3.4.3)
 - patchwork (v 1.1.3)
 - viridis (v 0.6.4)
 - cluster (v 2.1.4)
 - here (v 1.01)
 - forcats (v 1.0.0)
 - stringr (v 1.5.0)
 - dplyr (v 1.1.3)
 - pracma (v 2.4.2)
 
Instructions:

 1. Please clone the repository and open the Rproj file in RStudio. 
 2. Source the Rscript `src/install_required_pkgs.R` to install the required packages. (Note: this will install more than the required packages above.)
 3. Open the RMD file `doc/main-manuscript.rmd`. Be sure that the project knits to the project directory instead of the document directory, or you may get an error message.
     - If you would like to reproduce only the figures, proceed to step 5.
     - If you would like to reproduce the entire analysis from raw data, please uncomment line 49 to source `kinematics_paper_analysis.R` and update line 52 so that `analysis_date` is the current date. 
 5. Knit the document. If you are not able to knit to PDF, you can change line 7 from `bookdown::pdf_document2:` to `bookdown::word_document2:` or `bookdown::html_document2:` to knit to a Microscoft Word docx or HTML file, respectively.


If you encounter problems reproducing the figures, please contact waldrop@chapman.edu 


