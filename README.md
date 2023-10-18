# odor-plume-search

This is the repository with code and data to support the manuscript ``Physical properties of odorants affect behavior of trained detection dogs during close-quarters searches" currently in review. 

A compiled copy of the manuscript is located in `doc/` as `main-manuscript.pdf` and the supplementary info is also in the `doc/` folder as `supplementary-info.pdf`. 

Data from the processed kinematic tracks required to reproduce figures in the manuscript is located in `results/csv-files` as `kinematic-final-tracks_2023-08-25.csv`. The results of other analyses are also in this folder. 

## Instructions to reproduce figures of manuscript

This repository is design to work with RStudio as an RStudio project. It was most recently knitted under R version 4.3.1. It uses the following packages (manuscript figures only): 

 - bookdown (v 0.35)
 - knitr (v 1.43)
 - tidyr (v 1.3.0)
 - kableExtra (v 1.3.4)
 - ggplot2 (v 3.4.3)
 - patchwork (v 1.1.3)
 - viridis (v 0.6.4)
 - cluster (v 2.1.4)
 
Instructions:

 1. Please clone the repository and open the Rproj file in RStudio. Once opened, 
 2. Source the Rscript `src/install_required_pkgs.R` to install the required packages. (Note: this will install more than the required packages above.)
 3. Open the RMD file `doc/main-manuscript.rmd` and knit. Be sure that the project knits to the Project directory. 
 4. If you are not able to knit to PDF, you can change line 7 from `bookdown::pdf_document2:` to `bookdown::word_document2:` or `bookdown::html_document2:` to knit to a Microscoft Word docx or HTML file, respectively. 
 
If you encounter problems reproducing the figures, please contact waldrop@chapman.edu 


