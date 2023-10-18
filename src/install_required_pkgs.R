# Package installation script

packages <- c("markdown", "rmarkdown", "bookdown", "knitr", "pracma", "tidyr", 
              "forcats", "dplyr", "ggplot2", "patchwork", "kableExtra", "cluster",
               "viridis", "signal", "stringr", "here", "scatterplot3d")

package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)