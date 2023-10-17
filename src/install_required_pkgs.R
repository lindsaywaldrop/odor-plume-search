# Package installation script

packages <- c("markdown", "rmarkdown", "pracma", "tidyr", "forcats", "dplyr", "ggplot2", "patchwork",
               "viridis", "signal", "stringr", "here")

package.check <- lapply(
  packages,
  FUN <- function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE, repos = "http://cran.us.r-project.org")
      library(x, character.only = TRUE)
    }
  }
)