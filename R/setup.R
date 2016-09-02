## Libraries

# if(!"dplyr" %in% installed.packages()){
#   install.packages("dplyr")
# }
# if(!"tidyr" %in% installed.packages()){
#   install.packages("tidyr")
# }
# if(!"ggplot2" %in% installed.packages()){
#   install.packages("ggplot2")
# }
# if(!"binom" %in% installed.packages()){
#   install.packages("binom")
# }
# if(!"pander" %in% installed.packages()){
#   install.packages("pander")
# }
# if(!"grid" %in% installed.packages()){
#   install.packages("grid")
# }
# if(!"withr" %in% installed.packages()){
#   install.packages("withr")
# }

library(dplyr)
library(tidyr)
library(ggplot2)
#library(binom)
library(pander)
library(grid)
#library(withr)

## Install knitcitations in this repo if not already installed
# if(!"knitcitations" %in% installed.packages(lib.loc="../repoPackageLibrary")){
#   install.packages("knitcitations",lib="../repoPackageLibrary")
# }

# ## Install knitauthors in this repo if not already installed (requires devtools)
# if(!"knitauthors" %in% installed.packages(lib.loc="../repoPackageLibrary")){
#   if(!"devtools" %in% installed.packages()){
#     install.packages("devtools")
#   }
#   withr::with_libpaths(new = "../repoPackageLibrary/", 
#                           devtools::install_github("jamesware/knitauthors",force=T),
#                           action= "prefix")
# }
## Set up citations
# library("knitcitations",lib.loc="../repoPackageLibrary")
library(knitcitations)
cleanbib()
options("citation_format" = "pandoc")

## set knitr chunk options
# if(!"knitr" %in% installed.packages()){
#   install.packages("knitr")
# }
library(knitr)
opts_chunk$set(fig.path="../figures/", collapse = TRUE, echo=FALSE, warning=TRUE, message=TRUE)#, cache=TRUE)#, results="asis")
library(knitauthors)
# library("knitauthors",lib.loc="../repoPackageLibrary")
