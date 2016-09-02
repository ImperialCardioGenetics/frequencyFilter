### Exploring ExAC AC of variants from 6179 HCM cases (Walsh et al)

summariseCaseSeriesVariants <- function(){
  
  referenceVariants <<- read.delim("../data-raw/LMM_OXF_path_likelyPath_VUS_variants_allExACpopulations_withLowestClass.txt") %>%
    filter(Disease=="HCM") %>%
    filter(!(Gene %in% c('GLA','DMD','TAZ','EMD','LAMP2','VBP1','FHL1','KCNE1L','ZIC3','COL4A5','FLNA','MED12')))
  
  myVariantsN <<- nrow(referenceVariants)
  
  myFailPath <<- referenceVariants %>%
    filter(Class %in% c("Pathogenic","Likely Pathogenic")) %>%
    filter(ExACCountALL>9) %>%
    nrow
  myPassPath <<- referenceVariants %>%
    filter(Class %in% c("Pathogenic","Likely Pathogenic")) %>%
    filter(ExACCountALL<=9) %>%
    nrow
  myAbsentPath <<- referenceVariants %>%
    filter(Class %in% c("Pathogenic","Likely Pathogenic")) %>%
    filter(ExACCountALL==0) %>%
    nrow
  myAllPath <<- referenceVariants %>%
    filter(Class %in% c("Pathogenic","Likely Pathogenic")) %>%
    nrow
  myFailVus <<- referenceVariants %>%
    filter(Class=="VUS") %>%
    filter(ExACCountALL>9) %>%
    nrow
  myPassVus <<- referenceVariants %>%
    filter(Class=="VUS") %>%
    filter(ExACCountALL<=9) %>%
    nrow
  myAllVus <<- referenceVariants %>%
    filter(Class=="VUS") %>%
    nrow
  
  
  
}