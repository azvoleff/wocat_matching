library(rgdal)
library(rgeos)
library(foreach)
library(dplyr)
library(doParallel)
registerDoParallel(4)

data_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Data'

# Read and combine point data files
d_control <- read.csv(file.path(data_folder, 'control_lpd_covariates_20180601.csv'))
d_wocat <- read.csv(file.path(data_folder, 'wocat_lpd_covariates_20180601.csv'))

d <- bind_rows(d_control, d_wocat) %>%
    rename(iso = ISO,
           perf_initial = per,
           access = acc,
           land_cover = lco,
           slope = slo,
           elevation = dem,
           climate = cli) %>%
    filter(lpd >= 1,
           perf_initial %in% c(-1, 0, 1))
d$climate <- factor(d$climate,
                    levels=sequence(12),
                    labels = c("c01_WarmTemperateMoist",
                               "c02_WarmTemperateDry",
                               "c03_CoolTemperateMoist",
                               "c04_CoolTemperateDry",
                               "c05_PolarMoist",
                               "c06_PolarDry",
                               "c07_BorealMoist",
                               "c08_BorealDry",
                               "c09_TropicalMontane",
                               "c10_TropicalWet",
                               "c11_TropicalMoist",
                               "c12_TropicalDry"))
d$land_cover <- factor(d$land_cover)
d$perf_initial <- ordered(d$perf_initial)
d$lpd <- ordered(d$lpd)
d$iso <- factor(d$iso)
d$implementation_approximate_fill <- ordered(d$implementation_approximate_fill,
                                                 levels=c("",
                                                          "3_>50 years",
                                                          "2_>10 years",
                                                          "1_<10 years"),
                                                 labels=c('No date',
                                                          '> 50 years',
                                                          '10 - 50 years',
                                                          '0 - 10 years'))
d$treatment <- as.character(d$treatment)
d$treatment[d$treatment == 'control'] <- FALSE
d$treatment[d$treatment == 'wocat'] <- TRUE
d$treatment <- as.logical(d$treatment)

save(d, file=file.path(data_folder, 'input_data.RData'))
