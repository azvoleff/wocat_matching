library(rgdal)
library(rlang)
library(rgeos)
library(raster)
library(foreach)
library(gdalUtils)
library(dplyr)
library(geoMatch)
library(ggplot2)
library(optmatch)
library(doParallel)
library(RItools)
registerDoParallel(4)

options("optmatch_max_problem_size"=Inf)


data_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Data'
plot_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Plots'

###############################################################################
# Setup rasters

# First combine the tiled GEE outputs into a set of 1 band VRTs
in_files <- list.files(data_folder,
                       'wocat_covariates_1km-',
                       full.names=TRUE)
vrts <- c()
for (n in 1:8) {
    vrt <- tempfile(fileext='vrt')
    gdalbuildvrt(in_files, b=n, vrt, resolution='highest')
    vrts <- c(vrts, vrt)
}
# Now combine any other single band tifs with the mosaiced GEE output
vrt_file <- tempfile(fileext='vrt')
gdalbuildvrt(c(vrts, file.path(data_folder, 'lp_perf_globe_1986_2000_avhrr_v2.tif')),
             vrt_file, resolution='highest', separate=TRUE)
r <- stack(vrt_file)
names(r) <- c('lpd_te7cl_v2', 'te_lulc', 'te_eleva', 'te_slop',
              'te_ppt', 'te_clima', 'te_acces', 'te_pop', 'initial_perf')


###############################################################################
# Pull control points from rasters
set.seed(932)
TOTAL_PTS <- 100000

# Pull a number of points from each country proportional to the number of WOCAT
# observations from that country. To do this need a count of WOCAT obs per
# country ISO code
d <- read.csv(file.path(data_folder, 'wocat_database_trendsearth_indicators_20180527_clean.csv'))
d <- SpatialPointsDataFrame(cbind(d$lon, d$lat), d, proj4string=CRS(proj4string(r)))

names(d@data)

d@data <- select(d@data, -lpd_te7cl_v2, -starts_with('te'))

# Pull ISOs for the WOCAT data
adm <- readOGR(data_folder, 'ne_50m_admin_0_countries_split')
adm$id <- sequence(nrow(adm))
d_codes <- over(d, adm[, 'ADM0_A3'])
d$iso <- d_codes$ADM0_A3
d$treatment <- TRUE

# TODO: Need to make holes in the ADM layer to remove a buffer around each
# WOCAT observation.

# Drop countries with very few (say less than 2 or 3?) observations
d_filtered <- group_by(d@data, iso) %>%
    mutate(n=n()) %>%
    filter(n >= 10)
d <- d[d$X %in% d_filtered$X, ]

# Select potential matching points within the polygon part each WOCAT point
# comes from.
d_poly_ids <- over(d, adm[, 'id'])
d$poly_id <- d_poly_ids$id

# Filter out WOCAT observations that don't fall within a polygon at all
d <- d[!is.na(d$poly_id), ]

# Add covariates for matching
e <- extract(r, d, df=TRUE)
e <- select(e, -ID)
stopifnot(nrow(d@data) == nrow(e))
d@data <- cbind(d@data, e)
d$iso <- droplevels(d$iso)

save(d, file=file.path(data_folder, 'data_treatment.RData'))

s <- foreach (poly_id=unique(d$poly_id), .combine=rbind) %do% {
    n_obs <- sum(d$poly_id == poly_id, na.rm=TRUE)
    poly <- adm[adm$id == poly_id, ]
    r_cropped <- crop(r, poly)
    # Ensure enough points are drawn to have at least a few within the
    # country polygon - so draw at least 100.
    n_pts <- max(c(ceiling(n_obs/nrow(d) * TOTAL_PTS), 100))
    print(paste0(poly$ADM0_A3, ' (', poly_id, ', ', n_obs, '): ', n_pts))
    s <- data.frame()
    while (nrow(s) < 10) {
        s <- sampleRandom(r_cropped, n_pts, sp=TRUE)
        in_pts = gWithin(s, poly, byid=TRUE)[1, ]
        # Handle the case of no points being within the polygon
        if (sum(in_pts) == 0) {
          s <- data.frame()
          next
        }
        s <- s[in_pts, ]
    }
    s$iso <- poly$ADM0_A3
    s$treatment <- FALSE
    return(as.data.frame(s))
}
s_filtered <- select(s, treatment, x, y, lpd_te7cl_v2,
                     iso, te_lulc, te_eleva, te_slop,
                     te_ppt, te_clima, te_acces, te_pop, initial_perf)
save(s_filtered, file=file.path(data_folder, 'data_controls.RData'))

###############################################################################
# Perform matching

# TODO: Filter by date as well
d_all <-  rename(d@data, x=lon, y=lat) %>%
    bind_rows(s_filtered) %>%
    filter(lpd_te7cl_v2 >= 1)

d_all$lpd_te7cl_v2[d_all$lpd_te7cl_v2 == 3] <- 4
d_all$lpd_te7cl_v2[d_all$lpd_te7cl_v2 == 5] <- 4
                                                           
d_all$lpd_te7cl_v2 <- ordered(d_all$lpd_te7cl_v2, 
                              levels=c(1, 2, 4, 6, 7),
                              labels=c('Declining',
                                        'Early decline',
                                        'Stable',
                                        'Early increase',
                                        'Increasing'))
d_all$te_lulc <- factor(d_all$te_lulc)
d_all$iso <- factor(d_all$iso)
d_all$initial_perf <- ordered(d_all$initial_perf)
d_all$te_clima[d_all$te_clima == 1]  <- "c01_WarmTemperateMoist"
d_all$te_clima[d_all$te_clima == 2]  <- "c02_WarmTemperateDry"
d_all$te_clima[d_all$te_clima == 3]  <- "c03_CoolTemperateMoist"
d_all$te_clima[d_all$te_clima == 4]  <- "c04_CoolTemperateDry"
d_all$te_clima[d_all$te_clima == 5]  <- "c05_PolarMoist"
d_all$te_clima[d_all$te_clima == 6]  <- "c06_PolarDry"
d_all$te_clima[d_all$te_clima == 7]  <- "c07_BorealMoist"
d_all$te_clima[d_all$te_clima == 8]  <- "c08_BorealDry"
d_all$te_clima[d_all$te_clima == 9]  <- "c09_TropicalMontane"
d_all$te_clima[d_all$te_clima == 10] <- "c10_TropicalWet"
d_all$te_clima[d_all$te_clima == 11] <- "c11_TropicalMoist"
d_all$te_clima[d_all$te_clima == 12] <- "c12_TropicalDry"
d_all$te_clima <- factor(d_all$te_clima)

d_all$implementation_approximate_fill <- ordered(d_all$implementation_approximate_fill,
                                                 levels=c("",
                                                          "3_>50 years",
                                                          "2_>10 years",
                                                          "1_<10 years"),
                                                 labels=c('No date',
                                                          '> 50 years',
                                                          '10 - 50 years',
                                                          '0 - 10 years'))

save(d_all, file=file.path(data_folder, 'data_all.RData'))
write.csv(filter(d_all, treatment) %>% select(x, y), file.path(plot_folder, 'treatment_points.csv'))

#load(file=file.path(data_folder, 'data_all.RData'))

table(d_all$lpd_te7cl_v2)
table(d_all$initial_perf, d_all$treatment)
table(d_all$treatment)

# Plot color coded with date
d_all %>%
    filter(treatment) %>%
    ggplot() +
    geom_histogram(aes(lpd_te7cl_v2, fill=implementation_approximate_fill), stat='count') +
    xlab('Land productivity') +
    ylab('Frequency') +
    guides(fill=guide_legend('Time since\nimplementation')) +
    theme_bw(base_size=8)
ggsave(file.path(plot_folder, 'all_data_lpd_by_implementation_year.png'), width=5, height=3)
    
plot_adjacent <- function(m) {
    n <- filter(m, treatment) %>%
        summarise(n=n())
    p <- group_by(m, treatment) %>%
        mutate(n=n()) %>%
        group_by(lpd_te7cl_v2, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd_te7cl_v2, frac,
                           fill=lpd_te7cl_v2,
                           colour=treatment,
                           size=treatment),
                       width=.5,
                       stat='identity',
                       position=position_dodge(.7)) +
        geom_text(data=n, aes(x=1, y=Inf, label=paste0('\n    n = ', n)), size=2) +
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
        scale_size_manual(element_blank(),
                          values=c(.15, .4),
                          breaks=c(TRUE, FALSE),
                          labels=c('WOCAT', 'Control')) +
        scale_colour_manual(element_blank(),
                            values=c('grey', 'black'),
                            breaks=c(TRUE, FALSE),
                            labels=c('WOCAT', 'Control')) +
        guides(fill=guide_legend('Land productivity'),
               color=guide_legend(override.aes=list(fill=NA))) +
        theme_bw(base_size=8) +
        theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
    return(p)
}

# TODO: drop covariates to avoid overfitting?

###############################################################################
# Approaches

# All approaches
d_filt_all <- select(d_all, treatment, iso, te_lulc, te_eleva, te_slop,
            te_ppt, te_clima, te_acces, te_pop, lpd_te7cl_v2, initial_perf) %>%
    filter(complete.cases(.))



match_wocat <- function(d) {
    foreach (this_iso=unique(d$iso), .packages=c('optmatch', 'dplyr'),
             .combine=rbind, .inorder=FALSE) %dopar% {
        d <- filter(d, iso == this_iso)
        d$te_clima <- droplevels(d$te_clima)
        d$te_lulc <- droplevels(d$te_lulc)
        if ((nlevels(d$te_lulc) == 1) | (nlevels(d$te_clima) == 1)) {
            # Can't stratify by them if there is only one level for each
            model <- glm(treatment ~ te_lulc + te_clima +
                         te_eleva + te_slop + te_ppt + te_acces + te_pop +
                         initial_perf, data=d)
        if ((nlevels(d$te_lulc) == 1) | (nlevels(d$te_clima) == 1)) {
            # Can't stratify by them if there is only one level for each
            model <- glm(treatment ~ te_lulc + te_clima +
                         te_eleva + te_slop + te_ppt + te_acces + te_pop +
                         initial_perf, data=d)
        } else {
            model <- glm(treatment ~ strata(te_lulc, te_clima) +
                         te_eleva + te_slop + te_ppt + te_acces + te_pop +
                         initial_perf, data=d)
        }
        dists <- match_on(model, data=d)
        dists <- caliper(dists, 2)
        m <- pairmatch(dists, data=d)
        d$m <- m
        return(d)
    }
    return(d)
}

m_all <- match_wocat(d_filt_all)

summary(m_all)
d_filt_all$m_all <- m_all
plot_adjacent(d_filt_all[matched(m_all), ])
plot_adjacent(d_filt_all)

summary(m_all)

m_all <- match.data(m)

m_all$group <- 'all'
ggsave(file.path(plot_folder, 'approaches_all.png'), width=4, height=3)


# Just land deg and improvement
d_filt_ld_imp <- filter(d_all, (p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    select(treatment, iso, te_lulc, te_eleva, te_slop,
            te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
    filter(complete.cases(.))
m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                         te_ppt + te_clima + te_acces + te_pop + initial_perf,
             exact=c('iso'),
             data=d_filt_ld_imp,
             method = "optimal")
m_ld_imp <- match.data(m)
m_ld_imp$group <- 'ld_imp'
plot_adjacent(m_ld_imp)
ggsave(file.path(plot_folder, 'approaches_ld_and_prod.png'), width=4, height=3)

###############################################################################
# By slm_group

n <- 1
for (slm_group in sequence(25)) {
    var_name <- names(d_all)[grepl(sprintf('s%02d', slm_group), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
                te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                             te_ppt + te_acces + te_pop + initial_perf,
                 exact=c('iso'),
                 data=d_filt,
                 method = "optimal")
    m_data <- match.data(m)
    m_data$group <- 'slm_group'
    m_data$variable <- var_name
    if (n == 1) {
        m_slmgroup <- m_data
    } else {
        m_slmgroup <- bind_rows(m_slmgroup, m_data)
    }
    n <- n + 1
}

###############################################################################
# By management measure

n <- 1
for (man_measure in sequence(4)) {
    var_name <- names(d_all)[grepl(sprintf('m%02d', man_measure), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
                te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                             te_ppt + te_acces + te_pop + initial_perf,
                 exact=c('iso'),
                 data=d_filt,
                 method = "optimal")
    m_data <- match.data(m)
    m_data$group <- 'man_measure'
    m_data$variable <- var_name
    if (n == 1) {
        m_manmeasure <- m_data
    } else {
        m_manmeasure <- bind_rows(m_manmeasure, m_data)
    }
    n <- n + 1
}

###############################################################################
# By objective

n <- 1
for (deg_objective in sequence(6)) {
    var_name <- names(d_all)[grepl(sprintf('d%02d', deg_objective), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
               te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                     te_ppt + te_acces + te_pop + initial_perf,
                 exact=c('iso'),
                 data=d_filt,
                 method = "optimal")
    m_data <- match.data(m)
    m_data$group <- 'deg_objective'
    m_data$variable <- var_name
    if (n == 1) {
        m_degobjective <- m_data
    } else {
        m_degobjective <- bind_rows(m_degobjective, m_data)
    }
    n <- n + 1
}

###############################################################################
# By whether tech is preventing, reducing avoiding, restoring

n <- 1
for (deg_objective in sequence(4)) {
    var_name <- names(d_all)[grepl(sprintf('r%02d', deg_objective), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
               te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                     te_ppt + te_acces + te_pop + initial_perf,
                 exact=c('iso'),
                 data=d_filt,
                 method = "optimal")
    m_data <- match.data(m)
    m_data$group <- 'prevention'
    m_data$variable <- var_name
    if (n == 1) {
        m_prevention <- m_data
    } else {
        m_prevention <- bind_rows(m_prevention, m_data)
    }
    n <- n + 1
}

###############################################################################
# Make combined plots for paper
plot_combined <- function(m) {
    n <- group_by(m, variable) %>%
        filter(treatment) %>%
        summarise(n=n())
    p <- group_by(m, variable, treatment) %>%
        mutate(n=n()) %>%
        group_by(variable, lpd_te7cl_v2, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd_te7cl_v2, frac,
                           fill=lpd_te7cl_v2,
                           colour=treatment,
                           size=treatment),
                       width=.5,
                       stat='identity',
                       position=position_dodge(.7)) +
        facet_wrap(~ variable) +
        geom_text(data=n, aes(x=1, y=Inf, label=paste0('\n    n = ', n)), size=2) +
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
        scale_size_manual(element_blank(),
                          values=c(.15, .4),
                          breaks=c(TRUE, FALSE),
                          labels=c('WOCAT', 'Control')) +
        scale_colour_manual(element_blank(),
                            values=c('grey', 'black'),
                            breaks=c(TRUE, FALSE),
                            labels=c('WOCAT', 'Control')) +
        guides(fill=guide_legend('Land productivity'),
               color=guide_legend(override.aes=list(fill=NA))) +
        theme_bw(base_size=8) +
        theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
    return(p)
}

m_manmeasure_relabeled <- m_manmeasure
m_manmeasure_relabeled$variable <- ordered(m_manmeasure$variable,
                                  levels=c('m01_agronomm',
                                           'm02_vegetatm',
                                           'm03_structum',
                                           'm04_managemm'),
                                  labels=c('Agronomic',
                                           'Vegetative',
                                           'Structural',
                                           'Management'))
plot_combined(m_manmeasure_relabeled)
ggsave(file.path(plot_folder, 'manmeasures.png'), width=4.5, height=4)

m_slmgroup_relabeled <- m_slmgroup
m_slmgroup_relabeled$variable <- ordered(m_slmgroup$variable,
                                  levels=c('s03_agrofore',
                                           's07_grazingm',
                                           's09_impcover',
                                           's10_minsoild',
                                           's11_soilfert',
                                           's12_slopeman',
                                           's15_waterhar',
                                           's16_irrmanag',
                                           's17_drainage',
                                           's18_swateman',
                                           's19_gwateman'),
                                  labels=c('Agroforestry',
                                           'Grazing Management',
                                           'Improved cover',
                                           'Minimal soil disturbance',
                                           'Soil fertility management',
                                           'Cross-slope measures',
                                           'Water harvesting',
                                           'Irrigation management',
                                           'Drainage',
                                           'Surface water management',
                                           'Groundwater management'))
plot_combined(m_slmgroup_relabeled)
ggsave(file.path(plot_folder, 'slmgroup.png'), width=6.5, height=6.5)

m_degobjective_relabeled <- m_degobjective
m_degobjective_relabeled$variable <- ordered(m_degobjective$variable,
                                  levels=c('d01_winderos',
                                           'd02_waterero',
                                           'd03_chemical',
                                           'd04_physical',
                                           'd05_biologic',
                                           'd06_waterdeg'),
                                  labels=c('Wind erosion',
                                           'Water erosion',
                                           'Chemical deterioration',
                                           'Physical deterioration',
                                           'Biological degradation',
                                           'Water degradation'))
plot_combined(m_degobjective_relabeled)
ggsave(file.path(plot_folder, 'degobjectives.png'), width=6, height=4)

m_prevention_relabeled <- m_prevention
m_prevention_relabeled$variable <- ordered(m_prevention$variable,
                                  levels=c('r01_preventl',
                                           'r02_reduceld',
                                           'r03_restorel'),
                                  labels=c('Prevent', 'Reduce', 'Restore'))
plot_combined(m_prevention_relabeled)
ggsave(file.path(plot_folder, 'prevention.png'), width=6, height=3)

###############################
# BELOW IS NOT USED
                 
# TODO: Temporarily reverse coords until geoMatch is fixed:
# d_land_deg_sp <- SpatialPointsDataFrame(cbind(d_land_deg$y, d_land_deg$x),
#                                      d_land_deg,
#                                      proj4string=CRS('+init=epsg:4326'))
# 
#match <- geoMatch(treatment ~ iso + te_eleva + te_slope + te_clima +
#                  te_acces + te_rainf + te_lulc + te_pop
# 
# 
# match <- geoMatch(treatment ~ iso + te_eleva,
#                   method = "nearest", 
#                   caliper=0.25, 
#                   data = d_land_deg_sp, 
#                   outcome.variable="lpd_te7cl_v2", 
#                   outcome.suffix="_adjusted")