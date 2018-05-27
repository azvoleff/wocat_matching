library(rgdal)
library(rlang)
library(rgeos)
library(raster)
library(foreach)
library(gdalUtils)
library(dplyr)
library(geoMatch)
library(ggplot2)

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
gdalbuildvrt(c(vrts, file.path(data_folder, 'lp_perf_globe_1986_2000_avhrr_v2.tif'))
             vrt_file, resolution='highest', separate=TRUE)
r <- stack(vrt_file)
names(r) <- c('lpd_te7cl_v2', 'te_lulc', 'te_eleva', 'te_slop',
              'te_ppt', 'te_clima', 'te_acces', 'te_pop', 'initial_perf')

###############################################################################
# Pull control points from rasters
set.seed(932)
TOTAL_PTS <- 50000

# Pull a number of points from each country proportional to the number of WOCAT
# observations from that country. To do this need a count of WOCAT obs per
# country ISO code
d <- read.csv(file.path(data_folder, 'wocat_database_trendsearth_indicators_20180523_clean.csv'))
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
d_all$lpd_te7cl_v2 <- ordered(d_all$lpd_te7cl_v2, labels=c('Declining',
                                                           'Early signs of decline',
                                                           'Stable (low)',
                                                           'Stable (moderate)',
                                                           'Stable (high)',
                                                           'Early signs of increase',
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


table(d_all$lpd_te7cl_v2)
table(d_all$initial_perf, d_all$treatment)
table(d_land_deg$treatment)


options("optmatch_max_problem_size" = Inf)






# Plot color coded with date
d_all %>%
    filter(treatment) %>%
    ggplot() +
    geom_histogram(aes(lpd_te7cl_v2, fill=implementation_approximate_fill), stat='count') +
    xlab('Land productivity') +
    ylab('Frequency')
ggsave(file.path(plot_folder, 'all_data_lpd_by_implementation_year.png'))
    
plot_relative <- function(m) {
    p <- group_by(match.data(m), treatment) %>%
        mutate(n=n()) %>%
        group_by(lpd_te7cl_v2, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd_te7cl_v2, frac,
                           fill=lpd_te7cl_v2), stat='identity') +
        facet_grid(treatment~.) +
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
        guides(fill=FALSE)
    return(p)
}
                 
plot_difference <- function(m) {
# Plot the diff in relative frequencies within each class
    p <- group_by(match.data(m), treatment) %>%
        mutate(n=n()) %>%
        group_by(lpd_te7cl_v2, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        group_by(lpd_te7cl_v2) %>%
        summarise(frac_diff=frac[2] - frac[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd_te7cl_v2, frac_diff, fill=lpd_te7cl_v2), stat='identity') +
        xlab('Land productivity') +
        ylab('Difference in Relative Frequency\n(treatment - control)') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
        guides(fill=FALSE)
    return(p)
}

# TODO: drop covariates to avoid overfitting?



###############################################################################
# Approaches

# All approaches
d_filt <- select(d_all, treatment, iso, te_lulc, te_eleva, te_slop,
            te_ppt, te_clima, te_acces, te_pop, lpd_te7cl_v2, initial_perf) %>%
    filter(complete.cases(.))
m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                         te_ppt + te_clima + te_acces + te_pop + initial_perf,
             data=d_filt,
             method = "optimal")
plot_relative(m)
ggsave(file.path(plot_folder, 'approaches_all_relative.png'))
plot_difference(m)
ggsave(file.path(plot_folder, 'approaches_all_difference.png'))


# Just land deg and improvement
d_filt <- filter(d_all, (p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    select(treatment, iso, te_lulc, te_eleva, te_slop,
            te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
    filter(complete.cases(.))
m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                         te_ppt + te_clima + te_acces + te_pop + initial_perf,
             data=d_filt,
             method = "optimal")
plot_relative(m)
ggsave(file.path(plot_folder, 'approaches_ld_and_prod_relative.png'))
plot_difference(m)
ggsave(file.path(plot_folder, 'approaches_ld_and_prod_difference.png'))


###############################################################################
# By climate zone

for (zone in unique(d_all$te_clima)) {
    d_filt <- filter(d_all, (te_clima == zone) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
                te_ppt, te_acces, te_pop, lpd_te7cl_v2, initial_perf) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                             te_ppt + te_acces + te_pop + initial_perf,
                 data=d_filt,
                 method = "optimal")
    plot_relative(m)
    ggsave(file.path(plot_folder, paste0('climate_relative_zone_', zone, '.png')))
    plot_difference(m)
    ggsave(file.path(plot_folder, paste0('climate_difference_zone_', zone, '.png')))
}

###############################################################################
# By slm_group

for (slm_group in sequence(25)) {
    var_name <- names(d_all)[grepl(sprintf('s%02d', slm_group), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
                te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                             te_ppt + te_acces + te_pop + initial_perf,
                 data=d_filt,
                 method = "optimal")
    plot_relative(m)
    ggsave(file.path(plot_folder, paste0('slmgroup_relative_', var_name, '.png')))
    plot_difference(m)
    ggsave(file.path(plot_folder, paste0('slmgroup_difference_', var_name, '.png')))
}

###############################################################################
# By management measure

for (man_measure in sequence(4)) {
    var_name <- names(d_all)[grepl(sprintf('m%02d', man_measure), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
                te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                             te_ppt + te_acces + te_pop + initial_perf,
                 data=d_filt,
                 method = "optimal")
    plot_relative(m)
    ggsave(file.path(plot_folder, paste0('manmeasure_relative_', var_name, '.png')))
    plot_difference(m)
    ggsave(file.path(plot_folder, paste0('manmeasure_difference_', var_name, '.png')))
}

###############################################################################
# By objective

for (deg_objective in sequence(6)) {
    var_name <- names(d_all)[grepl(sprintf('d%02d', deg_objective), names(d_all))]
    d_filt <- filter(d_all, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, te_lulc, te_eleva, te_slop,
               te_ppt, te_clima, te_acces, te_pop, initial_perf, lpd_te7cl_v2) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + te_lulc + te_eleva + te_slop +
                     te_ppt + te_acces + te_pop + initial_perf,
                 data=d_filt,
                 method = "optimal")
    plot_relative(m)
    ggsave(file.path(plot_folder, paste0('deg_objective_relative_', var_name, '.png')))
    plot_difference(m)
    ggsave(file.path(plot_folder, paste0('deg_objective_difference_', var_name, '.png')))
}

###############################
# BELOW IS NOT USED
                 
# TODO: Temporarily reverse coords until geoMatch is fixed:
d_land_deg_sp <- SpatialPointsDataFrame(cbind(d_land_deg$y, d_land_deg$x),
                                     d_land_deg,
                                     proj4string=CRS('+init=epsg:4326'))

#match <- geoMatch(treatment ~ iso + te_eleva + te_slope + te_clima +
#                  te_acces + te_rainf + te_lulc + te_pop


match <- geoMatch(treatment ~ iso + te_eleva,
                  method = "nearest", 
                  caliper=0.25, 
                  data = d_land_deg_sp, 
                  outcome.variable="lpd_te7cl_v2", 
                  outcome.suffix="_adjusted")