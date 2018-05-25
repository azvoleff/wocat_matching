library(rgdal)
library(rgeos)
library(raster)
library(foreach)
library(gdalUtils)
library(dplyr)
library(geoMatch)
library(ggplot2)

data_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Data'

###############################################################################
# Setup rasters

in_files <- list.files(data_folder, 'wocat_covariates_1km-',
                       full.names=TRUE)
vrt_file <- tempfile(fileext='vrt')
gdalbuildvrt(in_files, vrt_file)
r <- stack(vrt_file)
names(r) <- c('lpd_te7cl_v2', 'te_lulc', 'te_eleva', 'te_slop',
              'te_te_clima', 'te_ppt', 'te_acces', 'te_pop')

###############################################################################
# Pull control points from rasters
set.seed(932)
TOTAL_PTS <- 30000

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
                     te_te_clima, te_ppt, te_acces, te_pop)

###############################################################################
# Perform matching

# TODO: Filter by date as well

d_land_deg <- filter(d@data, (p01_imprprod == 1) | (p02_redldegr == 1)) %>%
    select(treatment, x=lon, y=lat, lpd_te7cl_v2,
           iso, te_lulc, te_eleva, te_slop,
           te_te_clima, te_ppt, te_acces, te_pop) %>%
    bind_rows(s_filtered) %>%
    filter(complete.cases(.), lpd_te7cl_v2 >= -3)
d_land_deg$lpd_te7cl_v2 <- factor(d_land_deg$lpd_te7cl_v2)

table(d_land_deg$treatment)

m <- matchit(treatment ~ iso, te_lulc, te_eleva, te_slop,
                         te_te_clima, te_ppt, te_acces, te_pop,
             method = "optimal")
ggplot(match.data(m)) +
    geom_histogram(aes(lpd_te7cl_v2, ..count../sum(..count..),
                       fill=lpd_te7cl_v2), stat='count') +
    facet_grid(treatment~.) +
    xlab('Land productivity') +
    ylab('Frequency') +
    scale_fill_discrete() +
    guides(fill=FALSE)
                 




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
