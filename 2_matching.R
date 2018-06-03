library(foreach)
library(dplyr)
library(ggplot2)
library(optmatch)
library(doParallel)
library(RItools)
library(rgeos)
library(sp)
registerDoParallel(4)

options("optmatch_max_problem_size"=Inf)

data_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Data'
plot_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Plots'

###############################################################################
# Load point data

load(file.path(data_folder, 'input_data.RData'))

wocat_rows <- filter(d, treatment)

###############################################################################
# Clean data for matching

# Combine all stable categories into a single category
d$lpd[d$lpd == 3] <- 4
d$lpd[d$lpd == 5] <- 4
                                                           
d$lpd <- ordered(d$lpd, 
                 levels=c(1, 2, 4, 6, 7),
                 labels=c('Declining',
                          'Early decline',
                          'Stable',
                          'Early increase',
                          'Increasing'))

table(d$lpd)


# Plot color coded with date
d %>%
    filter(treatment) %>%
    ggplot() +
    geom_histogram(aes(lpd, fill=implementation_approximate_fill), stat='count') +
    xlab('Land productivity') +
    ylab('Frequency') +
    guides(fill=guide_legend('Time since\nimplementation')) +
    theme_bw(base_size=8)
ggsave(file.path(plot_folder, 'all_data_lpd_by_implementation_year.png'),
       width=5, height=3)
    
plot_adjacent <- function(m) {
    n <- filter(m, treatment, matched) %>%
        summarise(n=n())
    p <- group_by(m, treatment) %>%
        mutate(n=n()) %>%
        group_by(lpd, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd, frac,
                           fill=lpd,
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

# TODO: Filter by date as well

match_wocat <- function(d) {
    # Filter out countries without at least one treatment unit
    d <- d %>%
        filter(complete.cases(.)) %>%
        group_by(iso) %>%
        mutate(n_treatment=sum(treatment)) %>%
        filter(n_treatment >= 1)
    ret <- foreach (this_iso=unique(d$iso), .packages=c('optmatch', 'dplyr'),
             .combine=rbind) %dopar% {
        this_d <- filter(d, iso == this_iso)
        d_wocat <- filter(this_d, treatment)
        # Filter out climates and land covers that don't appear in the wocat
        # sample, and drop these levels from the factors
        this_d <- filter(this_d,
                    climate %in% unique(d_wocat$climate),
                    land_cover %in% unique(d_wocat$land_cover),
                    land_cover %in% unique(d_wocat$land_cover))
        this_d$climate <- droplevels(this_d$climate)
        this_d$land_cover <- droplevels(this_d$land_cover)
        f <- treatment ~ elevation + slope + ppt + access + pop + perf_initial
        # Can't stratify by land cover or climate if they only have one level
        if (nlevels(this_d$land_cover) >= 2) {
            f <- update(f, ~ . + strata(land_cover))
        } else {
            f <- update(f, ~ . - land_cover)
        }
        if (nlevels(this_d$climate) >= 2) {
            f <- update(f, ~ . + strata(climate))
        } else {
            f <- update(f, ~ . - climate)
        }
        if (nrow(d_wocat) > 2) {
            model <- glm(f, data=this_d)
            dists <- match_on(model, data=this_d)
        } else {
            # Use Mahalanobis distance if there aren't enought points to run a
            # glm
            dists <- match_on(f, data=this_d)
        }
        dists <- caliper(dists, 2)
        m <- fullmatch(dists, max.controls=1, data=this_d)
        this_d$matched <- matched(m)
        return(this_d)
    }
    return(ret)
}

###############################################################################
# Match by approaches

# Match on initial performance as a numeric so the ordering is accounted for
d$perf_initial <- as.numeric(d$perf_initial)

# All approaches
d_filt_all <- select(d, treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, lpd, perf_initial) %>%
    filter(complete.cases(.))
m_all <- match_wocat(d_filt_all)
plot_adjacent(m_all)
ggsave(file.path(plot_folder, 'approaches_all.png'), width=4, height=3)


# Just land deg and improvement
d_filt_ld_imp <- filter(d, (p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    select(treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
m_ld_imp <- match_wocat(d_filt_ld_imp)
plot_adjacent(m_ld_imp)
ggsave(file.path(plot_folder, 'approaches_ld_and_prod.png'), width=4, height=3)

###############################################################################
# Run matches by slm_group

n <- 1
for (slm_group in sequence(25)) {
    var_name <- names(d)[grepl(sprintf('s%02d', slm_group), names(d))]
    d_filt <- filter(d, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, land_cover, elevation, slope,
                ppt, climate, access, pop, perf_initial, lpd) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + land_cover + elevation + slope +
                             ppt + access + pop + perf_initial,
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
    var_name <- names(d)[grepl(sprintf('m%02d', man_measure), names(d))]
    d_filt <- filter(d, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, land_cover, elevation, slope,
                ppt, climate, access, pop, perf_initial, lpd) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + land_cover + elevation + slope +
                             ppt + access + pop + perf_initial,
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
    var_name <- names(d)[grepl(sprintf('d%02d', deg_objective), names(d))]
    d_filt <- filter(d, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, land_cover, elevation, slope,
               ppt, climate, access, pop, perf_initial, lpd) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + land_cover + elevation + slope +
                     ppt + access + pop + perf_initial,
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
    var_name <- names(d)[grepl(sprintf('r%02d', deg_objective), names(d))]
    d_filt <- filter(d, (!!sym(var_name) == 1) | !treatment) %>%
        select(treatment, iso, land_cover, elevation, slope,
               ppt, climate, access, pop, perf_initial, lpd) %>%
        filter(complete.cases(.))
    if (sum(d_filt$treatment) < 50) next
    m <- matchit(treatment ~ iso + land_cover + elevation + slope +
                     ppt + access + pop + perf_initial,
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
        group_by(variable, lpd, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ggplot() +
        geom_histogram(aes(lpd, frac,
                           fill=lpd,
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
#match <- geoMatch(treatment ~ iso + elevation + slopee + climate +
#                  access + te_rainf + land_cover + pop
# 
# 
# match <- geoMatch(treatment ~ iso + elevation,
#                   method = "nearest", 
#                   caliper=0.25, 
#                   data = d_land_deg_sp, 
#                   outcome.variable="lpd", 
#                   outcome.suffix="_adjusted")