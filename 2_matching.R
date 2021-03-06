library(MASS)
library(foreach)
library(dplyr)
library(ggplot2)
library(optmatch)
library(doParallel)
library(RItools)
library(rgeos)
library(rlang)
library(sp)
library(rgdal)
library(tidyr)
registerDoParallel(4)

options("optmatch_max_problem_size"=Inf)

data_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Data'
plot_folder <- 'C:/Users/azvol/Code/LandDegradation/WOCAT/matching/Plots'

###############################################################################
# Load point data

load(file.path(data_folder, 'input_data.RData'))

# Match on initial performance as a numeric so the ordering is accounted for
d$perf_initial <- as.numeric(d$perf_initial)

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
    
get_odds <- function(m) {
    r <- polr(lpd ~ treatment, data=m, Hess=TRUE)
    bounds <- round(exp(confint(r)), 2)
    return(sprintf('odds = %.2f (%.2f, %2.2f)',
                   round(exp(coef(r)), 2), bounds[1], bounds[2]))
}

plot_adjacent <- function(m) {
    p <- group_by(m, treatment) %>%
        mutate(n=n()) %>%
        group_by(lpd, treatment) %>%
        summarise(frac=n()/n[1]) %>%
        ungroup() %>%
        ggplot() +
        geom_histogram(aes(lpd, frac,
                           fill=lpd,
                           colour=treatment,
                           size=treatment),
                       width=.5,
                       stat='identity',
                       position=position_dodge(.7)) +
        geom_text(aes(x=.85, y=.5,
                      label=paste0('n = ',
                                   nrow(filter(m, treatment)))),
                  size=2, hjust=0) +
        geom_text(aes(x=.85, y=.45, label=get_odds(m)),
                      size=2, hjust=0)+
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#ffffe0', '#73c374', '#45a146')) +
        #scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
        scale_size_manual(element_blank(),
                          values=c(.15, .5),
                          breaks=c(TRUE, FALSE),
                          labels=c('SLM', 'Control')) +
        scale_colour_manual(element_blank(),
                            values=c('grey', 'black'),
                            breaks=c(TRUE, FALSE),
                            labels=c('SLM', 'Control')) +
        guides(fill=guide_legend('Land productivity'),
               color=guide_legend(override.aes=list(fill=NA))) +
        ylim(c(0, .55)) +
        theme_bw(base_size=8) +
        theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              legend.key.size=unit(.13, 'in'),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
    return(p)
}

# Function to allow rbinding dataframes with foreach even when some dataframes 
# may not have any rows
foreach_rbind <- function(d1, d2) {
    if (is.null(d1) & is.null(d2)) {
        return(NULL)
    } else if (!is.null(d1) & is.null(d2)) {
        return(d1)
    } else if (is.null(d1) & !is.null(d2)) {
        return(d2)
    } else  {
        return(rbind(d1, d2))
    }
}

match_wocat <- function(d) {
    # Filter out countries without at least one treatment unit or without at
    # least one control unit
    d <- d %>%
        filter(complete.cases(.)) %>%
        group_by(iso) %>%
        mutate(n_treatment=sum(treatment),
               n_control=sum(!treatment)) %>%
        filter(n_treatment >= 1, n_control >= 1)

    # Note custom combine to handle iterations that don't return any value
    ret <- foreach (this_iso=unique(d$iso),
                    .packages=c('optmatch', 'dplyr'),
                    .combine=foreach_rbind, .inorder=FALSE) %dopar% {
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
        # If the controls are too far from the treatments (due to the caliper) 
        # then the matching may fail. Can test for this by seeing if subdim 
        # runs successfully
        subdim_works <- tryCatch(is.data.frame(subdim(dists)),
                                 error=function(e) return(FALSE))
        if (subdim_works) {
            m <- fullmatch(dists, min.controls=1, max.controls=1, data=this_d)
            this_d <- this_d[matched(m), ]
        } else {
            this_d <- data.frame()
        }
        # Need to handle the possibility that there were no matches for this 
        # treatment, meadning this_d will be an empty data.frame
        if (nrow(this_d) == 0) {
            return(NULL)
        } else {
            return(this_d)
        }
    }
    return(ret)
}

plot_combined <- function(m) {
    # Note the below is a loop as dplyr wasn't working for some reason...
    m <- ungroup(m)
    texts <- foreach(this_grouping=unique(m$grouping),
                     .combine=rbind) %do% {
        filter(m, grouping == this_grouping) %>%
            summarise(grouping=this_grouping,
                      n=sum(treatment),
                      odds=get_odds(.))
    }
    texts$grouping <- factor(texts$grouping, labels=levels(m$grouping))
    d <- group_by(m, grouping, treatment) %>%
        mutate(n=n()) %>%
        group_by(grouping, lpd, treatment) %>%
        summarise(frac=n()/n[1])
    p <- ggplot(d) +
        geom_histogram(aes(lpd, frac,
                           fill=lpd,
                           colour=treatment,
                           size=treatment),
                       width=.5,
                       stat='identity',
                       position=position_dodge(.7)) +
        facet_wrap(~ grouping) +
        geom_text(data=texts, aes(x=.85, y=.5, label=paste0('n = ', n, '; ', odds)),
                  size=2, hjust=0) +
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#ffffe0', '#73c374', 
                                   '#45a146')) +
        scale_size_manual(element_blank(),
                          values=c(.15, .5),
                          breaks=c(TRUE, FALSE),
                          labels=c('SLM', 'Control')) +
        scale_colour_manual(element_blank(),
                            values=c('grey', 'black'),
                            breaks=c(TRUE, FALSE),
                            labels=c('SLM', 'Control')) +
        guides(fill=guide_legend('Land productivity'),
               color=guide_legend(override.aes=list(fill=NA))) +
        ylim(c(0, .55)) +
        theme_bw(base_size=8) +
        theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
    return(p)
}

###############################################################################
# Match by approaches

d_wocat_all <- filter(d, treatment)
d_wocat_all_sp <- SpatialPointsDataFrame(select(d_wocat_all, lon, lat),
                                         d_wocat_all,
                                         proj4string=CRS('+init=epsg:4326'))
d_wocat_all_sp$lpd <- as.character(d_wocat_all_sp$lpd)
d_wocat_all_sp$implementation_approximate_fill <- as.character(d_wocat_all_sp$implementation_approximate_fill)
writeOGR(d_wocat_all_sp, data_folder, 'wocat_points', 'ESRI Shapefile',
         overwrite=TRUE)

# All approaches
d_filt_all <- select(d, treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
m_all <- match_wocat(d_filt_all)

# Just land deg and improvement
d_filt_ld_imp <- filter(d, (p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    select(treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
m_ld_imp <- match_wocat(d_filt_ld_imp)

# Just land deg and improvement
d_filt_ld_imp_last10 <- d%>%
    filter((p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    filter((implementation_approximate_fill == '0 - 10 years') | !treatment) %>%
    select(treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
m_ld_imp_last10 <- match_wocat(d_filt_ld_imp_last10)

m_all$grouping <- 'All SLM technologies'
m_ld_imp$grouping <- 'SLM technologies addressing\nland degradation or productivity'
m_ld_imp_last10$grouping <- 'SLM technologies addressing\nland degradation or productivity\n(last 10 years only)'
m_main_plot <-bind_rows(m_all, m_ld_imp, m_ld_imp_last10)
m_main_plot$grouping <- factor(m_main_plot$grouping)
plot_combined(m_main_plot)
ggsave(file.path(plot_folder, 'approaches_all_combined.png'), width=6.5, height=2.4)

write.csv(select(m_main_plot, -n_treatment, -n_control),
                 file=file.path(data_folder, 'approaches_all_combined.csv'),
                 row.names=FALSE)

# Calculate summary stats on all data
m_all_long <- m_all %>%
    ungroup() %>%
    select(treatment, iso, land_cover, elevation, slope, 
           ppt, climate, access, pop, perf_initial, lpd) %>%
    gather(key='variable', value='value', -treatment)
m_all_long$variable <- factor(m_all_long$variable,
                              levels=c('access',
                                       'climate',
                                       'elevation',
                                       'iso',
                                       'land_cover',
                                       'lpd',
                                       'perf_initial',
                                       'pop',
                                       'ppt',
                                       'slope'),
                              labels=c('log(Accessibility (min))',
                                       'Climate zone',
                                       'log(Elevation (m))',
                                       'Country',
                                       'Land cover',
                                       'LPD',
                                       'Initial performance',
                                       'log(Population (1000s))',
                                       'log(Precipitation (mm))',
                                       'log(Slope (degrees))'))
m_all_long$treatment<- ordered(m_all_long$treatment,
                               levels=c(TRUE, FALSE),
                               labels=c('SLM', 'Control'))

# Plot summary histograms for the variables that are not matched on exactly
excluded_vars <- c('Country',
                   'Land cover',
                   'Initial performance',
                   'LPD',
                   'Climate zone')
filter(m_all_long, !(variable %in% excluded_vars)) %>%
    mutate(value=as.numeric(value)) %>%
    ggplot() +
    facet_grid(treatment~variable) +
    geom_histogram(aes(value, ..count../sum(..count..)*5)) +
    ylab('Frequency') +
    xlab('Value') +
    theme_bw(base_size=8) +
    ggsave(file.path(plot_folder, 'summary_histograms_logged_vars.png'),
           width=6.5, height=3)

scalings <- data.frame(variable=c('Accessibility (min)',
                                  'Elevation (m)',
                                  'Population (1000s)',
                                  'Precipitation (mm)',
                                  'Slope (degrees)'),
                       scaling=c(10,
                                 100,
                                 1,
                                 100,
                                 1),
                       stringsAsFactors=FALSE)

# Plot summary histograms on original scale for the variables that are not
# matched on exactly
filter(m_all_long, !(variable %in% excluded_vars)) %>%
    mutate(value=exp(as.numeric(value)),
           variable=gsub('\\)$', '', gsub('log\\(', '', as.character(variable)))) %>%
    group_by(variable) %>%
    mutate(value=value * scalings$scaling[scalings$variable == variable[1]]) %>%
    ungroup() %>%
    ggplot() +
    facet_grid(treatment~variable, scales='free') +
    geom_histogram(aes(value, ..count../sum(..count..)*5)) +
    ylab('Frequency') +
    xlab('Value') +
    theme_bw(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    ggsave(file.path(plot_folder, 'summary_histograms.png'),
           width=6.5, height=3)

# Make summary table of values after logging
filter(m_all_long, !(variable %in% excluded_vars)) %>%
    mutate(value=exp(as.numeric(value)),
           variable=gsub('\\)$', '', gsub('log\\(', '', as.character(variable)))) %>%
    group_by(variable) %>%
    mutate(value=value * scalings$scaling[scalings$variable == variable[1]]) %>%
    summarise(min(value), max(value), mean(value), median(value), sd(value)) %>%
    write.csv(file.path(plot_folder, file='summary_stats.csv'),
              row.names=FALSE)

# Make summary table of values in original_units
filter(m_all_long, !(variable %in% excluded_vars)) %>%
    mutate(value=as.numeric(value)) %>%
    group_by(variable) %>%
    summarise(min(value), max(value), mean(value), median(value), sd(value)) %>%
    write.csv(file.path(plot_folder, file='summary_stats_logged.csv'),
              row.names=FALSE)

###############################################################################
# Run matches by group

match_by_group <- function(d, pattern, group_name) {
    var_indices <- which(grepl(pattern, names(d)))
    ret <- foreach (i=var_indices, .combine=foreach_rbind,
                    .inorder=FALSE) %do% {
        var_name <- names(d)[i]
        d_filt <- filter(d, (!!sym(var_name) == 1) | !treatment) %>%
            select(treatment, iso, land_cover, elevation, slope,
                   ppt, climate, access, pop, perf_initial, lpd)
        if (sum(d_filt$treatment) < 50) {
            return(NULL)
        } else {
            this_m <- match_wocat(d_filt)
            this_m$group <- group_name
            this_m$grouping <- var_name
            return(this_m)
        }
    }
    return(ret)
}

# By slm_group
m_slmgroup <- match_by_group(d, 's[0-9]{1,2}', 'slm_group')
# By management measure
m_mgtmeasure <- match_by_group(d, 'm[0-9]{1,2}', 'man_measure')
# By objective
m_degobjective <- match_by_group(d, 'd[0-9]{1,2}', 'deg_objective')
# By whether tech is preventing, reducing avoiding, restoring
m_prevention <- match_by_group(d, 'r[0-9]{1,2}', 'prevention')

# #### Now only last 10 years
# d_last10 <- filter(d, (implementation_approximate_fill == '0 - 10 years') | !treatment)
# # By slm_group
# m_slmgroup_last10 <- match_by_group(d_last10, 's[0-9]{1,2}', 'slm_group')
# # By management measure
# m_mgtmeasure_last10 <- match_by_group(d_last10, 'm[0-9]{1,2}', 'man_measure')
# # By objective
# m_degobjective_last10 <- match_by_group(d_last10, 'd[0-9]{1,2}', 'deg_objective')
# # By whether tech is preventing, reducing avoiding, restoring
# m_prevention_last10 <- match_by_group(d_last10, 'r[0-9]{1,2}', 'prevention')

###############################################################################
# Make combined plots for paper

m_slmgroup_relabeled <- filter(m_slmgroup,
                               grouping %in% c('s03_agrofore',
                                               's07_grazingm',
                                               's09_impcover',
                                               's10_minsoild',
                                               's11_soilfert',
                                               's12_slopeman'))
m_slmgroup_relabeled$grouping <- ordered(m_slmgroup_relabeled$grouping,
                                  levels=c('s03_agrofore',
                                           's07_grazingm',
                                           's09_impcover',
                                           's10_minsoild',
                                           's11_soilfert',
                                           's12_slopeman'),
                                  labels=c('Agroforestry',
                                           'Grazing management',
                                           'Improved cover',
                                           'Minimal soil disturbance',
                                           'Soil fertility management',
                                           'Slope management'))
plot_combined(m_slmgroup_relabeled)
ggsave(file.path(plot_folder, 'slmgroup.png'), width=6.5, height=3)

m_mgtmeasure_relabeled <- m_mgtmeasure
m_mgtmeasure_relabeled$grouping <- ordered(m_mgtmeasure$grouping,
                                  levels=c('m01_agronomm',
                                           'm02_vegetatm',
                                           'm03_structum',
                                           'm04_managemm'),
                                  labels=c('Agronomic',
                                           'Vegetative',
                                           'Structural',
                                           'Management'))
plot_combined(m_mgtmeasure_relabeled)
ggsave(file.path(plot_folder, 'mgtmeasures.png'), width=4.5, height=3)

m_degobjective_relabeled <- m_degobjective
m_degobjective_relabeled$grouping <- ordered(m_degobjective$grouping,
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
ggsave(file.path(plot_folder, 'degobjectives.png'), width=6.5, height=3)

m_prevention_relabeled <- m_prevention
m_prevention_relabeled$grouping <- ordered(m_prevention$grouping,
                                  levels=c('r01_preventl',
                                           'r02_reduceld',
                                           'r03_restorel'),
                                  labels=c('Prevent', 'Reduce', 'Restore'))
plot_combined(m_prevention_relabeled)
ggsave(file.path(plot_folder, 'prevention.png'), width=6.5, height=2.25)


# ###
# ### Last 10
# ###
# m_slmgroup_relabeled_last10 <- filter(m_slmgroup_last10,
#                                grouping %in% c('s03_agrofore',
#                                                's07_grazingm',
#                                                's09_impcover',
#                                                's10_minsoild',
#                                                's11_soilfert',
#                                                's12_slopeman',
#                                                's15_waterhar'))
# m_slmgroup_relabeled_last10$grouping <- ordered(m_slmgroup_relabeled_last10$grouping,
#                                   levels=c('s03_agrofore',
#                                            's07_grazingm',
#                                            's09_impcover',
#                                            's10_minsoild',
#                                            's11_soilfert',
#                                            's12_slopeman',
#                                            's15_waterhar'),
#                                   labels=c('Agroforestry',
#                                            'Grazing management',
#                                            'Improved cover',
#                                            'Minimal soil disturbance',
#                                            'Soil fertility management',
#                                            'Slope management',
#                                            'Water harvesting'))
# plot_combined(m_slmgroup_relabeled_last10)
# ggsave(file.path(plot_folder, 'slmgroup_last10.png'), width=6.5, height=6.5)
# 
# m_mgtmeasure_relabeled_last10 <- m_mgtmeasure_last10
# m_mgtmeasure_relabeled_last10$grouping <- ordered(m_mgtmeasure$grouping,
#                                   levels=c('m01_agronomm',
#                                            'm02_vegetatm',
#                                            'm03_structum',
#                                            'm04_managemm'),
#                                   labels=c('Agronomic',
#                                            'Vegetative',
#                                            'Structural',
#                                            'Management'))
# plot_combined(m_mgtmeasure_relabeled_last10)
# ggsave(file.path(plot_folder, 'mgtmeasures_last10.png'), width=4.5, height=4)
# 
# m_degobjective_relabeled_last10 <- m_degobjective_last10
# m_degobjective_relabeled_last10$grouping <- ordered(m_degobjective$grouping,
#                                   levels=c('d01_winderos',
#                                            'd02_waterero',
#                                            'd03_chemical',
#                                            'd04_physical',
#                                            'd05_biologic',
#                                            'd06_waterdeg'),
#                                   labels=c('Wind erosion',
#                                            'Water erosion',
#                                            'Chemical deterioration',
#                                            'Physical deterioration',
#                                            'Biological degradation',
#                                            'Water degradation'))
# plot_combined(m_degobjective_relabeled_last10)
# ggsave(file.path(plot_folder, 'degobjectives_last10.png'), width=6, height=4)
# 
# m_prevention_relabeled_last10 <- m_prevention_last10
# m_prevention_relabeled_last10$grouping <- ordered(m_prevention$grouping,
#                                   levels=c('r01_preventl',
#                                            'r02_reduceld',
#                                            'r03_restorel'),
#                                   labels=c('Prevent', 'Reduce', 'Restore'))
# plot_combined(m_prevention_relabeled_last10)
# ggsave(file.path(plot_folder, 'prevention_last10.png'), width=6, height=3)


save.image(file='2_matching_output.Rdata')
