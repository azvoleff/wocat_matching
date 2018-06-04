library(MASS)
library(foreach)
library(dplyr)
library(ggplot2)
library(optmatch)
library(doParallel)
library(RItools)
library(rgeos)
library(rlang)
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
    
plot_adjacent <- function(m) {
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
        geom_text(aes(x=1, y=Inf,
                      label=paste0('\n    n = ',
                                   nrow(filter(m, treatment)))), size=2) +
        xlab('Land productivity') +
        ylab('Frequency') +
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#ffffe0', '#73c374', '#45a146')) +
        #scale_fill_manual(values=c('#ab2727', '#ed7428', '#fee9a9', '#ffffe0', '#c4e6c4', '#73c374', '#45a146')) +
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
        ylim(c(0, .4)) +
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

###############################################################################
# Match by approaches

# All approaches
d_filt_all <- select(d, treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
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

# Just land deg and improvement
d_filt_ld_imp_last10 <- d%>%
    filter((p01_imprprod == 1) | (p02_redldegr == 1) | !treatment) %>%
    filter((implementation_approximate_fill == '0 - 10 years') | !treatment) %>%
    select(treatment, iso, land_cover, elevation, slope,
            ppt, climate, access, pop, perf_initial, lpd)
m_ld_imp_last10 <- match_wocat(d_filt_ld_imp_last10)
plot_adjacent(m_ld_imp_last10)
ggsave(file.path(plot_folder, 'approaches_ld_and_prod_last10.png'), width=4, height=3)

###############################################################################
# Run matches by group

match_by_group <- function(pattern, group_name) {
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
            this_m$variable <- var_name
            return(this_m)
        }
    }
    return(ret)
}

# By slm_group
m_slmgroup <- match_by_group('s[0-9]{1,2}', 'slm_group')

# By management measure
m_mgtmeasure <- match_by_group('m[0-9]{1,2}', 'man_measure')

# By objective
m_degobjective <- match_by_group('d[0-9]{1,2}', 'deg_objective')

# By whether tech is preventing, reducing avoiding, restoring
m_prevention <- match_by_group('r[0-9]{1,2}', 'deg_objective')

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
        scale_fill_manual(values=c('#ab2727', '#ed7428', '#ffffe0', '#73c374', 
                                   '#45a146')) +
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
        ylim(c(0, .5)) +
        theme_bw(base_size=8) +
        theme(panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
    return(p)
}

m_slmgroup_relabeled <- filter(m_slmgroup,
                               variable %in% c('s03_agrofore',
                                               's07_grazingm',
                                               's09_impcover',
                                               's10_minsoild',
                                               's11_soilfert',
                                               's12_slopeman',
                                               's15_waterhar'))
m_slmgroup_relabeled$variable <- ordered(m_slmgroup_relabeled$variable,
                                  levels=c('s03_agrofore',
                                           's07_grazingm',
                                           's09_impcover',
                                           's10_minsoild',
                                           's11_soilfert',
                                           's12_slopeman',
                                           's15_waterhar'),
                                  labels=c('Agroforestry',
                                           'Grazing management',
                                           'Improved cover',
                                           'Minimal soil disturbance',
                                           'Soil fertility management',
                                           'Slope management',
                                           'Water harvesting'))
plot_combined(m_slmgroup_relabeled)
ggsave(file.path(plot_folder, 'slmgroup.png'), width=6.5, height=6.5)

m_mgtmeasure_relabeled <- m_mgtmeasure
m_mgtmeasure_relabeled$variable <- ordered(m_mgtmeasure$variable,
                                  levels=c('m01_agronomm',
                                           'm02_vegetatm',
                                           'm03_structum',
                                           'm04_managemm'),
                                  labels=c('Agronomic',
                                           'Vegetative',
                                           'Structural',
                                           'Management'))
plot_combined(m_mgtmeasure_relabeled)
ggsave(file.path(plot_folder, 'mgtmeasures.png'), width=4.5, height=4)

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

dim(m_prevention_relabeled)
