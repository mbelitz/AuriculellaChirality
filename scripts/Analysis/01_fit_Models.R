library(tidyverse)
library(brms)
library(sf)
library(terra)

## read in data
mdf_scaled <- read.csv("data/mdf_scaled_July2026Update.csv") %>% 
    mutate(n_sc = scale(n),
           mpd_sc = scale(MPD)) %>% 
    filter(correctedspeciesname != "Auriculella cerea")

## generate site IDs
uniqueSites <- distinct(mdf_scaled, decimallatitudeUSE, decimallongitudeUSE)

uniqueSites <- uniqueSites %>% 
    mutate(siteID = 1:nrow(uniqueSites))

mdf_scaled <- mdf_scaled %>% 
    left_join(uniqueSites)

mdf_scaled <- mdf_scaled %>% 
    mutate(is_monomorphic = if_else(condition = dex == 0 | dex == n,
                                    true = 1, false = 0))

## remove any NAs in predictors:::
nas <- mdf_scaled %>% 
    filter(is.na(elev_sc),
           is.na(snailRichness_sc),
           is.na(birdRichness_sc),
           is.na(shellDim_sc),
           is.na(relativeMass_sc),
           is.na(ratio_sc),
           is.na(n_sc),
           is.na(Species),
           is.na(siteID)
           )

mdf_scaled2 <- mdf_scaled %>% 
    filter(!is.na(elev_sc),
           !is.na(snailRichness_sc),
           !is.na(birdRichness_sc),
           !is.na(shellDim_sc),
           !is.na(relativeMass_sc),
           !is.na(ratio_sc),
           !is.na(n_sc)
    )


### set priors
priors_step1 <- c(
    # Intercept: on logit scale, Normal(0, 1.5) keeps most prior probability
    # between ~5% and ~95% probability of monomorphism — weakly informative
    prior(normal(0, 1.5), class = Intercept),
    
    # Fixed effects: scaled predictors, so Normal(0, 1) allows effects up to
    # ~±2 on logit scale (odds ratio ~0.14 to 7.4) — permissive but regularizing
    prior(normal(0, 1), class = b),
    
    # Random effect SDs: half-Normal(0, 1) is the standard weakly informative
    # choice for variance components on logit scale — allows substantial
    # among-species/among-site variation but pulls toward zero
    prior(normal(0, 1), class = sd)
)

## fit model
step1 <- brm(
    formula = bf(
        is_monomorphic ~ 
            bio12_sc + 
            elev_sc + 
            snailRichness_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc +
            n_sc + # account for the number of individuals collected in a lot
            (1 | Species) + (1 | siteID) # site-specific intercept probably makes sense too
        ),
    family = bernoulli(),
    data = mdf_scaled,
    cores = 4,
    prior = priors_step1,
    threads = threading(4),
    control = list(adapt_delta = 0.98),
    backend = "cmdstanr",
    file = "brmsFittedModels/manuscript/step1.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step1)

## mpd for snailRichness
step1_mpd <- brm(
    formula = bf(
        is_monomorphic ~ 
            bio12_sc + 
            elev_sc + 
            mpd_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc +
            n_sc + # account for the number of individuals collected in a lot
            (1 | Species) + (1 | siteID) 
    ),
    family = bernoulli(),
    data = mdf_scaled,
    cores = 4,
    prior = priors_step1,
    threads = threading(4),
    control = list(adapt_delta = 0.98),
    backend = "cmdstanr",
    file = "brmsFittedModels/manuscript/step1_mpd.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step1_mpd)


######## now for step 2 model ####################
mdf_scaled_poly <- filter(mdf_scaled, is_monomorphic == 0)

priors_step2 <- c(
    # Intercept: Normal(0, 1.5) keeps most prior mass between ~5-95%
    # probability of being dextral — same logic as step1, neutral and weakly informative
    prior(normal(0, 1.5), class = Intercept),
    
    # Fixed effects: Normal(0, 1) — same as step1, allows odds ratios up to ~7x
    # per SD change in predictors, which is permissive for ecological data
    prior(normal(0, 1), class = b),
    
    # Random effect SDs: half-Normal(0, 1) — same reasoning as step1
    prior(normal(0, 1), class = sd)
)

step2 <- brm(
    formula = bf(
        dex | trials(n) ~ 
            bio12_sc + 
            elev_sc + 
            snailRichness_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc + 
            (1 | Species) + (1 | siteID) 
    ),
    family = binomial(),
    data = mdf_scaled_poly,
    cores = 4,
    prior = priors_step2,
    threads = threading(4),
    control = list(adapt_delta = 0.98),
    backend = "cmdstanr",
    file = "brmsFittedModels/manuscript/step2.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step2)

### mpd
step2_mpd <- brm(
    formula = bf(
        dex | trials(n) ~ 
            bio12_sc + 
            elev_sc + 
            mpd_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc +
            (1 | Species) + (1 | siteID) 
    ),
    family = binomial(),
    data = mdf_scaled_poly,
    cores = 4,
    prior = priors_step2,
    threads = threading(4),
    control = list(adapt_delta = 0.98),
    backend = "cmdstanr",
    file = "brmsFittedModels/manuscript/step2_mpd.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step2_mpd)
