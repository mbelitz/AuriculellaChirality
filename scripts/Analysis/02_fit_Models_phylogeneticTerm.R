library(ape)
library(tidyverse)
library(brms)

# Read in and prepare data

mdf_scaled <- read.csv("data/mdf_scaled_July2026Update.csv") %>% 
    mutate(n_sc  = scale(n),
           mpd_sc = scale(MPD)) %>% 
    filter(correctedspeciesname != "Auriculella cerea")

# Create monomorphic indicator
mdf_scaled <- mdf_scaled %>% 
    mutate(is_monomorphic = if_else(dex == 0 | dex == n, 1, 0))

# Standardize species names to match tree format
mdf_scaled <- mdf_scaled %>%
    mutate(Species_full = paste0("Auriculella_", Species))

# 2. Prepare phylogeny

tree <- read.tree("data/auriculella_species_collapsed.tre")

# --- Diagnose mismatches ---
data_species_full <- unique(mdf_scaled$Species_full)

not_in_tree <- setdiff(data_species_full, tree$tip.label)
not_in_data  <- setdiff(tree$tip.label, data_species_full)
in_both      <- intersect(data_species_full, tree$tip.label)

cat("Species in data but NOT in tree:", length(not_in_tree), "\n")
print(not_in_tree)
cat("\nSpecies in tree but NOT in data:", length(not_in_data), "\n")
print(not_in_data)
cat("\nSpecies in both:", length(in_both), "\n")
print(in_both)

# --- Prune tree to remove tips not in data (outgroups + unsampled Auriculella) ---
tree_pruned <- drop.tip(tree, not_in_data)

cat("\nPruned tree tips:\n")
print(tree_pruned$tip.label)

# --- Build phylogenetic covariance matrix ---
A_tree <- vcv.phylo(tree_pruned, corr = TRUE)

# --- Expand matrix to include ALL species in data ---
# Species not in tree get 0 covariance with others (phylogenetically uninformative)
# but diagonal = 1 (required for valid correlation matrix)
all_species <- sort(unique(mdf_scaled$Species_full))
n_spp       <- length(all_species)

A_full <- matrix(0,
                 nrow = n_spp, ncol = n_spp,
                 dimnames = list(all_species, all_species))

# Fill in known phylogenetic covariances
A_full[rownames(A_tree), colnames(A_tree)] <- A_tree

# Ensure diagonal is exactly 1 for all species
diag(A_full) <- 1

# --- Validate matrix ---
cat("\nMatrix dimensions:", dim(A_full), "\n")
cat("All species represented:", all(all_species %in% rownames(A_full)), "\n")

# Check positive definiteness (required for brms)
eigenvalues <- eigen(A_full)$values
cat("Min eigenvalue:", min(eigenvalues), "\n")

# If any eigenvalues <= 0, add small jitter to diagonal
if (any(eigenvalues <= 0)) {
    cat("Warning: matrix not positive definite — adding jitter\n")
    A_full <- A_full + diag(1e-6, n_spp)
}

# Define priors

priors_step1 <- c(
    prior(normal(0, 1.5), class = Intercept),
    prior(normal(0, 1),   class = b),
    prior(normal(0, 1),   class = sd)
)

priors_step2 <- c(
    prior(normal(0, 1.5), class = Intercept),
    prior(normal(0, 1),   class = b),
    prior(normal(0, 1),   class = sd)
)

# Prepare data subsets

# Step 2 uses only polymorphic lots
mdf_scaled_poly <- mdf_scaled %>% filter(is_monomorphic == 0)

cat("\nStep 1 data rows:", nrow(mdf_scaled), "\n")
cat("Step 2 data rows (polymorphic only):", nrow(mdf_scaled_poly), "\n")

# Fit models; snail richness versions

# Step 1: probability of monomorphism
step1_phylo <- brm(
    formula = bf(
        is_monomorphic ~ 
            bio12_sc + 
            elev_sc + 
            snailRichness_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc +
            n_sc + 
            (1 | gr(Species_full, cov = A_full)) +
            (1 | siteID)
    ),
    family   = bernoulli(),
    data     = mdf_scaled,
    data2    = list(A_full = A_full),
    prior    = priors_step1,
    cores    = 4,
    threads  = threading(4),
    control  = list(adapt_delta = 0.98),
    backend  = "cmdstanr",
    file     = "brmsFittedModels/manuscript/step1_phylo.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step1_phylo)

# Step 2: proportion dextral among polymorphic lots
step2_phylo <- brm(
    formula = bf(
        dex | trials(n) ~ 
            bio12_sc + 
            elev_sc + 
            snailRichness_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc + 
            (1 | gr(Species_full, cov = A_full)) +
            (1 | siteID)
    ),
    family   = binomial(),
    data     = mdf_scaled_poly,
    data2    = list(A_full = A_full),
    prior    = priors_step2,
    cores    = 4,
    threads  = threading(4),
    control  = list(adapt_delta = 0.98),
    backend  = "cmdstanr",
    file     = "brmsFittedModels/manuscript/step2_phylo.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step2_phylo)

# Fit models; mpd version

# Filter out MPD NAs explicitly
mdf_scaled_mpd      <- mdf_scaled %>% filter(!is.na(MPD))
mdf_scaled_poly_mpd <- mdf_scaled_mpd %>% filter(is_monomorphic == 0)

cat("\nStep 1 MPD data rows:", nrow(mdf_scaled_mpd), "\n")
cat("Step 2 MPD data rows:", nrow(mdf_scaled_poly_mpd), "\n")

step1_phylo_mpd <- brm(
    formula = bf(
        is_monomorphic ~ 
            bio12_sc + 
            elev_sc + 
            mpd_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc + 
            n_sc + 
            (1 | gr(Species_full, cov = A_full)) +
            (1 | siteID)
    ),
    family   = bernoulli(),
    data     = mdf_scaled_mpd,
    data2    = list(A_full = A_full),
    prior    = priors_step1,
    cores    = 4,
    threads  = threading(4),
    control  = list(adapt_delta = 0.98),
    backend  = "cmdstanr",
    file     = "brmsFittedModels/manuscript/step1_phylo_mpd.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change" 
)

summary(step1_phylo_mpd)

step2_phylo_mpd <- brm(
    formula = bf(
        dex | trials(n) ~ 
            bio12_sc + 
            elev_sc + 
            mpd_sc + 
            birdRichness_sc + 
            shellDim_sc + 
            relativeMass_sc + 
            relativeRatio_sc + 
            (1 | gr(Species_full, cov = A_full)) +
            (1 | siteID)
    ),
    family   = binomial(),
    data     = mdf_scaled_poly_mpd,
    data2    = list(A_full = A_full),
    prior    = priors_step2,
    cores    = 4,
    threads  = threading(4),
    control  = list(adapt_delta = 0.98),
    backend  = "cmdstanr",
    file     = "brmsFittedModels/manuscript/step2_phylo_mpd.rds",
    save_pars = save_pars(all = TRUE),       
    file_refit = "on_change"    
)

summary(step2_phylo_mpd)