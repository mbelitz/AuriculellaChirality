library(tidyverse)

# =============================================================================
# READ AND PREPARE DATA — matching exactly what 01_fit_Models_v3.R does
# =============================================================================

mdf_scaled <- read.csv("data/mdf_scaled_July2026Update.csv") %>%
    mutate(
        n_sc   = scale(n),
        mpd_sc = scale(MPD)
    ) %>%
    filter(correctedspeciesname != "Auriculella cerea")

mdf_scaled <- mdf_scaled %>%
    mutate(is_monomorphic = if_else(dex == 0 | dex == n, 1L, 0L))

# Polymorphic subset (Stage 2)
mdf_scaled_poly <- mdf_scaled %>% filter(is_monomorphic == 0)

# MPD subsets (exclude NAs)
mdf_scaled_mpd      <- mdf_scaled %>% filter(!is.na(MPD))
mdf_scaled_poly_mpd <- mdf_scaled_mpd %>% filter(is_monomorphic == 0)

# =============================================================================
# VERIFICATION
# =============================================================================

cat("================================================================\n")
cat("STAGE 1 — Species richness variant (M1a / M1a-P)\n")
cat("================================================================\n")
cat("Collection events (n rows)  :", nrow(mdf_scaled), "\n")
cat("Total individuals (sum of n):", sum(mdf_scaled$n, na.rm = TRUE), "\n")
cat("Species                     :", n_distinct(mdf_scaled$Species), "\n")
cat("Sites                       :", n_distinct(mdf_scaled$siteID), "\n")

cat("\n================================================================\n")
cat("STAGE 2 — Species richness variant (M2a / M2a-P)\n")
cat("================================================================\n")
cat("Collection events (n rows)  :", nrow(mdf_scaled_poly), "\n")
cat("Total individuals (sum of n):", sum(mdf_scaled_poly$n, na.rm = TRUE), "\n")
cat("Species                     :", n_distinct(mdf_scaled_poly$Species), "\n")
cat("Sites                       :", n_distinct(mdf_scaled_poly$siteID), "\n")

cat("\n================================================================\n")
cat("STAGE 1 — MPD variant (M1b / M1b-P)\n")
cat("================================================================\n")
cat("Collection events (n rows)  :", nrow(mdf_scaled_mpd), "\n")
cat("Total individuals (sum of n):", sum(mdf_scaled_mpd$n, na.rm = TRUE), "\n")
cat("Species                     :", n_distinct(mdf_scaled_mpd$Species), "\n")
cat("Sites                       :", n_distinct(mdf_scaled_mpd$siteID), "\n")

cat("\n================================================================\n")
cat("STAGE 2 — MPD variant (M2b / M2b-P)\n")
cat("================================================================\n")
cat("Collection events (n rows)  :", nrow(mdf_scaled_poly_mpd), "\n")
cat("Total individuals (sum of n):", sum(mdf_scaled_poly_mpd$n, na.rm = TRUE), "\n")
cat("Species                     :", n_distinct(mdf_scaled_poly_mpd$Species), "\n")
cat("Sites                       :", n_distinct(mdf_scaled_poly_mpd$siteID), "\n")

cat("\n================================================================\n")
cat("EXPECTED (from methods text)\n")
cat("================================================================\n")
cat("Stage 1 SR  : 3,863 events | 88,403 individuals | 18 species | 1,990 sites\n")
cat("Stage 2 SR  : 1,546 events | 39,022 individuals | 18 species |   712 sites\n")
cat("Stage 1 MPD : 3,250 events\n")
cat("Stage 2 MPD : 1,311 events\n")

cat("\n================================================================\n")
cat("MATCH CHECK\n")
cat("================================================================\n")
cat("Stage 1 SR  events   match:", nrow(mdf_scaled)          == 3863, "\n")
cat("Stage 1 SR  indivs   match:", sum(mdf_scaled$n)         == 88403, "\n")
cat("Stage 1 SR  species  match:", n_distinct(mdf_scaled$Species) == 18, "\n")
cat("Stage 1 SR  sites    match:", n_distinct(mdf_scaled$siteID)  == 1990, "\n")
cat("Stage 2 SR  events   match:", nrow(mdf_scaled_poly)     == 1546, "\n")
cat("Stage 2 SR  indivs   match:", sum(mdf_scaled_poly$n)    == 39022, "\n")
cat("Stage 2 SR  sites    match:", n_distinct(mdf_scaled_poly$siteID) == 712, "\n")
cat("Stage 1 MPD events   match:", nrow(mdf_scaled_mpd)      == 3250, "\n")
cat("Stage 2 MPD events   match:", nrow(mdf_scaled_poly_mpd) == 1311, "\n")