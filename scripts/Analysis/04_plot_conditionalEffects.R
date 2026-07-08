###############################################################################
# 04_plot_conditionalEffects.R
#
# Publication-quality conditional-effects figures + coefficient tables for the
# two-step hurdle-style brms models (step1 = P(monomorphic), step2 = P(dextral)
# among polymorphic lots). Handles all 8 fitted models:
#
#   step1.rds, step2.rds                         (snail richness, non-phylo)
#   step1_mpd.rds, step2_mpd.rds                 (MPD, non-phylo)
#   step1_phylo.rds, step2_phylo.rds             (snail richness, phylogenetic)
#   step1_phylo_mpd.rds, step2_phylo_mpd.rds     (MPD, phylogenetic)
#
# Only predictors whose 95% posterior CI excludes zero are plotted. A
# coefficient table (all fixed effects) is written per model plus a single
# combined table across all models.
#
# Run from the project root (i.e. the directory containing brmsFittedModels/).
###############################################################################

library(tidyverse)
library(brms)

# ----- setup ----------------------------------------------------------------

dir.create("figures/manuscript", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/manuscript",  recursive = TRUE, showWarnings = FALSE)

# Pretty labels for predictors. Update `ratio_sc` / `n_sc` if the variable
# names mean something different in your manuscript.
label_map <- c(
    bio12_sc         = "Annual precipitation",
    elev_sc          = "Elevation",
    snailRichness_sc = "Snail richness",
    mpd_sc           = "Snail MPD",
    birdRichness_sc  = "Bird richness",
    shellDim_sc      = "Shell size",
    relativeMass_sc  = "Relative mass",
    ratio_sc         = "Shell ratio",
    n_sc             = "Collection sample size"
)

# Control variables: never plotted in conditional-effects figures even when
# their CIs exclude zero. They still appear in the coefficient tables.
control_vars <- c("bio12_sc", "elev_sc", "n_sc")

# ----- helpers --------------------------------------------------------------

# Return fixed-effect predictor names whose 95% CI excludes zero (dropping
# the intercept and any control variables). Uses fixef() so it works for any
# brmsfit family.
get_sig_terms <- function(fit, prob = 0.95) {
    tail <- (1 - prob) / 2
    fe <- fixef(fit, probs = c(tail, 1 - tail))
    fe <- as.data.frame(fe)
    fe$term <- rownames(fe)
    lo <- fe[, 3]
    hi <- fe[, 4]
    sig <- (lo > 0 & hi > 0) | (lo < 0 & hi < 0)
    setdiff(fe$term[sig], c("Intercept", control_vars))
}

# Build a tidy conditional_effects data frame for a set of predictors.
# For binomial-with-trials models we condition on n = 1 so estimate__ is
# interpretable as a per-individual probability.
make_ce_df <- function(fit, vars, prob = 0.95) {
    is_binom <- isTRUE(fit$family$family == "binomial")
    ce <- conditional_effects(
        fit,
        effects    = vars,
        re_formula = NA,
        prob       = prob,
        conditions = if (is_binom) data.frame(n = 1) else NULL
    )
    ce_df <- bind_rows(
        lapply(names(ce), function(v) {
            as_tibble(ce[[v]]) %>% mutate(predictor = v)
        })
    )
    ce_df$predictor_lab <- factor(
        ce_df$predictor,
        levels = vars,
        labels = label_map[vars]
    )
    ce_df
}

# Panel-letter data frame in the same order as `vars`.
make_panel_df <- function(vars) {
    tibble(
        predictor     = vars,
        predictor_lab = factor(label_map[vars], levels = label_map[vars]),
        panel         = LETTERS[seq_along(vars)]
    )
}

# Core plotting function. Auto-sizes the panel grid by number of significant
# predictors, but caller can override.
plot_ce <- function(fit, out_png, y_label,
                    prob = 0.95,
                    ncol_wrap = 2,
                    width = NULL, height = NULL) {
    vars <- get_sig_terms(fit, prob = prob)
    if (length(vars) == 0) {
        message("No significant predictors for ", out_png, " — skipping figure.")
        return(invisible(NULL))
    }
    ce_df    <- make_ce_df(fit, vars, prob = prob)
    panel_df <- make_panel_df(vars)
    
    # adaptive grid: >4 significant predictors -> 2 cols (3x2 portrait layout),
    # otherwise up to 3 across in a single row.
    if (is.null(ncol_wrap)) {
        ncol_wrap <- if (length(vars) > 4) 2 else min(length(vars), 3)
    }
    nrow_wrap <- ceiling(length(vars) / ncol_wrap)
    if (is.null(width))  width  <- 3.0 * ncol_wrap + 0.5
    if (is.null(height)) height <- 2.6 * nrow_wrap + 0.8
    
    p <- ggplot(ce_df, aes(x = effect1__, y = estimate__)) +
        geom_ribbon(
            aes(ymin = lower__, ymax = upper__),
            alpha = 0.25
        ) +
        geom_line(linewidth = 1) +
        geom_text(
            data = panel_df,
            aes(x = -Inf, y = Inf, label = panel),
            hjust = -0.5,
            vjust = 2.5,
            inherit.aes = FALSE,
            fontface = "bold",
            size = 4
        ) +
        facet_wrap(
            ~ predictor_lab,
            scales = "free_x",
            ncol = ncol_wrap,
            strip.position = "bottom"
        ) +
        labs(x = NULL, y = y_label) +
        theme_classic(base_size = 12) +
        theme(
            strip.background = element_blank(),
            strip.placement  = "outside",
            strip.text       = element_text(face = "bold"),
            panel.spacing    = unit(1.2, "lines"),
            axis.title.y     = element_text(size = 12)
        )
    ggsave(out_png, p, width = width, height = height, dpi = 450)
    message("Wrote ", out_png, " (", length(vars), " predictors: ",
            paste(vars, collapse = ", "), ")")
    invisible(p)
}

# Coefficient table: all fixed effects with 95% CI, flagged for significance.
make_coef_table <- function(fit, out_csv, model_name, prob = 0.95) {
    tail <- (1 - prob) / 2
    fe <- fixef(fit, probs = c(tail, 1 - tail))
    lo_col <- paste0("Q", formatC(tail * 100, format = "f", digits = 1))
    hi_col <- paste0("Q", formatC((1 - tail) * 100, format = "f", digits = 1))
    
    df <- as.data.frame(fe) %>%
        rownames_to_column("Term") %>%
        rename(
            Lower = all_of(lo_col),
            Upper = all_of(hi_col)
        ) %>%
        mutate(
            Term_label = ifelse(Term %in% names(label_map),
                                label_map[Term], Term),
            Significant = (Lower > 0 & Upper > 0) | (Lower < 0 & Upper < 0),
            Model = model_name
        ) %>%
        select(Model, Term, Term_label,
               Estimate, Est.Error, Lower, Upper, Significant)
    
    write_csv(df, out_csv)
    message("Wrote ", out_csv)
    invisible(df)
}

# ----- model inventory -------------------------------------------------------
#
# For every model: rds path, y-axis label for conditional effects, file stems
# for figure + table, and a short display name for the combined table.

models <- list(
    list(
        file = "brmsFittedModels/manuscript/step1.rds",
        fig  = "figures/manuscript/step1_ce.png",
        tab  = "tables/manuscript/step1_coef.csv",
        ylab = "Probability of monomorphic lot",
        name = "Step 1 (snail richness)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step2.rds",
        fig  = "figures/manuscript/step2_ce.png",
        tab  = "tables/manuscript/step2_coef.csv",
        ylab = "Probability of dextral individual\n(polymorphic lots)",
        name = "Step 2 (snail richness)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step1_mpd.rds",
        fig  = "figures/manuscript/step1_mpd_ce.png",
        tab  = "tables/manuscript/step1_mpd_coef.csv",
        ylab = "Probability of monomorphic lot",
        name = "Step 1 (MPD)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step2_mpd.rds",
        fig  = "figures/manuscript/step2_mpd_ce.png",
        tab  = "tables/manuscript/step2_mpd_coef.csv",
        ylab = "Probability of dextral individual\n(polymorphic lots)",
        name = "Step 2 (MPD)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step1_phylo.rds",
        fig  = "figures/manuscript/step1_phylo_ce.png",
        tab  = "tables/manuscript/step1_phylo_coef.csv",
        ylab = "Probability of monomorphic lot",
        name = "Step 1 phylogenetic (snail richness)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step2_phylo.rds",
        fig  = "figures/manuscript/step2_phylo_ce.png",
        tab  = "tables/manuscript/step2_phylo_coef.csv",
        ylab = "Probability of dextral individual\n(polymorphic lots)",
        name = "Step 2 phylogenetic (snail richness)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step1_phylo_mpd.rds",
        fig  = "figures/manuscript/step1_phylo_mpd_ce.png",
        tab  = "tables/manuscript/step1_phylo_mpd_coef.csv",
        ylab = "Probability of monomorphic lot",
        name = "Step 1 phylogenetic (MPD)"
    ),
    list(
        file = "brmsFittedModels/manuscript/step2_phylo_mpd.rds",
        fig  = "figures/manuscript/step2_phylo_mpd_ce.png",
        tab  = "tables/manuscript/step2_phylo_mpd_coef.csv",
        ylab = "Probability of dextral individual\n(polymorphic lots)",
        name = "Step 2 phylogenetic (MPD)"
    )
)

# ----- run everything --------------------------------------------------------

all_coefs <- list()

for (m in models) {
    if (!file.exists(m$file)) {
        warning("Missing model file: ", m$file, " — skipping.")
        next
    }
    message("\n==== ", m$name, " ====")
    fit <- readRDS(m$file)
    
    plot_ce(fit, m$fig, m$ylab, prob = 0.95)
    all_coefs[[m$name]] <- make_coef_table(fit, m$tab, m$name, prob = 0.95)
}

# Combined wide-format-ish table across all models.
if (length(all_coefs) > 0) {
    bind_rows(all_coefs) %>%
        write_csv("tables/manuscript/all_models_coefficients.csv")
    message("\nWrote tables/manuscript/all_models_coefficients.csv")
}
