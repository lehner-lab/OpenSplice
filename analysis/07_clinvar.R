## 07_clinvar.R
## ClinVar benchmarking: bar charts and ROC curves.
## All plots are produced twice — once for all status_simple levels ("all")
## and once restricted to high-confidence annotations (status_simple ≥ 2).

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(scales)
  library(pROC)
  library(tidyr)
  library(here)
})

source(here("analysis", "config.R"))

# ── Paths ─────────────────────────────────────────────────────────────────────
MASTER_FILE <- MASTER_TABLE
PLOT_DIR    <- here("figures", "07_clinvar")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

results_dir <- here("results", "analysis", "07_clinvar")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ── Colour palettes ───────────────────────────────────────────────────────────
effect_pal <- c(
  "Not sig (FDR ≥ 0.1)" = "grey80",
  "Sig ↓ (ΔPSI < 0)"    = "#2171B5",
  "Sig ↑ (ΔPSI > 0)"    = "firebrick"
)

delta_bin_pal <- c(
  pos5 = "#8B0000",
  pos4 = "firebrick",
  pos3 = "#F88379",
  pos2 = "#F4C2C2",
  pos1 = "#FFD1DC",
  mid  = "gray70",
  neg1 = "#E0FFFF",
  neg2 = "#AFDBF5",
  neg3 = "#7AAED4",
  neg4 = "steelblue",
  neg5 = "#00416A"
)

delta_bin_labels <- c( "x > 3","2 < x ≤ 3", "1 < x ≤ 2","0.5 < x ≤ 1", "0 < x ≤ 0.5",
                 "FDR ≥ 0.1",
                 "-0.5 ≤ x < 0","-1 ≤ x < -0.5", "-2 ≤ x < -1","-3 ≤ x < -2", "x < -3" )

# ── Shared ggplot theme ───────────────────────────────────────────────────────
theme_cv <- function() {
  theme_minimal() +
    theme(
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      axis.line        = element_line(color = "black", linewidth = 0.3),
      axis.text        = element_text(size = 10, color = "black", family = "Helvetica"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.title       = element_text(size = 10, color = "black", family = "Helvetica"),
      strip.text       = element_text(size = 10, color = "black", family = "Helvetica"),
      plot.title       = element_text(size = 10, color = "black", family = "Helvetica"),
      plot.background  = element_rect(color = NA),
      panel.grid       = element_blank(),
      plot.margin      = margin(4,4,4,4)
    )
}

# ── Load & prepare base data ──────────────────────────────────────────────────
message("Loading master table...")
master <- fread(MASTER_FILE)

clnsig_levels <- c("Pathogenic", "Likely Pathogenic", "VUS/Conflicting",
                   "Likely Benign", "Benign")

clinvar_filter <- master %>%
  filter(
    !is.na(delta_logit),
    mc_simple %in% c("Synonymous", "Splice site", "Intronic"),
    !clnsig_simple %in% c("Not provided", "Other"),
    wt_psi >= 5,
    !grepl("frameshift|nonsense", mc)
  ) %>%
  select(nt_seq, delta_psi, delta_logit, significant, mc_simple,
         status_simple, clnsig_simple, delta_bin) %>%
  distinct() %>%
  mutate(
    mc_simple     = factor(mc_simple,     levels = c("Splice site", "Intronic", "Synonymous")),
    clnsig_simple = factor(clnsig_simple, levels = clnsig_levels),
    delta_bin     = factor(delta_bin, levels = c("mid", "pos5", "pos4", "pos3", "pos2", "pos1",
                                                  "neg5", "neg4", "neg3", "neg2", "neg1")),
    delta_cat = case_when(
      is.na(significant)          ~ NA_character_,
      significant != "yes"        ~ "Not sig (FDR ≥ 0.1)",
      significant == "yes" & delta_psi < 0 ~ "Sig ↓ (ΔPSI < 0)",
      significant == "yes" & delta_psi > 0 ~ "Sig ↑ (ΔPSI > 0)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(delta_cat))

message("  Total ClinVar variants after base filter: ", nrow(clinvar_filter))

# ── Filter conditions ─────────────────────────────────────────────────────────
filter_levels <- list(
  all    = clinvar_filter,
  high   = clinvar_filter %>% filter(status_simple >= 2)
)
filter_labels <- c(all = "All status", high = "Status ≥ 2")

# ── Helper: count labels for bar charts ───────────────────────────────────────
make_bar_n <- function(df) {
  df %>% count(mc_simple, clnsig_simple, name = "n")
}

# ── Helper: ROC computation ───────────────────────────────────────────────────
compute_roc <- function(df, title_label) {
  df_roc <- df %>%
    filter(clnsig_simple != "VUS/Conflicting") %>%
    mutate(class = case_when(
      clnsig_simple %in% c("Likely Pathogenic", "Pathogenic") ~ 1L,
      clnsig_simple %in% c("Likely Benign",     "Benign")     ~ 0L
    )) %>%
    filter(!is.na(class))

  if (length(unique(df_roc$class)) < 2) {
    message("  Skipping ROC for '", title_label, "': only one class present")
    return(NULL)
  }

  roc_obj <- roc(
    response  = factor(df_roc$class, levels = c(0, 1)),
    predictor = df_roc$delta_logit,
    levels    = c(0, 1),
    direction = ">",
    quiet     = TRUE
  )

  tibble(
    threshold         = roc_obj$thresholds,
    specificity       = roc_obj$specificities,
    sensitivity       = roc_obj$sensitivities,
    `1 - specificity` = 1 - specificity,
    auc_roc           = as.numeric(auc(roc_obj)),
    title             = title_label,
    n_path            = sum(df_roc$class == 1),
    n_benign          = sum(df_roc$class == 0)
  ) %>%
    arrange(`1 - specificity`)
}

# ── Helper: ROC plot ──────────────────────────────────────────────────────────
plot_roc <- function(roc_df) {
  if (is.null(roc_df)) return(NULL)
  auc_label <- paste0("AUC = ", round(roc_df$auc_roc[1], 3))
  ggplot(roc_df, aes(x = `1 - specificity`, y = sensitivity)) +
    geom_line(linewidth = 0.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    annotate("text", x = 1, y = 0.05, label = auc_label, size = 2.8, hjust = 1) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = roc_df$title[1], x = "1 - Specificity", y = "Sensitivity") +
    theme_cv() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "none",
          plot.margin      = margin(10,10,10,10))
}

# ── Generate plots for each filter level ─────────────────────────────────────
for (flt in names(filter_levels)) {

  df   <- filter_levels[[flt]]
  lbl  <- filter_labels[[flt]]
  tag  <- flt  # "all" or "high"

  message("\n── Filter: ", lbl, "  (n = ", nrow(df), ") ──")
  title_suffix <- if (tag == "high") "  [star≥2]" else ""
  bar_n <- make_bar_n(df)

  # ── Plot 1: ΔLogitPSI bins ────────────────────────────────────────────────
  p_bin <- ggplot(df, aes(x = clnsig_simple, fill = delta_bin)) +
    geom_bar(position = "fill", color = NA) +
    geom_text(
      data = bar_n,
      aes(x = clnsig_simple, y = 1, label = n),
      inherit.aes = FALSE, vjust = -0.3, size = 2.8
    ) +
    scale_y_continuous(labels = label_percent(suffix = "")) +
    scale_fill_manual(
      values = delta_bin_pal,
      breaks = names(delta_bin_pal),
      labels = delta_bin_labels,
      name   = expression(Delta * LogitPSI)
    ) +
    facet_wrap(~mc_simple, ncol = 4, scales = "free_y") +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("ΔLogitPSI by ClinVar classification", title_suffix),
      x = "", y = "Percentage of variants"
    ) +
    theme_cv()+
    theme(strip.text = element_text(margin = margin(b = 8)))

  ggsave(file.path(PLOT_DIR, paste0("clinvar_delta_bin_", tag, ".png")),
         p_bin, width = 10, height = 5, dpi = 300)

  # ── Plot 2: Simplified directionality ────────────────────────────────────
  p_cat <- ggplot(df, aes(x = clnsig_simple, fill = delta_cat, color = delta_cat)) +
    geom_bar(position = "fill", alpha = 0.8) +
    geom_text(
      data = bar_n,
      aes(x = clnsig_simple, y = 1, label = n),
      inherit.aes = FALSE, vjust = -0.3, size = 3
    ) +
    scale_y_continuous(labels = label_percent(suffix = "")) +
    scale_fill_manual(values  = effect_pal) +
    scale_color_manual(values = effect_pal) +
    facet_wrap(~mc_simple, ncol = 4, scales = "free_y") +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0("Splicing effect by ClinVar classification", title_suffix),
      x = "", y = "Percentage of variants", fill = "", color = ""
    ) +
    theme_cv()+
    theme(strip.text = element_text(margin = margin(b = 8)))

  ggsave(file.path(PLOT_DIR, paste0("clinvar_delta_cat_", tag, ".png")),
         p_cat, width = 10, height = 5, dpi = 300)

  # ── ROC: intronic only ───────────────────────────────────────────────────
  roc_intron <- compute_roc(
    df %>% filter(mc_simple == "Intronic"),
    paste0("Intronic", title_suffix)
  )
  p_roc_intron <- plot_roc(roc_intron)
  if (!is.null(p_roc_intron)) {
    ggsave(file.path(PLOT_DIR, paste0("roc_intronic_", tag, ".png")),
           p_roc_intron, width = 4, height = 4, dpi = 300)
  }

  # ── ROC: all variant classes ─────────────────────────────────────────────
  roc_all <- compute_roc(df, paste0("All variants", title_suffix))
  p_roc_all <- plot_roc(roc_all)
  if (!is.null(p_roc_all)) {
    ggsave(file.path(PLOT_DIR, paste0("roc_all_", tag, ".png")),
           p_roc_all, width = 4, height = 4, dpi = 300)
  }

  # ── Summary table: effect categories ─────────────────────────────────────
  filter_status <- ifelse(flt =='all',0,1)

  res <- df %>%
    filter(status_simple > filter_status) %>%
    mutate(
      clnsig_simple = case_when(
        clnsig_simple == "Likely Pathogenic" ~ "Pathogenic",
        clnsig_simple == "Likely Benign"     ~ "Benign",
        TRUE ~ as.character(clnsig_simple)
      ),
      effect = case_when(
        significant == "no"    ~ "neutral",
        delta_logit <= -1      ~ "decrease_strong",
        delta_logit <   0      ~ "decrease",
        delta_logit >   0      ~ "increase",
        TRUE                   ~ "neutral"
      )
    ) %>%
    group_by(mc_simple, clnsig_simple, effect) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(mc_simple, clnsig_simple) %>%
    mutate(total = sum(n), percent = round(100 * n / total, 1)) %>%
    ungroup() %>%
    pivot_wider(names_from = effect, values_from = c("percent", "n"), values_fill = 0)

  fwrite(res, file.path(results_dir, paste0("summary_effect_", tag, ".tsv")), sep = "\t")
  message("  Summary table written (", tag, "): ", nrow(res), " rows")
}

# ── Per-mc_simple variant counts ─────────────────────────────────────────────
message("\nVariant counts per molecular consequence (base filter, no status filter):")
print(table(clinvar_filter$mc_simple))
message("\nVariant counts per molecular consequence (status ≥ 2):")
print(table(filter_levels$high$mc_simple))

message("\n07_clinvar.R complete. Plots saved to: ", PLOT_DIR)
