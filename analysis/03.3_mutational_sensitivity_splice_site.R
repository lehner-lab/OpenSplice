## 03.3_mutational_sensitivity_splice_site.R
##
## Sections:
##   1. 5' SS — ∆LogitPSI vs ∆MaxEntScan per position and overall
##   2. 3' SS — ∆LogitPSI vs ∆MaxEntScan per position and overall
##   3. SS strength vs VastDB average PSI (genome-wide z-score scale)
##   4. SS strength vs WT PSI (mutagenesis libraries)
##   5. SS strength vs overall splicing sensitivity (median |∆LogitPSI|)
##   6. Save per-exon SS info table (EXON_SS_INFO_FILE, used by downstream scripts)
##
## Inputs:
##   MASTER_TABLE, MAXENTSCAN_SS5, MAXENTSCAN_SS3,
##   MAXENTSCAN_GENOME_SS5, MAXENTSCAN_GENOME_SS3,
##   VASTDB_SS_FILE  (gitignored — run 01_wt_screening_libraries pipeline step 03),
##   METADATA_FILE, VASTDB_FILE
##
## Outputs:
##   figures/03_mutational_sensitivity/03.3_mutational_sensitivity_splice_site/
##     ss5_delta_logit_by_position.png
##     ss5_delta_logit_overall.png
##     ss3_delta_logit_by_position.png
##     ss3_delta_logit_overall.png
##     ss3_vs_vastdb_psi.png
##     ss5_vs_vastdb_psi.png
##     zmean_vs_vastdb_psi.png
##     ss3_vs_wt_psi.png
##     ss5_vs_wt_psi.png
##     zmean_vs_wt_psi.png
##     ss3_vs_sensitivity.png
##     ss5_vs_sensitivity.png
##     zmean_vs_sensitivity.png
##   results/analysis/03_mutational_sensitivity/03.3_splice_site/
##     cor_ss5_by_position.tsv
##     cor_ss5_overall.tsv
##     cor_ss3_by_position.tsv
##     cor_ss3_overall.tsv
##     cor_vastdb_psi_vs_ss.tsv
##     cor_wt_psi_vs_ss.tsv
##     cor_sensitivity_vs_ss.tsv
##   data/databases/exon_ss_info.tsv  (EXON_SS_INFO_FILE — committed)

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

source(here("analysis", "config.R"))

plot_dir    <- here("figures", "03_mutational_sensitivity",
                    "03.3_mutational_sensitivity_splice_site")
results_dir <- here("results", "analysis", "03_mutational_sensitivity", "03.3_splice_site")

dir.create(plot_dir,    showWarnings = FALSE, recursive = TRUE)
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)


# ── Shared theme ──────────────────────────────────────────────────────────────
theme_ss <- theme_minimal() +
  theme(
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    axis.line       = element_line(color = "black", linewidth = 0.5),
    panel.grid      = element_blank(),
    panel.border    = element_blank(),
    axis.text       = element_text(size = 10, color = "black", family = 'Helvetica'),
    axis.title      = element_text(size = 10, color = "black", family = 'Helvetica'),
    plot.background = element_rect(color = NA),
    legend.position = "none",
    plot.margin     = margin(6, 6, 6, 6)
  )


# ── Helper: tidy correlation result into a data frame row ─────────────────────
tidy_cor <- function(x, y, label = "") {
  ct <- cor.test(x, y, use = "complete.obs")
  data.frame(
    comparison = label,
    r          = round(unname(ct$estimate), 2),
    p_value    = ct$p.value,
    n          = sum(complete.cases(x, y))
  )
}


# ── Helper functions: map MaxEntScan scores to genome-wide z-scores ───────────
# Sections 3 and 4 each call one variant (differ only in the exon_id column name).

map_mes_to_genome <- function(your_mes, genome_mes5, genome_mes3) {
  ecdf5 <- ecdf(genome_mes5$ss_5_strength)
  mu5   <- mean(genome_mes5$ss_5_strength); sd5 <- sd(genome_mes5$ss_5_strength)
  ecdf3 <- ecdf(genome_mes3$ss_3_strength)
  mu3   <- mean(genome_mes3$ss_3_strength); sd3 <- sd(genome_mes3$ss_3_strength)
  tibble::tibble(
    exon_id = your_mes$COORD,
    MES5    = your_mes$ss5_strength,
    MES3    = your_mes$ss3_strength,
    p5      = ecdf5(MES5), p3 = ecdf3(MES3),
    z5      = (MES5 - mu5) / sd5, z3 = (MES3 - mu3) / sd3,
    p_geom  = sqrt(p5 * p3),
    z_min   = pmin(z5, z3),
    z_mean  = (z5 + z3) / 2
  )
}

map_mes_to_mut <- function(your_mes, genome_mes5, genome_mes3) {
  ecdf5 <- ecdf(genome_mes5$ss_5_strength)
  mu5   <- mean(genome_mes5$ss_5_strength); sd5 <- sd(genome_mes5$ss_5_strength)
  ecdf3 <- ecdf(genome_mes3$ss_3_strength)
  mu3   <- mean(genome_mes3$ss_3_strength); sd3 <- sd(genome_mes3$ss_3_strength)
  tibble::tibble(
    exon_id = your_mes$exon_id,
    MES5    = your_mes$ss5_strength,
    MES3    = your_mes$ss3_strength,
    p5      = ecdf5(MES5), p3 = ecdf3(MES3),
    z5      = (MES5 - mu5) / sd5, z3 = (MES3 - mu3) / sd3,
    p_geom  = sqrt(p5 * p3),
    z_min   = pmin(z5, z3),
    z_mean  = (z5 + z3) / 2
  )
}


# ── 0. Shared data loading ────────────────────────────────────────────────────

# Raw master table (unfiltered — used for wt_df and per-exon sensitivity)
df_raw <- fread(MASTER_TABLE, sep = '\t')

# Substitutions only, with valid ∆PSI — used in Sections 1 & 2
df_subs <- df_raw %>% filter(!is.na(delta_psi) & !grepl('∆', mut))

# MaxEntScan per-sequence scores
maxentscan_ss5 <- fread(MAXENTSCAN_SS5, sep = '\t')
colnames(maxentscan_ss5) <- c('ss5', 'MAXENT')
maxentscan_ss5 <- maxentscan_ss5 %>% unique()

maxentscan_ss3 <- fread(MAXENTSCAN_SS3, sep = '\t')
colnames(maxentscan_ss3) <- c('ss3', 'MAXENT')
maxentscan_ss3 <- maxentscan_ss3 %>% unique()

# Genome-wide MaxEntScan distributions (for z-score normalisation)
genome_mes5 <- fread(MAXENTSCAN_GENOME_SS5, sep = '\t')
colnames(genome_mes5) <- c('ss5', 'ss_5_strength')

genome_mes3 <- fread(MAXENTSCAN_GENOME_SS3, sep = '\t')
colnames(genome_mes3) <- c('ss3', 'ss_3_strength')

# WT sequences with MaxEntScan scores per exon (shared by Sections 1–5)
wt_df <- df_raw %>%
  filter(grepl('wt', variant_id)) %>%
  select(exon_id, variant_id, nt_seq, psi) %>%
  unique() %>%
  mutate(
    ss5   = substr(nt_seq, nchar(nt_seq) - 27, nchar(nt_seq) - 19),
    ss5_u = gsub('T', 'U', ss5),
    ss3   = substr(nt_seq, 51, 73),
    ss5_u = gsub('T', 'U', ss3)
  ) %>%
  merge(maxentscan_ss5, by = 'ss5') %>%
  rename(ss5_strength = MAXENT) %>%
  merge(maxentscan_ss3, by = 'ss3') %>%
  rename(ss3_strength = MAXENT)

wt_df <- as.data.table(wt_df)
setkey(wt_df, exon_id)


# ── 1. 5' SS — ∆LogitPSI vs ∆MaxEntScan ─────────────────────────────────────

df_ss5 <- df_subs %>%
  filter(!grepl('wt', variant_id) &
           start >= nchar(nt_seq) - 27 & start <= nchar(nt_seq) - 19) %>%
  unique() %>%
  mutate(
    ss5   = substr(nt_seq, nchar(nt_seq) - 27, nchar(nt_seq) - 19),
    ss5_u = gsub('T', 'U', ss5)
  ) %>%
  merge(maxentscan_ss5, by = 'ss5')

df_ss5$MAXENT_wt    <- wt_df[.(df_ss5$exon_id)]$ss5_strength
df_ss5$MAXENT_delta <- df_ss5$MAXENT - df_ss5$MAXENT_wt

df_ss5 <- df_ss5 %>%
  rowwise() %>%
  mutate(ss5_pos = start - (nchar(nt_seq) - 27) + 1)  # relative position within 5'SS window
df_ss5$ss5_pos <- factor(df_ss5$ss5_pos, levels = 1:9)

# Correlation per position
cor_ss5_by_pos <- df_ss5 %>%
  group_by(ss5_pos) %>%
  summarise(
    r       = round(cor(delta_logit, MAXENT_delta, use = "complete.obs"), 2),
    n       = sum(complete.cases(delta_logit, MAXENT_delta)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, "\nn = ", n), x = Inf, y = -Inf)

fwrite(cor_ss5_by_pos %>% select(ss5_pos, r, n),
       file.path(results_dir, "cor_ss5_by_position.tsv"), sep = '\t')

p_ss5_delta_by_pos <- ggplot(df_ss5, aes(x = delta_logit, y = MAXENT_delta)) +
  geom_point(alpha = 0.1, size = 0.1) +
  geom_text(data = cor_ss5_by_pos, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3, color = "black") +
  facet_grid(cols = vars(ss5_pos)) +
  theme_ss +
  ylab("∆MaxEntScan score") + xlab("∆LogitPSI")

ggsave(file.path(plot_dir, "ss5_delta_logit_by_position.png"),
       plot = p_ss5_delta_by_pos, height = 3, width = 10)

# Correlation overall
cor_ss5_overall <- df_ss5 %>%
  mutate(ss5_pos2 = 'All') %>%
  group_by(ss5_pos2) %>%
  summarise(
    r       = round(cor(delta_logit, MAXENT_delta, use = "complete.obs"), 2),
    n       = sum(complete.cases(delta_logit, MAXENT_delta)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, "\nn = ", n), x = Inf, y = -Inf)

fwrite(cor_ss5_overall %>% select(r, n),
       file.path(results_dir, "cor_ss5_overall.tsv"), sep = '\t')

p_ss5_delta_overall <- ggplot(df_ss5, aes(x = delta_logit, y = MAXENT_delta)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_text(data = cor_ss5_overall, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3.5, color = "black") +
  theme_ss +
  ylab("∆MaxEntScan score") + xlab("∆LogitPSI")

ggsave(file.path(plot_dir, "ss5_delta_logit_overall.png"),
       plot = p_ss5_delta_overall, height = 3, width = 3)


# ── 2. 3' SS — ∆LogitPSI vs ∆MaxEntScan ─────────────────────────────────────

df_ss3 <- df_subs %>%
  filter(!grepl('wt', variant_id) & !is.na(delta_logit) &
           mut_type == 'sub' & start >= 51 & start <= 73) %>%
  unique() %>%
  mutate(ss3 = substr(nt_seq, 51, 73)) %>%
  merge(maxentscan_ss3, by = 'ss3')

df_ss3$MAXENT_wt    <- wt_df[.(df_ss3$exon_id)]$ss3_strength
df_ss3$MAXENT_delta <- df_ss3$MAXENT - df_ss3$MAXENT_wt

df_ss3 <- df_ss3 %>%
  rowwise() %>%
  mutate(ss3_pos = start - 50)  # relative position within 3'SS window
df_ss3$ss3_pos <- factor(df_ss3$ss3_pos, levels = 1:23)

# Correlation per position
cor_ss3_by_pos <- df_ss3 %>%
  group_by(ss3_pos) %>%
  summarise(
    r       = round(cor(delta_logit, MAXENT_delta, use = "complete.obs"), 2),
    n       = sum(complete.cases(delta_logit, MAXENT_delta)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, "\nn = ", n), x = Inf, y = -Inf)

fwrite(cor_ss3_by_pos %>% select(ss3_pos, r, n),
       file.path(results_dir, "cor_ss3_by_position.tsv"), sep = '\t')

p_ss3_delta_by_pos <- ggplot(df_ss3, aes(x = delta_logit, y = MAXENT_delta)) +
  geom_point(alpha = 0.1, size = 0.1) +
  geom_text(data = cor_ss3_by_pos, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3, color = "black") +
  facet_wrap(~ss3_pos, ncol = 13) +
  theme_ss +
  ylab("∆MaxEntScan score") + xlab("∆LogitPSI")

ggsave(file.path(plot_dir, "ss3_delta_logit_by_position.png"),
       plot = p_ss3_delta_by_pos, height = 5, width = 18)

# Correlation overall
cor_ss3_overall <- df_ss3 %>%
  mutate(ss3_pos2 = 'All') %>%
  group_by(ss3_pos2) %>%
  summarise(
    r       = round(cor(delta_logit, MAXENT_delta, use = "complete.obs"), 2),
    n       = sum(complete.cases(delta_logit, MAXENT_delta)),
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", r, "\nn = ", n), x = Inf, y = -Inf)

fwrite(cor_ss3_overall %>% select(r, n),
       file.path(results_dir, "cor_ss3_overall.tsv"), sep = '\t')

p_ss3_delta_overall <- ggplot(df_ss3, aes(x = delta_logit, y = MAXENT_delta)) +
  geom_point(alpha = 0.1, size = 0.2) +
  geom_text(data = cor_ss3_overall, aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 1.1, vjust = -0.5, size = 3, color = "black") +
  theme_ss +
  ylab("∆MaxEntScan score") + xlab("∆LogitPSI")

ggsave(file.path(plot_dir, "ss3_delta_logit_overall.png"),
       plot = p_ss3_delta_overall, height = 3, width = 3)


# ── 3. SS strength vs VastDB average PSI ──────────────────────────────────────
# Requires VASTDB_SS_FILE (gitignored — run 01_wt_screening_libraries step 03).

df_vastdb <- fread(VASTDB_SS_FILE, sep = '\t') %>%
  select(COORD, ss3_strength, ss5_strength, Average, Min, Max, Range) %>%
  filter(!is.na(Average)) %>%
  unique()

ss_scaled <- map_mes_to_genome(df_vastdb, genome_mes5, genome_mes3)
ss_scaled  <- merge(ss_scaled, df_vastdb, by.x = 'exon_id', by.y = 'COORD')

cor_vastdb <- bind_rows(
  tidy_cor(ss_scaled$Average, ss_scaled$ss3_strength, "3'SS strength vs VastDB PSI"),
  tidy_cor(ss_scaled$Average, ss_scaled$ss5_strength, "5'SS strength vs VastDB PSI"),
  tidy_cor(ss_scaled$Average, ss_scaled$z_mean,       "mean SS z-score vs VastDB PSI")
)
fwrite(cor_vastdb, file.path(results_dir, "cor_vastdb_psi_vs_ss.tsv"), sep = '\t')

p_ss3_vs_vastdb <- ggplot(df_vastdb, aes(x = Average, y = ss3_strength)) +
  geom_point(size = 0.1, alpha = 0.1) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_vastdb$r[1], 2), ", p < 2.2e-16"),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  geom_smooth(method = 'lm', linewidth = 0.7, linetype = 'dashed') +
  theme_ss +
  labs(x = "Average PSI VastDB", y = "3' SS strength")

ggsave(file.path(plot_dir, "ss3_vs_vastdb_psi.png"),
       plot = p_ss3_vs_vastdb, height = 3, width = 3)

p_ss5_vs_vastdb <- ggplot(df_vastdb, aes(x = Average, y = ss5_strength)) +
  geom_point(size = 0.1, alpha = 0.1) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_vastdb$r[2], 2), ", p < 2.2e-16"),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  geom_smooth(method = 'lm', linewidth = 0.7, linetype = 'dashed') +
  theme_ss +
  labs(x = "Average PSI VastDB", y = "5' SS strength")

ggsave(file.path(plot_dir, "ss5_vs_vastdb_psi.png"),
       plot = p_ss5_vs_vastdb, height = 3, width = 3)

p_zmean_vs_vastdb <- ggplot(ss_scaled, aes(x = Average, y = z_mean)) +
  geom_point(size = 0.1, alpha = 0.1) +
  geom_smooth(method = 'lm', linewidth = 0.7, linetype = 'dashed') +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_vastdb$r[3], 2), ", p < 2.2e-16"),
           hjust = -0.05, vjust = -0.5, size = 3.4, family = 'Helvetica') +
  theme_ss +
  labs(x = "Average PSI VastDB", y = "Mean SS strength")

ggsave(file.path(plot_dir, "zmean_vs_vastdb_psi.png"),
       plot = p_zmean_vs_vastdb, height = 3, width = 3)


# ── 4 & 5. SS strength vs WT PSI and sensitivity (mutagenesis) ────────────────

# Per-exon summary: median |∆LogitPSI|, WT PSI, and SS scores
psi <- df_raw %>%
  filter(!grepl('wt', variant_id) & !is.na(delta_logit)) %>%
  select(exon_id, delta_logit, wt_psi, nt_seq) %>%
  distinct() %>%
  group_by(exon_id, wt_psi) %>%
  summarise(median_delta_logit = median(abs(delta_logit)), .groups = "drop")

psi <- merge(psi, wt_df, by = 'exon_id')

ss_scaled_mut <- map_mes_to_mut(psi, genome_mes5, genome_mes3)
psi           <- merge(psi, ss_scaled_mut, by = 'exon_id')

# ── 4. SS strength vs WT PSI ──────────────────────────────────────────────────

cor_wt_psi <- bind_rows(
  tidy_cor(psi$wt_psi, psi$ss3_strength, "3'SS strength vs WT PSI"),
  tidy_cor(psi$wt_psi, psi$ss5_strength, "5'SS strength vs WT PSI"),
  tidy_cor(psi$wt_psi, psi$z_mean,       "mean SS z-score vs WT PSI")
)
fwrite(cor_wt_psi, file.path(results_dir, "cor_wt_psi_vs_ss.tsv"), sep = '\t')

p_ss3_vs_wt_psi <- ggplot(psi, aes(x = wt_psi, y = ss3_strength)) +
  geom_point(size = 0.5, alpha = 0.5) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_wt_psi$r[1], 2),
                          ", p = ", signif(cor_wt_psi$p_value[1], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "WT PSI — Mutagenesis", y = "3' SS strength")

ggsave(file.path(plot_dir, "ss3_vs_wt_psi.png"),
       plot = p_ss3_vs_wt_psi, height = 3, width = 3)

p_ss5_vs_wt_psi <- ggplot(psi, aes(x = wt_psi, y = ss5_strength)) +
  geom_point(size = 0.5, alpha = 0.3) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_wt_psi$r[2], 2),
                          ", p = ", signif(cor_wt_psi$p_value[2], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "WT PSI — Mutagenesis", y = "5' SS strength")

ggsave(file.path(plot_dir, "ss5_vs_wt_psi.png"),
       plot = p_ss5_vs_wt_psi, height = 3, width = 3)

p_zmean_vs_wt_psi <- ggplot(psi, aes(x = wt_psi, y = z_mean)) +
  geom_point(size = 0.5, alpha = 0.3) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_wt_psi$r[3], 2),
                          ", p = ", signif(cor_wt_psi$p_value[3], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "WT PSI — Mutagenesis", y = "Mean SS strength")

ggsave(file.path(plot_dir, "zmean_vs_wt_psi.png"),
       plot = p_zmean_vs_wt_psi, height = 3, width = 3)


# ── 5. SS strength vs overall splicing sensitivity ────────────────────────────

cor_sensitivity <- bind_rows(
  tidy_cor(psi$median_delta_logit, psi$ss3_strength, "3'SS strength vs sensitivity"),
  tidy_cor(psi$median_delta_logit, psi$ss5_strength, "5'SS strength vs sensitivity"),
  tidy_cor(psi$median_delta_logit, psi$z_mean,       "mean SS z-score vs sensitivity")
)
fwrite(cor_sensitivity, file.path(results_dir, "cor_sensitivity_vs_ss.tsv"), sep = '\t')

p_ss3_vs_sensitivity <- ggplot(psi, aes(x = median_delta_logit, y = ss3_strength)) +
  geom_point(size = 0.5, alpha = 0.3) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_sensitivity$r[1], 2),
                          ", p = ", signif(cor_sensitivity$p_value[1], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "Median |∆LogitPSI|", y = "3' SS strength")

ggsave(file.path(plot_dir, "ss3_vs_sensitivity.png"),
       plot = p_ss3_vs_sensitivity, height = 3, width = 3)

p_ss5_vs_sensitivity <- ggplot(psi, aes(x = median_delta_logit, y = ss5_strength)) +
  geom_point(size = 0.5, alpha = 0.3) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_sensitivity$r[2], 2),
                          ", p = ", signif(cor_sensitivity$p_value[2], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "Median |∆LogitPSI|", y = "5' SS strength")

ggsave(file.path(plot_dir, "ss5_vs_sensitivity.png"),
       plot = p_ss5_vs_sensitivity, height = 3, width = 3)

p_zmean_vs_sensitivity <- ggplot(psi, aes(x = median_delta_logit, y = z_mean)) +
  geom_point(size = 0.5, alpha = 0.3) +
  annotate("text", x = -Inf, y = -Inf,
           label = paste0("r = ", round(cor_sensitivity$r[3], 2),
                          ", p = ", signif(cor_sensitivity$p_value[3], 3)),
           hjust = -0.05, vjust = -0.5, size = 3.5, family = 'Helvetica') +
  theme_ss +
  labs(x = "Median |∆LogitPSI|", y = "Mean SS strength")

ggsave(file.path(plot_dir, "zmean_vs_sensitivity.png"),
       plot = p_zmean_vs_sensitivity, height = 3, width = 3)


# ── 6. Save per-exon splice-site info (EXON_SS_INFO_FILE) ─────────────────────
# This table is used by downstream scripts (03.2, 06.x, etc.).
# It joins per-exon SS scores + z-scores with VastDB average PSI.

meta   <- fread(METADATA_FILE, select = c("exon_id", "vastdb_event"))
vastdb <- fread(VASTDB_FILE,   select = c("EVENT", "Average"))

exon_ss_info <- psi %>%
  select(exon_id, ss5, ss3, ss5_strength, ss3_strength,
         z5, z3, z_mean, p5, p3, p_geom, z_min) %>%
  distinct(exon_id, .keep_all = TRUE) %>%
  left_join(meta,   by = "exon_id") %>%
  left_join(vastdb, by = c("vastdb_event" = "EVENT")) %>%
  rename(
    ss_5_seq        = ss5,
    ss_3_seq        = ss3,
    ss_5_strength   = ss5_strength,
    ss_3_strength   = ss3_strength,
    psi_vastdb_mean = Average
  ) %>%
  select(exon_id, psi_vastdb_mean, ss_3_seq, ss_3_strength,
         ss_5_seq, ss_5_strength, z5, z3, z_mean, p5, p3, p_geom, z_min)

message("Saving per-exon SS info to ", EXON_SS_INFO_FILE, " ...")
fwrite(exon_ss_info, EXON_SS_INFO_FILE, sep = "\t")
message("Written: ", EXON_SS_INFO_FILE)
