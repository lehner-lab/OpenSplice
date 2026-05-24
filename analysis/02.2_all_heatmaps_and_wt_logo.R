library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggh4x)
library(patchwork)
library(stringr)
library(ggseqlogo)
library(cowplot)
library(here)

source(here("analysis", "config.R"))

heatmaps_df <- fread(MASTER_TABLE, sep = "\t")

heatmaps_df[, region_start := factor(region_start, levels = c("Intron up", "Exon", "Intron down"))]
heatmaps_df[, mut := factor(mut, levels = c("∆21nt","∆6nt","∆3nt","∆1nt","A","G","C","U"))]

wt_logo_df <- heatmaps_df %>%
  filter(mut_type == "sub", start == end) %>%
  mutate(delta_logit = ifelse(abs(delta_psi) < 2.5,0.1,delta_logit)) %>%
  transmute(
    exon_id,
    wt_psi,
    pos = start,
    wt = str_replace_all(wt, "U", "T"),
    delta_logit
  ) %>%
  group_by(exon_id, wt_psi, pos, wt) %>%
  summarise(
    wt_score = -mean(delta_logit, na.rm = TRUE),
    .groups = "drop"
  )

info_ss <- fread(file.path(SUP_TABLES_DIR, "Supplementary_Table1.tsv"), sep = "\t")

delta_colours <- colorRampPalette(
  c("#2171B5", "#4292C6", "#A6CEE3", "white", "#FDAEAE", "#FB6A4A", "#EF3B2C")
)(n = 200)

nt_colours <- c(A = "#1A9641", C = "#4575B4", G = "#FDAE61", T = "#D73027")

mut_list <- data.frame(
  y = c("∆21nt", "∆6nt", "∆3nt", "∆1nt", "A", "G", "C", "U"),
  label = c("del21","del6","del3","del1","A","G","C","U")
)

# Exclude exons without a WT PSI measurement before ordering
exon_list <- unique(heatmaps_df[, .(exon_id, wt_psi)])[!is.na(wt_psi)][order(wt_psi)]$exon_id

fill_scale <- scale_fill_gradientn(
  colours = delta_colours,
  limits = c(-100, 100),
  breaks = c(-100,-50,-25,-5,0,5,25,50,100),
  trans = scales::pseudo_log_trans(sigma = 10),
  oob = scales::squish,
  guide = guide_colourbar(
    title.position = "top",
    title.hjust = 0.5,
    barwidth = unit(3, "mm"),
    barheight = unit(45, "mm"),
    ticks = TRUE
  ),
  na.value = "grey60"
)

base_theme <- theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 12, color = "black", margin = margin(t = 3)),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size = 9, color = "black", margin = margin(b = 3)),
    strip.placement = "outside",
    legend.text = element_text(size = 7, color = "black"),
    legend.title = element_text(size = 8, color = "black"),
    legend.position = "right",
    legend.justification = "right",
    legend.box.just = "right",
    legend.direction = "vertical",
    panel.spacing = unit(0, "mm"),
    plot.margin = margin(t = 10, r = 7.5, b = 10, l = 7.5, unit = "pt")
  )

get_L_ex <- function(dt_ex) {
  wt_pos <- dt_ex %>%
    filter(grepl("wt", variant_id), region_start == "Exon") %>%
    summarise(L = max(start, na.rm = TRUE)) %>%
    pull(L)

  if (length(wt_pos) == 0 || is.na(wt_pos)) wt_pos <- 70
  max(71, min(250, wt_pos))
}

build_scheme_plot <- function(L_ex) {
  lvls <- c("Intron up", "Exon", "Intron down")
  up_rng <- 1:70
  ex_rng <- 71:L_ex
  down_rng <- (L_ex + 1):250

  seg_df <- data.frame(
    region_start = factor(c("Intron up", "Intron down"), levels = lvls),
    x = c(1, L_ex + 1),
    xend = c(70, L_ex + 25),
    y = 0,
    yend = 0
  )

  rect_df <- data.frame(
    region_start = factor("Exon", levels = lvls),
    xmin = 70.5,
    xmax = L_ex + 0.5,
    ymin = -1.2,
    ymax = 1.2
  )

  pad_df <- bind_rows(
    data.frame(x = c(min(up_rng), max(up_rng)), region_start = "Intron up"),
    data.frame(x = c(min(ex_rng), max(ex_rng)), region_start = "Exon"),
    data.frame(x = c(min(down_rng), max(down_rng)), region_start = "Intron down")
  ) %>%
    mutate(region_start = factor(region_start, levels = lvls))

  label_df <- data.frame(
    region_start = factor(c("Intron up", "Exon", "Intron down"), levels = lvls),
    x = c(35.5, (71 + L_ex) / 2, L_ex + 7),
    label = c("Intron", "Exon", "Intron"),
    hjust = c(0.5, 0.5, 0)
  )

  ggplot() +
    geom_blank(data = pad_df, aes(x = x)) +
    geom_segment(
      data = seg_df,
      aes(x = x, xend = xend, y = 0, yend = 0),
      color = "black",
      linewidth = 0.5,
      inherit.aes = FALSE
    ) +
    geom_rect(
      data = rect_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey70",
      color = "black",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_vline(
      xintercept = c(70.5, L_ex + 0.5),
      color = "black",
      linewidth = 0.5,
      linetype = "dashed"
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = -1.45, label = label, hjust = hjust),
      vjust = 1,
      size = 3,
      color = "black",
      inherit.aes = FALSE
    ) +
    facet_grid(cols = vars(region_start), scales = "free", space = "free_x") +
    ggh4x::facetted_pos_scales(x = list(
      region_start == "Intron up" ~ scale_x_continuous(
        limits = c(min(up_rng) - 0.5, max(up_rng) + 0.5), expand = c(0, 0)),
      region_start == "Exon" ~ scale_x_continuous(
        limits = c(min(ex_rng) - 0.5, max(ex_rng) + 0.5), expand = c(0, 0)),
      region_start == "Intron down" ~ scale_x_continuous(
        limits = c(min(down_rng) - 0.5, max(down_rng) + 0.5), expand = c(0, 0))
    )) +
    scale_y_continuous(limits = c(-1.5, 1.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      strip.text.x = element_blank(),
      panel.spacing = unit(0, "mm"),
      plot.margin = margin(t = 10, r = 7.5, b = 15, l = 7.5, unit = "pt")
    )
}

build_logo_plot <- function(exon_name, wt_logo_df, L_ex, base_theme) {
  logo_ex <- wt_logo_df %>% filter(exon_id == exon_name)
  if (nrow(logo_ex) == 0) return(NULL)

  lvls <- c("Intron up", "Exon", "Intron down")
  up_rng <- 1:70
  ex_rng <- 71:L_ex
  down_rng <- (L_ex + 1):250
  bases_rna <- c("A", "C", "G", "U")

  mat <- matrix(0, nrow = 4, ncol = 250,
                dimnames = list(bases_rna, as.character(1:250)))

  for (i in seq_len(nrow(logo_ex))) {
    w_rna <- ifelse(logo_ex$wt[i] == "T", "U", logo_ex$wt[i])
    p <- logo_ex$pos[i]

    if (w_rna %in% bases_rna && p >= 1 && p <= 250) {
      mat[w_rna, p] <- logo_ex$wt_score[i]
    }
  }

  cutoff <- min(L_ex + 25, 250)
  if (cutoff < 250) mat[, (cutoff + 1):250] <- 0

  y_lim <- c(
    min(0, logo_ex$wt_score, na.rm = TRUE),
    max(0, logo_ex$wt_score, na.rm = TRUE)
  )

  cs <- ggseqlogo::make_col_scheme(
    chars = c("A", "C", "G", "U"),
    cols = unname(nt_colours[c("A", "C", "G", "T")])
  )

  add_region <- function(dat) {
    if (!is.data.frame(dat) || !"x" %in% names(dat)) return(dat)

    dat$region_start <- factor(
      dplyr::case_when(
        dat$x <= 70 ~ "Intron up",
        dat$x <= L_ex ~ "Exon",
        TRUE ~ "Intron down"
      ),
      levels = lvls
    )

    dat
  }

  p <- suppressWarnings(
    ggseqlogo::ggseqlogo(mat, method = "custom", col_scheme = cs)
  )

  p$data <- add_region(p$data)
  p$layers <- lapply(p$layers, function(layer) {
    layer$data <- add_region(layer$data)
    layer
  })

  pad_df <- bind_rows(
    data.frame(x = c(min(up_rng), max(up_rng)), region_start = "Intron up"),
    data.frame(x = c(min(ex_rng), max(ex_rng)), region_start = "Exon"),
    data.frame(x = c(min(down_rng), max(down_rng)), region_start = "Intron down")
  ) %>%
    mutate(region_start = factor(region_start, levels = lvls))

  zero_line_df <- data.frame(
    region_start = factor(c("Intron up", "Exon", "Intron down"), levels = lvls),
    x = c(1, 71, L_ex + 1),
    xend = c(70, L_ex, cutoff),
    y = 0,
    yend = 0
  )

  vline_df <- data.frame(
    region_start = factor(c("Intron up", "Exon"), levels = lvls),
    xintercept = c(70.5, L_ex + 0.5)
  )

  p +
    geom_blank(data = pad_df, aes(x = x), inherit.aes = FALSE) +
    geom_segment(
      data = zero_line_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.3,
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_vline(
      data = vline_df,
      aes(xintercept = xintercept),
      color = "black",
      linewidth = 0.5,
      linetype = "dashed",
      inherit.aes = FALSE
    ) +
    facet_grid(cols = vars(region_start), scales = "free", space = "free_x") +
    ggh4x::facetted_pos_scales(x = list(
      region_start == "Intron up" ~ scale_x_continuous(
        limits = c(min(up_rng) - 0.5, max(up_rng) + 0.5), expand = c(0, 0)),
      region_start == "Exon" ~ scale_x_continuous(
        limits = c(min(ex_rng) - 0.5, max(ex_rng) + 0.5), expand = c(0, 0)),
      region_start == "Intron down" ~ scale_x_continuous(
        limits = c(min(down_rng) - 0.5, max(down_rng) + 0.5), expand = c(0, 0))
    )) +
    coord_cartesian(ylim = y_lim, clip = "off") +
    theme(
      strip.text.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(color = "black", linewidth = 0.3),
      axis.ticks.y = element_line(color = "black", linewidth = 0.3),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.title.y = element_text(size = 8, color = "black"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "mm"),
      legend.position = "none",
      plot.margin = margin(t = 0, r = 7.5, b = 10, l = 7.5, unit = "pt")
    ) +
    labs(y = expression(-mean(Delta*LogitPSI[sub])))
}

build_exon_plot <- function(exon_id, dt_ex, info_ss, mut_list, fill_scale, base_theme,
                            show_legend = TRUE, L_ex = NULL) {
  exon_key <- exon_id
  mut_df <- as.data.frame(dt_ex)
  if (nrow(mut_df) == 0) return(NULL)

  mut_df <- mut_df %>%
    mutate(
      region_start = factor(region_start, levels = c("Intron up","Exon","Intron down")),
      mut = factor(mut, levels = c("∆21nt","∆6nt","∆3nt","∆1nt","A","G","C","U")),
      mut_group = ifelse(mut %in% c("∆1nt","∆3nt","∆6nt","∆21nt"), "del", "sub"),
      mut_group = factor(mut_group, levels = c("sub","del")),
      start = as.integer(as.character(start)),
      wt = ifelse(wt == "T", "U", wt),
      delta_psi = as.numeric(delta_psi)
    )

  ex_title <- gsub("_e", " exon ", exon_key)

  wt_psi_raw <- mut_df %>%
    filter(grepl("wt", variant_id)) %>%
    distinct(wt_psi) %>%
    pull(wt_psi)

  if (length(wt_psi_raw) == 0 || all(is.na(wt_psi_raw))) {
    wt_psi_raw <- unique(mut_df$wt_psi)[1]
  } else {
    wt_psi_raw <- wt_psi_raw[!is.na(wt_psi_raw)][1]
  }

  wt_psi <- round(wt_psi_raw, 1)

  info_coord <- as.data.frame(info_ss) %>%
    filter(.data$exon_id == exon_key) %>%
    slice(1)

  coord <- paste0(
    "chr", info_coord$chr, ":",
    info_coord$coord_start_mutagenesis, "-",
    info_coord$coord_end_mutagenesis
  )

  wt_df <- mut_df %>%
    filter(grepl("wt", variant_id)) %>%
    select(start, wt, region_start) %>%
    distinct()

  wt_labels <- rep("", 250)
  ok <- wt_df$start[wt_df$start >= 1 & wt_df$start <= 250]
  wt_labels[ok] <- wt_df$wt[match(ok, wt_df$start)]
  names(wt_labels) <- as.character(1:250)

  if (is.null(L_ex)) {
    L_ex <- wt_df %>%
      filter(region_start == "Exon") %>%
      summarise(L = max(start, na.rm = TRUE)) %>%
      pull(L)

    if (is.na(L_ex)) L_ex <- 70
    L_ex <- max(71, min(250, L_ex))
  }

  up_rng <- 1:70
  ex_rng <- 71:L_ex
  down_rng <- (L_ex + 1):250

  tile_df <- mut_df %>%
    filter(!is.na(mut_group), !is.na(mut), !is.na(start)) %>%
    filter(start >= 1, start <= 250)

  pad_df <- bind_rows(
    data.frame(start = c(min(up_rng), max(up_rng)), mut = levels(mut_df$mut)[1],
               region_start = "Intron up", mut_group = levels(mut_df$mut_group)[1]),
    data.frame(start = c(min(ex_rng), max(ex_rng)), mut = levels(mut_df$mut)[1],
               region_start = "Exon", mut_group = levels(mut_df$mut_group)[1]),
    data.frame(start = c(min(down_rng), max(down_rng)), mut = levels(mut_df$mut)[1],
               region_start = "Intron down", mut_group = levels(mut_df$mut_group)[1])
  ) %>%
    mutate(
      region_start = factor(region_start, levels = levels(mut_df$region_start)),
      mut_group = factor(mut_group, levels = levels(mut_df$mut_group)),
      mut = factor(mut, levels = levels(mut_df$mut))
    )

  label_df <- data.frame(
    start = c(61, 71, L_ex + 1),
    mut = 9.2,
    label = c(
      paste0("ss3: ", info_ss$ss_3_strength[info_ss$exon_id == exon_key][1]),
      "",
      paste0("ss5: ", info_ss$ss_5_strength[info_ss$exon_id == exon_key][1])
    ),
    region_start = factor(c("Intron up","Exon","Intron down"),
                          levels = levels(mut_df$region_start)),
    mut_group = factor(c("sub","sub","sub"), levels = levels(mut_df$mut_group)),
    hjust = c(0,0,0),
    vjust = c(1,1.1,1.1)
  )

  wt_labels_local <- wt_labels

  lab_fun_local <- function(x) {
    x <- as.integer(round(x))
    out <- wt_labels_local[as.character(x)]
    out[is.na(out)] <- ""
    unname(out)
  }

  p <- ggplot() +
    geom_blank(data = pad_df, aes(x = start)) +
    geom_tile(
      color = "grey80",
      data = tile_df,
      aes(x = start, y = mut, fill = delta_psi)
    ) +
    fill_scale +
    geom_text(
      data = label_df,
      aes(x = start, y = mut, label = label, hjust = hjust, vjust = vjust),
      inherit.aes = FALSE,
      size = 3
    ) +
    scale_y_discrete(breaks = mut_list$y, labels = mut_list$label) +
    facet_grid(
      cols = vars(region_start),
      scales = "free",
      space = "free_x",
      labeller = labeller(region_start = c(
        "Intron up" = "Intron",
        "Exon" = "Exon",
        "Intron down" = "Intron"
      ))
    ) +
    ggh4x::facetted_pos_scales(x = list(
      region_start == "Intron up" ~ scale_x_continuous(
        limits = c(min(up_rng) - 0.5, max(up_rng) + 0.5),
        breaks = up_rng,
        labels = lab_fun_local,
        expand = c(0, 0)
      ),
      region_start == "Exon" ~ scale_x_continuous(
        limits = c(min(ex_rng) - 0.5, max(ex_rng) + 0.5),
        breaks = ex_rng,
        labels = lab_fun_local,
        expand = c(0, 0)
      ),
      region_start == "Intron down" ~ scale_x_continuous(
        limits = c(min(down_rng) - 0.5, max(down_rng) + 0.5),
        breaks = down_rng,
        labels = lab_fun_local,
        expand = c(0, 0)
      )
    )) +
    base_theme +
    theme(
      strip.text.x = element_blank(),
      plot.margin = margin(t = 10, r = 7.5, b = 10, l = 7.5, unit = "pt")
    ) +
    ylab("Mutation") +
    labs(fill = expression(Delta*PSI)) +
    ggtitle(paste0(ex_title, " (WT PSI = ", wt_psi, "% - ", coord, ")")) +
    geom_vline(
      xintercept = c(70.5, L_ex + 0.5),
      linetype = "dashed",
      color = "black",
      linewidth = 0.5
    )

  if (!show_legend) {
    p <- p + theme(
      legend.position = "none",
      axis.title.y = element_blank()
    )
  }

  p
}

dir.create(here("figures", "02_heatmaps"), showWarnings = FALSE, recursive = TRUE)
out_pdf <- here("figures", "02_heatmaps", "Heatmaps_all_exons.pdf")

pdf(out_pdf, width = 18, height = 4, onefile = TRUE)

for (i in seq_along(exon_list)) {
  exon_name <- exon_list[i]
  dt_ex <- heatmaps_df %>% filter(exon_id == exon_name)

  L_ex <- get_L_ex(dt_ex)

  p_heat <- build_exon_plot(
    exon_id = exon_name,
    dt_ex = dt_ex,
    info_ss = info_ss,
    mut_list = mut_list,
    fill_scale = fill_scale,
    base_theme = base_theme,
    show_legend = (i == 1),
    L_ex = L_ex
  )

  p_logo <- if (i <= 2) NULL else build_logo_plot(exon_name, wt_logo_df, L_ex, base_theme)
  p_scheme <- build_scheme_plot(L_ex)

  if (!is.null(p_heat) && !is.null(p_logo)) {
    aligned <- cowplot::align_plots(p_heat, p_scheme, p_logo, align = "v", axis = "lr")
    print(cowplot::plot_grid(plotlist = aligned, ncol = 1, rel_heights = c(1, 0.16, 0.4)))
  } else if (!is.null(p_heat)) {
    aligned <- cowplot::align_plots(p_heat, p_scheme, align = "v", axis = "lr")
    print(cowplot::plot_grid(plotlist = aligned, ncol = 1, rel_heights = c(1, 0.16)))
  }
}

dev.off()

