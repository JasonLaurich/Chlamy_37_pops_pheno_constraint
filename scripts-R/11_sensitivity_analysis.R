# Jason R Laurich

# July 2026

# lightly edited for re-submission on August 16th, 2026. 

# Sensitivity analysis: Li et al. (2019) empty-space test across data-inclusion thresholds.
# For ALL 15 panels (Figs 2 & 3), runs the polygonal empty-space test at 6 thresholds:
#   100%, 66.67%, 50%, 33.33% (current cut-off), 25%, and 10% of data retained,
#   where retention is based on scaled Euclidean distance from the minimum trait values
#   (i.e., always keeping the points closest to the Pareto front).

# Inputs:  processed-data/27_summary_table.csv
# Outputs:
#   processed-data/28_sensitivity_pvalues.csv         (long-format results)
#   processed-data/29_sensitivity_grid_table.csv      (grid table for supplement)
#   figures-supplemental/02_fig_s2_sensitivity.pdf  (15 panels)


# Packages & Functions --------------------------------------------------------

library(tidyverse)
library(pracma)    # polyarea()
library(cowplot)

par_frt <- function(df, xvar, yvar) {
  df <- df[order(-df[[xvar]], df[[yvar]]), ]
  pareto_points <- df[1, ]
  for (i in 2:nrow(df)) {
    if (df[i, yvar] >= tail(pareto_points[[yvar]], 1)) {
      pareto_points <- rbind(pareto_points, df[i, ])
    }
  }
  return(pareto_points)
}

# Li et al. (2019) polygonal empty-space test. Expects df with columns z.x and z.y.
li_test <- function(df, n_perm = 1000) {
  df <- df[!is.na(df$z.x) & !is.na(df$z.y), ]
  if (nrow(df) < 4) return(NA_real_)

  par.res <- par_frt(df, xvar = "z.x", yvar = "z.y") %>% arrange(z.x)
  x.max <- max(df$z.x); y.max <- max(df$z.y)
  poly  <- par.res[, c("z.x", "z.y")] %>% add_row(z.x = x.max, z.y = y.max)
  a.emp <- polyarea(poly$z.x, poly$z.y)

  null.areas <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    sh <- df %>% mutate(
      z.x.sim = sample(z.x, replace = FALSE),
      z.y.sim = sample(z.y, replace = FALSE)
    )
    par.n <- par_frt(sh, xvar = "z.x.sim", yvar = "z.y.sim") %>% arrange(z.x.sim)
    poly.n <- par.n[, c("z.x.sim", "z.y.sim")] %>%
      add_row(z.x.sim = max(sh$z.x.sim), z.y.sim = max(sh$z.y.sim))
    null.areas[i] <- polyarea(poly.n$z.x.sim, poly.n$z.y.sim)
  }

  mean(null.areas >= a.emp)
}

# Filter to top X% of data by scaled Euclidean distance from (x_min, y_min).
# df must already have a 'distance' column and be sorted ascending by distance.
filter_threshold <- function(df, keep_frac) {
  if (keep_frac >= 1.0) return(df)
  n_drop <- floor((1 - keep_frac) * nrow(df))
  df %>% arrange(distance) %>% slice((n_drop + 1):n())
}


# Load data ------------------------------------------------------------------

df <- read.csv("processed-data/27_summary_table.csv")

df$Evol.plt <- factor(df$Evol,
                      levels = c("none", "L", "N", "P", "S", "B", "BS", "C"),
                      labels = c("Ancestral", "Light limitation", "Nitrogen limitation",
                                 "Phosphorus limitation", "Salt stress",
                                 "Biotic depletion", "Biotic depletion x Salt", "Control"))


# Panel definitions ----------------------------------------------------------
# Each panel: panel_id, short label for plots, x-column, y-column

panels <- list(
  
  # Light
  
  list(panel_id = "L", row_grad = "L", col_grad = "L",
       label = "A — Light (L)",
       xcol = "I.µ.max", ycol = "I.comp"),
  
  list(panel_id = "LN", row_grad = "L", col_grad = "N",
       label = "B — L & N",
       xcol = "N.comp", ycol = "I.comp"),
  
  list(panel_id = "LP", row_grad = "L", col_grad = "P",
       label = "C — L & P",
       xcol = "P.comp", ycol = "I.comp"),
  
  list(panel_id = "LS", row_grad = "L", col_grad = "S",
       label = "D — L & S",
       xcol = "S.c", ycol = "I.comp"),
  
  list(panel_id = "LT", row_grad = "L", col_grad = "T",
       label = "E — L & T",
       xcol = "T.br", ycol = "I.comp"),
       
  # Nitrogen
  
  list(panel_id = "N", row_grad = "N", col_grad = "N",
       label = "F — Nitrogen (N)",
       xcol = "N.µ.max", ycol = "N.comp"),
  
  list(panel_id = "NP", row_grad = "N", col_grad = "P",
       label = "G — N & P",
       xcol = "P.comp", ycol = "N.comp"),
  
  list(panel_id = "NS", row_grad = "N", col_grad = "S",
       label = "H — N & S",
       xcol = "S.c", ycol = "N.comp"),
  
  list(panel_id = "NT", row_grad = "N", col_grad = "T",
       label = "I — N & T",
       xcol = "T.br", ycol = "N.comp"),
  
  # Phosphorus
  
  list(panel_id = "P", row_grad = "P", col_grad = "P",
       label = "J — Phosphorus (P)",
       xcol = "P.µ.max", ycol = "P.comp"),
  
  list(panel_id = "PS", row_grad = "P", col_grad = "S",
       label = "K — P & S",
       xcol = "S.c", ycol = "P.comp"),
  
  list(panel_id = "PT", row_grad = "P", col_grad = "T",
       label = "L — P & T",
       xcol = "T.br", ycol = "P.comp"),
  
  # Salt
  
  list(panel_id = "S", row_grad = "S", col_grad = "S",
       label = "M — Salt (S)",
       xcol = "S.µ.max", ycol = "S.c"),
  
  list(panel_id = "ST", row_grad = "S", col_grad = "T",
       label = "N — S & T",
       xcol = "T.br", ycol = "S.c"),
  
  # Temperature
  
  list(panel_id = "T", row_grad = "T", col_grad = "T",
       label = "O — Temperature (T)",
       xcol = "T.µ.max", ycol = "T.br")
  
)

# Thresholds -----------------------------------------------------------------

thresholds    <- c(1.0, 0.6667, 0.5, 0.3333, 0.25, 0.10)
thresh_labels <- c("100%", "66.67%", "50%", "33.33%", "25%", "10%")


# Run sensitivity analysis ---------------------------------------------------

set.seed(123)

results <- vector("list", length(panels) * length(thresholds))
idx <- 1

for (p in panels) {
  cat("Processing panel:", p$panel_id, "\n")

  # Prepare panel-specific dataset with z.x, z.y, and distance
  df.base <- df %>%
    select(all_of(c(p$xcol, p$ycol)), everything()) %>%
    filter(!is.na(.data[[p$xcol]]), !is.na(.data[[p$ycol]])) %>%
    mutate(
      z.x  = .data[[p$xcol]],
      z.y  = .data[[p$ycol]],
      z.x2 = scale(z.x)[, 1],
      z.y2 = scale(z.y)[, 1]
    ) %>%
    mutate(
      distance = sqrt((z.x2 - min(z.x2, na.rm = TRUE))^2 +
                      (z.y2 - min(z.y2, na.rm = TRUE))^2)
    ) %>%
    arrange(distance)

  for (t_idx in seq_along(thresholds)) {
    keep  <- thresholds[t_idx]
    tlabel <- thresh_labels[t_idx]

    df.filt <- filter_threshold(df.base, keep)
    n_obs   <- nrow(df.filt)
    p_val   <- li_test(df.filt, n_perm = 1000)

    cat("  ", tlabel, "(n =", n_obs, ") → P =",
        ifelse(is.na(p_val), "NA", round(p_val, 4)), "\n")

    results[[idx]] <- data.frame(
      panel_id   = p$panel_id,
      label      = p$label,
      row_grad   = p$row_grad,
      col_grad   = p$col_grad,
      threshold  = tlabel,
      pct_kept   = keep * 100,
      n_obs      = n_obs,
      p_value    = p_val,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

results_df <- bind_rows(results)
results_df$threshold <- factor(results_df$threshold, levels = thresh_labels)

write.csv(results_df, "processed-data/28_sensitivity_pvalues.csv", row.names = FALSE) # Save long-format results

# Per-panel P-value plots ----------------------------------------------------

results_df <- read.csv("processed-data/28_sensitivity_pvalues.csv")

gradients  <- c("L", "N", "P", "S", "T")
grad_names <- c("Light", "Nitrogen", "Phosphorus", "Salt", "Temperature")

# Build a lookup for panel titles with sequential letters (upper triangle, row by row)
panel_titles <- list()
letter_i <- 1
for (row_i in seq_along(gradients)) {
  for (col_j in seq_along(gradients)) {
    if (col_j >= row_i) {
      pid <- if (col_j == row_i) gradients[row_i] else paste0(gradients[row_i], gradients[col_j])
      label <- if (col_j == row_i) {
        paste0(LETTERS[letter_i], ") ", grad_names[row_i])
      } else {
        paste0(LETTERS[letter_i], ") ", grad_names[row_i], "\u2014", grad_names[col_j])
      }
      panel_titles[[pid]] <- label
      letter_i <- letter_i + 1
    }
  }
}

# Theme matching main figures
theme_grid <- theme(
  plot.title  = element_text(size = 10, face = "bold", hjust = 0.03),
  axis.title  = element_text(size = 10),
  axis.text   = element_text(size = 10),
  axis.text.x = element_text(angle = 35, hjust = 1),
  plot.margin = unit(c(2, 2, 2, 2), "pt")
)

# Plot function
make_grid_plot <- function(pid) {
  df.plot <- results_df %>% filter(panel_id == pid)
  ggplot(df.plot, aes(x = pct_kept, y = p_value)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
    geom_line(linewidth = 0.6, color = "grey40") +
    geom_point(aes(fill = p_value <= 0.05), shape = 21, size = 2,
               color = "grey30", stroke = 0.4) +
    scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = "none") +
    scale_x_reverse(breaks = c(100, 66.67, 50, 33.33, 25, 10),
                    labels = c("100", "67", "50", "33", "25", "10")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    labs(title = panel_titles[[pid]],
         x = "% data retained",
         y = expression(italic(P)*"-value")) +
    theme_classic() +
    theme_grid
}

# Build 5×5 grid
grid_list <- vector("list", 25)
k <- 1
for (row_i in seq_along(gradients)) {
  for (col_j in seq_along(gradients)) {
    if (col_j < row_i) {
      grid_list[[k]] <- ggplot() + theme_void()                    # blank lower triangle
    } else if (col_j == row_i) {
      grid_list[[k]] <- make_grid_plot(gradients[row_i])           # diagonal
    } else {
      grid_list[[k]] <- make_grid_plot(paste0(gradients[row_i], gradients[col_j]))  # off-diagonal
    }
    k <- k + 1
  }
}

combined <- plot_grid(plotlist = grid_list, ncol = 5, align = "hv",
                      rel_widths  = rep(1, 5),
                      rel_heights = rep(1, 5))

ggsave("figures-supplemental/02_fig_s2_sensitivity.pdf", plot = combined, width = 13, height = 13)

# Grid table of P-values -----------------------------------------------------
# Layout mirrors Table S3: rows = Gradient × Threshold, columns = L, N, P, S, T
# Diagonal = Fig2 intra-gradient panels; upper triangle = Fig3 inter-gradient panels
# Lower triangle = blank (same pair as upper triangle, different orientation not tested)

gradients <- c("L", "N", "P", "S", "T")

# Lookup: (row_grad, col_grad) → panel_id
panel_lookup <- results_df %>%
  distinct(panel_id, row_grad, col_grad) %>%
  { setNames(.$panel_id, paste(.$row_grad, .$col_grad, sep = "_")) }

table_rows <- list()

for (row_grad in gradients) {
  for (t_idx in seq_along(thresholds)) {
    tlabel <- thresh_labels[t_idx]

    row <- data.frame(Gradient = row_grad, Threshold = tlabel,
                      stringsAsFactors = FALSE)

    for (col_grad in gradients) {

      # Only fill diagonal and upper triangle
      grad_order <- match(c(row_grad, col_grad), gradients)
      if (grad_order[1] > grad_order[2]) {
        row[[col_grad]] <- ""   # lower triangle: blank
        next
      }

      key <- paste(row_grad, col_grad, sep = "_")
      pid <- panel_lookup[key]

      if (!is.na(pid)) {
        pval <- results_df %>%
          filter(panel_id == pid, threshold == tlabel) %>%
          pull(p_value)
        if (length(pval) > 0 && !is.na(pval)) {
          # Format: bold < 0.05
          row[[col_grad]] <- ifelse(pval == 0, "<0.001",
                                    format(round(pval, 3), nsmall = 3))
        } else {
          row[[col_grad]] <- NA_character_
        }
      } else {
        row[[col_grad]] <- NA_character_
      }
    }

    table_rows[[length(table_rows) + 1]] <- row
  }
}

grid_table <- bind_rows(table_rows) %>%
  select(Gradient, Threshold, all_of(gradients))

write.csv(grid_table, "processed-data/29_sensitivity_grid_table.csv", row.names = FALSE)