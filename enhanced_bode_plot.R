###############################################################################
# Enhanced Bode plot script with -3dB cutoff lines and improved annotations
# Formatted for PLOS Computational Biology submission
# - Minimal decimal places in CSV output
# - Bold text for key plot elements
###############################################################################

library(ggplot2)
library(patchwork)
library(pracma)
library(stringr)
library(svglite)
library(scales)

# --------------------------
# User parameters
# --------------------------
sequence1 <- "UGGCAGUGUCUUAGCUGGUUGU"
sequence2 <- "AGCUGCUGUUGACAGUGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCU"

# Base molecular weights (g/mol)
Adenine <- 135.13
Guanine <- 151.13
Cytosine <- 111.10
Uracil <- 112.09

omegan_base <- 0.1
zeta_base  <- 1.0

freq_range_rad <- c(1e-3, 1e3)
n_points <- 2000

out_root <- "bode_outputs"

make_dir <- function(...) {
  dir <- file.path(..., fsep = .Platform$file.sep)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

# --------------------------
# Derived base parameters
# --------------------------
ratioAG <- Adenine / Guanine
ratioAC <- Adenine / Cytosine
ratioAU <- Adenine / Uracil

omeganA <- omegan_base
omeganG <- omegan_base * ratioAG
omeganC <- omegan_base * ratioAC
omeganU <- omegan_base * ratioAU

zetaA <- zeta_base
zetaG <- zeta_base * ratioAG
zetaC <- zeta_base * ratioAC
zetaU <- zeta_base * ratioAU

get_base_params <- function(base) {
  switch(base,
         "A" = list(omega = omeganA, zeta = zetaA),
         "G" = list(omega = omeganG, zeta = zetaG),
         "C" = list(omega = omeganC, zeta = zetaC),
         "U" = list(omega = omeganU, zeta = zetaU),
         stop(paste("Unknown base:", base)))
}

unwrap_phase <- function(phase_rad) {
  n <- length(phase_rad)
  if (n <= 1) return(phase_rad)
  unwrapped <- numeric(n)
  unwrapped[1] <- phase_rad[1]
  for (i in 2:n) {
    delta <- phase_rad[i] - phase_rad[i-1]
    delta <- delta - 2*pi * round(delta / (2*pi))
    unwrapped[i] <- unwrapped[i-1] + delta
  }
  return(unwrapped)
}

parse_sequence_params <- function(seq_string) {
  bases <- str_split(seq_string, "")[[1]]
  omegas <- numeric(length(bases)); zetas <- numeric(length(bases))
  for (i in seq_along(bases)) {
    p <- get_base_params(bases[i])
    omegas[i] <- p$omega
    zetas[i] <- p$zeta
  }
  k2 <- 2 * zetas * omegas
  k3 <- omegas
  k4 <- omegas
  list(bases = bases, k2 = k2, k3 = k3, k4 = k4)
}

seq1_params <- parse_sequence_params(sequence1)
seq2_params <- parse_sequence_params(sequence2)

# --------------------------
# Stage responses
# --------------------------
stage_response_G <- function(k2, k3, k4, omega_vec) {
  s <- 1i * omega_vec
  numerator <- k3 * (s + k2)
  denominator <- s * (s + k2) + k3 * k4
  numerator / denominator
}

stage_response_H <- function(k2, k3, k4, omega_vec) {
  s <- 1i * omega_vec
  numerator <- k2 * k4
  denominator <- s * (s + k2) + k3 * k4
  numerator / denominator
}

cascade_response <- function(stages_type, k2_vec, k3_vec, k4_vec, omega_vec) {
  if (length(stages_type) != length(k2_vec)) stop("stages_type must match k vectors")
  H_total <- rep(1+0i, length(omega_vec))
  for (i in seq_along(stages_type)) {
    if (stages_type[i] == "G") {
      H_i <- stage_response_G(k2_vec[i], k3_vec[i], k4_vec[i], omega_vec)
    } else if (stages_type[i] == "H") {
      H_i <- stage_response_H(k2_vec[i], k3_vec[i], k4_vec[i], omega_vec)
    } else stop("Unknown stage type")
    H_total <- H_total * H_i
  }
  mag_db <- 20 * log10(Mod(H_total))
  phase_rad <- Arg(H_total)
  phase_unwrapped <- unwrap_phase(phase_rad)
  phase_deg <- phase_unwrapped * 180 / pi
  list(H = H_total, magnitude = mag_db, phase = phase_deg, raw_phase = phase_rad)
}

# --------------------------
# Configuration builder
# --------------------------
build_configs <- function() {
  configs <- list()

  stages1 <- rep("G", length(seq1_params$k2))
  configs[[1]] <- list(id = "config1_seq1_G_only",
                       desc = "Config 1: 22-nt miRNA (G-type)",
                       short = "Config 1",
                       stages = stages1,
                       k2 = seq1_params$k2,
                       k3 = seq1_params$k3,
                       k4 = seq1_params$k4)

  stages2 <- rep("H", length(seq1_params$k2))
  configs[[2]] <- list(id = "config2_seq1_H_only",
                       desc = "Config 2: 22-nt miRNA (H-type)",
                       short = "Config 2",
                       stages = stages2,
                       k2 = seq1_params$k2,
                       k3 = seq1_params$k3,
                       k4 = seq1_params$k4)

  stages3 <- rep("G", length(seq2_params$k2))
  configs[[3]] <- list(id = "config3_seq2_G_only",
                       desc = "Config 3: 140-nt ncRNA (G-type)",
                       short = "Config 3",
                       stages = stages3,
                       k2 = seq2_params$k2,
                       k3 = seq2_params$k3,
                       k4 = seq2_params$k4)

  stages4 <- c(rep("G", length(seq1_params$k2)), rep("H", length(seq2_params$k2)))
  k2_4 <- c(seq1_params$k2, seq2_params$k2)
  k3_4 <- c(seq1_params$k3, seq2_params$k3)
  k4_4 <- c(seq1_params$k4, seq2_params$k4)
  configs[[4]] <- list(id = "config4_seq1G_seq2H",
                       desc = "Config 4: Mixed (Seq1-G + Seq2-H)",
                       short = "Config 4",
                       stages = stages4,
                       k2 = k2_4,
                       k3 = k3_4,
                       k4 = k4_4)

  configs
}

configs <- build_configs()

omega_rad <- 10^seq(log10(freq_range_rad[1]), log10(freq_range_rad[2]), length.out = n_points)
freq_hz <- omega_rad / (2 * pi)

# --------------------------
# Enhanced metrics with gain/phase margins
# Formatted with minimal decimal places
# --------------------------
compute_enhanced_metrics <- function(df) {
  magnitude <- df$magnitude
  freq <- df$frequency
  phase <- df$phase

  dc_gain <- magnitude[1]
  hf_gain <- magnitude[length(magnitude)]

  # Peak magnitude
  peak_idx <- which.max(magnitude)
  peak_mag <- magnitude[peak_idx]
  peak_freq <- freq[peak_idx]

  # -3dB cutoff frequency
  cutoff <- dc_gain - 3
  idx_3db <- which(magnitude < cutoff)[1]
  freq_3db <- if (!is.na(idx_3db)) freq[idx_3db] else NA

  # Gain margin: freq where phase = -180°, then check gain
  phase_crossover_idx <- which.min(abs(phase + 180))
  phase_crossover_freq <- freq[phase_crossover_idx]
  gain_at_crossover <- magnitude[phase_crossover_idx]
  gain_margin_dB <- -gain_at_crossover  # GM = -|G(jω)| at phase crossover

  # Phase margin: freq where |G| = 0 dB, then check phase
  gain_crossover_idx <- which.min(abs(magnitude))
  gain_crossover_freq <- freq[gain_crossover_idx]
  phase_at_crossover <- phase[gain_crossover_idx]
  phase_margin_deg <- 180 + phase_at_crossover  # PM = 180° + phase

  # Bandwidth (-3dB)
  bandwidth <- freq_3db

  phase_initial <- phase[1]
  phase_final <- phase[length(phase)]

  # Return data frame with appropriately rounded values
  data.frame(
    config = df$config[1],
    dc_gain_dB = round(dc_gain, 2),
    hf_gain_dB = round(hf_gain, 2),
    peak_mag_dB = round(peak_mag, 2),
    peak_freq = signif(peak_freq, 3),
    freq_3dB = signif(freq_3db, 3),
    bandwidth = signif(bandwidth, 3),
    gain_margin_dB = round(gain_margin_dB, 2),
    phase_crossover_freq = signif(phase_crossover_freq, 3),
    phase_margin_deg = round(phase_margin_deg, 1),
    gain_crossover_freq = signif(gain_crossover_freq, 3),
    phase_initial_deg = round(phase_initial, 1),
    phase_final_deg = round(phase_final, 1),
    stringsAsFactors = FALSE
  )
}

# --------------------------
# Enhanced plotting functions with bold text
# Formatted for PLOS Computational Biology
# --------------------------
plot_mag_enhanced <- function(df, metrics, title_text, x_label, show_annotations = TRUE) {
  p <- ggplot(df, aes(x = frequency, y = magnitude)) +
    geom_line(color = "#2E86AB", size = 1) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = x_label, y = "Magnitude (dB)", title = title_text) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_line(color = "gray90", size = 0.3),
      panel.grid.major = element_line(color = "gray80", size = 0.5)
    )

  if (show_annotations && !is.na(metrics$freq_3dB)) {
    # Add -3dB horizontal line
    p <- p + geom_hline(yintercept = metrics$dc_gain_dB - 3,
                        linetype = "dotted", color = "gray50", size = 0.5)

    # Add vertical line at cutoff frequency
    p <- p + geom_vline(xintercept = metrics$freq_3dB,
                        linetype = "dashed", color = "#A23B72", size = 0.8)

    # Add label for cutoff frequency with bold text
    label_text <- sprintf("ω[c] == %.2g", metrics$freq_3dB)

    p <- p + annotate("text",
                      x = metrics$freq_3dB * 2.5,
                      y = max(df$magnitude) * 0.9,
                      label = label_text,
                      parse = TRUE,
                      hjust = 0, size = 4, color = "#A23B72", fontface = "bold")

    # Add DC gain annotation with bold text
    p <- p + annotate("text",
                      x = min(df$frequency) * 2,
                      y = metrics$dc_gain_dB,
                      label = sprintf("DC: %.1f dB", metrics$dc_gain_dB),
                      hjust = 0, vjust = -0.5, size = 3.5, color = "gray30", fontface = "bold")
  }

  return(p)
}

plot_phase_enhanced <- function(df, metrics, title_text, x_label, show_annotations = TRUE) {
  p <- ggplot(df, aes(x = frequency, y = phase)) +
    geom_line(color = "#F18F01", size = 1) +
    scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
    labs(x = x_label, y = "Phase (degrees)") +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      axis.text = element_text(face = "bold", size = 10),
      panel.grid.minor = element_line(color = "gray90", size = 0.3),
      panel.grid.major = element_line(color = "gray80", size = 0.5)
    )

  if (show_annotations && !is.na(metrics$freq_3dB)) {
    # Add vertical line at cutoff frequency
    p <- p + geom_vline(xintercept = metrics$freq_3dB,
                        linetype = "dashed", color = "#A23B72", size = 0.8)

    # Add -180° reference line if relevant
    if (min(df$phase) < -180 && max(df$phase) > -180) {
      p <- p + geom_hline(yintercept = -180,
                          linetype = "dotted", color = "gray50", size = 0.5)
    }

    # Add phase at cutoff annotation with bold text
    phase_at_cutoff <- df$phase[which.min(abs(df$frequency - metrics$freq_3dB))]
    p <- p + annotate("text",
                      x = metrics$freq_3dB * 2.5,
                      y = phase_at_cutoff,
                      label = sprintf("φ = %.0f°", phase_at_cutoff),
                      hjust = 0, size = 3.5, color = "#A23B72", fontface = "bold")
  }

  return(p)
}

# --------------------------
# Create 4-panel combined figure (like Figure 3.2)
# --------------------------
create_combined_figure <- function(configs_list, omega_vec, x_label) {
  plots_list <- list()

  for (i in seq_along(configs_list)) {
    cfg <- configs_list[[i]]

    # Generate data
    res <- cascade_response(cfg$stages, cfg$k2, cfg$k3, cfg$k4, omega_vec)
    df <- data.frame(frequency = omega_vec,
                     magnitude = res$magnitude,
                     phase = res$phase,
                     config = cfg$desc,
                     stringsAsFactors = FALSE)

    # Compute metrics
    metrics <- compute_enhanced_metrics(df)

    # Create magnitude plot
    p_mag <- plot_mag_enhanced(df, metrics,
                               paste0("(", LETTERS[i], ") ", cfg$short),
                               x_label, show_annotations = TRUE)

    # Create phase plot
    p_phase <- plot_phase_enhanced(df, metrics, "", x_label, show_annotations = TRUE)

    # Combine magnitude and phase vertically
    plots_list[[i]] <- p_mag / p_phase + plot_layout(heights = c(1.5, 1))
  }

  # Arrange all 4 configs in 2x2 grid
  combined <- (plots_list[[1]] | plots_list[[2]]) / (plots_list[[3]] | plots_list[[4]])

  return(combined)
}

# --------------------------
# Main execution
# --------------------------
cat("Starting enhanced Bode plot generation...\n")
cat("Formatted for PLOS Computational Biology\n")

rad_s_dir <- make_dir(out_root, "rad_s")
Hz_dir <- make_dir(out_root, "Hz")
summ_dir <- make_dir(out_root, "summaries")

metrics_rad_list <- list()
metrics_Hz_list <- list()

# Generate individual plots for each config
for (i in seq_along(configs)) {
  cfg <- configs[[i]]
  cat(sprintf("Processing %s: %s\n", cfg$id, cfg$desc))

  # rad/s
  res_rad <- cascade_response(cfg$stages, cfg$k2, cfg$k3, cfg$k4, omega_rad)
  df_rad <- data.frame(frequency = omega_rad,
                       magnitude = res_rad$magnitude,
                       phase = res_rad$phase,
                       config = cfg$desc,
                       stringsAsFactors = FALSE)

  metrics_rad <- compute_enhanced_metrics(df_rad)
  metrics_rad_list[[i]] <- metrics_rad

  # Create enhanced plots
  p_mag_rad <- plot_mag_enhanced(df_rad, metrics_rad, cfg$desc, "Frequency (rad/s)")
  p_phase_rad <- plot_phase_enhanced(df_rad, metrics_rad, cfg$desc, "Frequency (rad/s)")
  p_full_rad <- p_mag_rad / p_phase_rad + plot_layout(heights = c(1.5, 1))

  # Save rad/s plots
  png_dir <- make_dir(rad_s_dir, "PNG")
  pdf_dir <- make_dir(rad_s_dir, "PDF")
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", cfg$id)

  ggsave(file.path(png_dir, paste0(safe_name, "_enhanced.png")),
         p_full_rad, width = 10, height = 7, dpi = 300)
  ggsave(file.path(pdf_dir, paste0(safe_name, "_enhanced.pdf")),
         p_full_rad, width = 10, height = 7)

  # Hz
  res_Hz <- cascade_response(cfg$stages, cfg$k2, cfg$k3, cfg$k4, omega_rad)
  df_Hz <- data.frame(frequency = freq_hz,
                      magnitude = res_Hz$magnitude,
                      phase = res_Hz$phase,
                      config = cfg$desc,
                      stringsAsFactors = FALSE)

  metrics_Hz <- compute_enhanced_metrics(df_Hz)
  metrics_Hz_list[[i]] <- metrics_Hz
}

# Generate combined 4-panel figure (Figure 3.2 style)
cat("\nGenerating combined 4-panel figure...\n")
combined_rad <- create_combined_figure(configs, omega_rad, "Frequency (rad/s)")
combined_Hz <- create_combined_figure(configs, freq_hz, "Frequency (Hz)")

# Save combined figures at PLOS-recommended resolution
ggsave(file.path(out_root, "Figure_3.2_combined_rad_s.png"),
       combined_rad, width = 14, height = 12, dpi = 300)
ggsave(file.path(out_root, "Figure_3.2_combined_rad_s.pdf"),
       combined_rad, width = 14, height = 12)
ggsave(file.path(out_root, "Figure_3.2_combined_Hz.png"),
       combined_Hz, width = 14, height = 12, dpi = 300)
ggsave(file.path(out_root, "Figure_3.2_combined_Hz.pdf"),
       combined_Hz, width = 14, height = 12)

# Save metrics with minimal decimal places
metrics_rad_df <- do.call(rbind, metrics_rad_list)
metrics_Hz_df <- do.call(rbind, metrics_Hz_list)

write.csv(metrics_rad_df, file.path(summ_dir, "enhanced_metrics_rad_s.csv"), row.names = FALSE)
write.csv(metrics_Hz_df, file.path(summ_dir, "enhanced_metrics_Hz.csv"), row.names = FALSE)

cat("\n=== Enhanced Bode plot generation complete ===\n")
cat(sprintf("Output directory: %s\n", out_root))
cat("Combined 4-panel figures saved as Figure_3.2_combined_*.png/pdf\n")
cat("Individual enhanced plots saved in rad_s/ and Hz/ subdirectories\n")
cat("Enhanced metrics (with minimal decimal places) saved in summaries/\n")
cat("\nFormatting applied:\n")
cat("  - CSV files: 2 decimal places for dB values, 1 for degrees, 3 significant figures for frequencies\n")
cat("  - Plots: Bold text for all axis labels, titles, and annotations\n")
cat("  - Resolution: 300 dpi (PLOS standard)\n")
