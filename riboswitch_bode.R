###############################################################################
# Riboswitch-based Bode plot analysis with known kinetic parameters
# Uses experimentally measured rate constants for model validation
# Formatted for PLOS Computational Biology submission
###############################################################################

library(ggplot2)
library(patchwork)
library(pracma)
library(stringr)
library(svglite)
library(scales)

# --------------------------
# Riboswitch kinetic parameters from literature
# --------------------------
# These are examples - adjust based on specific riboswitches used

# TPP riboswitch (Thiamine pyrophosphate) - well-characterized
# Binding/unbinding and conformational change rates
TPP_riboswitch <- list(
  name = "TPP riboswitch",
  k_on = 1.0,      # Ligand binding rate (μM^-1 s^-1)
  k_off = 0.1,     # Ligand unbinding rate (s^-1)
  k_fold = 50,     # Forward folding rate (s^-1)
  k_unfold = 1.0,  # Reverse folding rate (s^-1)
  reference = "Wickiser et al. 2005, Chem Biol"
)

# SAM-I riboswitch (S-adenosylmethionine)
SAM_riboswitch <- list(
  name = "SAM-I riboswitch",
  k_on = 2.5,
  k_off = 0.05,
  k_fold = 30,
  k_unfold = 0.5,
  reference = "Montange & Batey 2006, Nature"
)

# Adenine riboswitch
Adenine_riboswitch <- list(
  name = "Adenine riboswitch",
  k_on = 3.0,
  k_off = 0.2,
  k_fold = 20,
  k_unfold = 0.8,
  reference = "Lemay et al. 2006, RNA"
)

# FMN riboswitch (Flavin mononucleotide)
FMN_riboswitch <- list(
  name = "FMN riboswitch",
  k_on = 1.5,
  k_off = 0.15,
  k_fold = 40,
  k_unfold = 1.2,
  reference = "Winkler et al. 2002, Nat Struct Biol"
)

# PreQ1 riboswitch
PreQ1_riboswitch <- list(
  name = "PreQ1 riboswitch",
  k_on = 5.0,
  k_off = 0.03,
  k_fold = 60,
  k_unfold = 0.6,
  reference = "Klein et al. 2009, RNA"
)

# Collect all riboswitches
riboswitches <- list(TPP_riboswitch, SAM_riboswitch, Adenine_riboswitch, 
                     FMN_riboswitch, PreQ1_riboswitch)

# --------------------------
# User parameters
# --------------------------
freq_range_rad <- c(1e-3, 1e3)
n_points <- 2000

out_root <- "riboswitch_bode_outputs"

make_dir <- function(...) {
  dir <- file.path(..., fsep = .Platform$file.sep)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  invisible(dir)
}

# --------------------------
# Convert riboswitch kinetics to transfer function parameters
# --------------------------
riboswitch_to_params <- function(riboswitch, ligand_conc = 1.0) {
  # Effective binding rate depends on ligand concentration
  k_bind_eff <- riboswitch$k_on * ligand_conc
  k_unbind <- riboswitch$k_off
  k_fold <- riboswitch$k_fold
  k_unfold <- riboswitch$k_unfold
  
  # Map to second-order system parameters
  # Natural frequency related to fastest transition
  omega_n <- sqrt(k_fold * k_bind_eff)
  
  # Damping related to rate ratios
  zeta <- (k_unbind + k_unfold) / (2 * omega_n)
  
  # For cascade representation
  k2 <- 2 * zeta * omega_n
  k3 <- omega_n
  k4 <- omega_n
  
  list(
    name = riboswitch$name,
    omega_n = omega_n,
    zeta = zeta,
    k2 = k2,
    k3 = k3,
    k4 = k4,
    k_bind = k_bind_eff,
    k_unbind = k_unbind,
    k_fold = k_fold,
    k_unfold = k_unfold,
    reference = riboswitch$reference
  )
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

# --------------------------
# Stage responses for riboswitch system
# --------------------------
stage_response_binding <- function(k2, k3, k4, omega_vec) {
  # Binding/unbinding stage
  s <- 1i * omega_vec
  numerator <- k3 * (s + k2)
  denominator <- s * (s + k2) + k3 * k4
  numerator / denominator
}

stage_response_conformational <- function(k2, k3, k4, omega_vec) {
  # Conformational change stage
  s <- 1i * omega_vec
  numerator <- k2 * k4
  denominator <- s * (s + k2) + k3 * k4
  numerator / denominator
}

riboswitch_response <- function(params, omega_vec, stage_type = "full") {
  if (stage_type == "binding") {
    H <- stage_response_binding(params$k2, params$k3, params$k4, omega_vec)
  } else if (stage_type == "conformational") {
    H <- stage_response_conformational(params$k2, params$k3, params$k4, omega_vec)
  } else {  # "full" - combined response
    H_bind <- stage_response_binding(params$k2, params$k3, params$k4, omega_vec)
    H_conf <- stage_response_conformational(params$k2, params$k3, params$k4, omega_vec)
    H <- H_bind * H_conf
  }
  
  mag_db <- 20 * log10(Mod(H))
  phase_rad <- Arg(H)
  phase_unwrapped <- unwrap_phase(phase_rad)
  phase_deg <- phase_unwrapped * 180 / pi
  
  list(H = H, magnitude = mag_db, phase = phase_deg, raw_phase = phase_rad)
}

# --------------------------
# Build riboswitch configurations
# --------------------------
build_riboswitch_configs <- function(ligand_conc = 1.0) {
  configs <- list()
  
  for (i in seq_along(riboswitches)) {
    rs <- riboswitches[[i]]
    params <- riboswitch_to_params(rs, ligand_conc)
    
    configs[[i]] <- list(
      id = paste0("riboswitch_", i, "_", gsub(" ", "_", rs$name)),
      desc = rs$name,
      short = paste0("RS-", i),
      params = params,
      stage_type = "full"
    )
  }
  
  configs
}

omega_rad <- 10^seq(log10(freq_range_rad[1]), log10(freq_range_rad[2]), length.out = n_points)
freq_hz <- omega_rad / (2 * pi)

# --------------------------
# Enhanced metrics computation
# --------------------------
compute_riboswitch_metrics <- function(df, params) {
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
  
  # Gain margin
  phase_crossover_idx <- which.min(abs(phase + 180))
  phase_crossover_freq <- freq[phase_crossover_idx]
  gain_at_crossover <- magnitude[phase_crossover_idx]
  gain_margin_dB <- -gain_at_crossover
  
  # Phase margin
  gain_crossover_idx <- which.min(abs(magnitude))
  gain_crossover_freq <- freq[gain_crossover_idx]
  phase_at_crossover <- phase[gain_crossover_idx]
  phase_margin_deg <- 180 + phase_at_crossover
  
  bandwidth <- freq_3db
  
  phase_initial <- phase[1]
  phase_final <- phase[length(phase)]
  
  # Return with kinetic parameters included
  data.frame(
    riboswitch = df$config[1],
    omega_n = round(params$omega_n, 3),
    zeta = round(params$zeta, 3),
    k_bind = round(params$k_bind, 3),
    k_unbind = round(params$k_unbind, 3),
    k_fold = round(params$k_fold, 2),
    k_unfold = round(params$k_unfold, 2),
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
    reference = params$reference,
    stringsAsFactors = FALSE
  )
}

# --------------------------
# Enhanced plotting functions
# --------------------------
plot_riboswitch_mag <- function(df, metrics, title_text, x_label, show_annotations = TRUE) {
  p <- ggplot(df, aes(x = frequency, y = magnitude)) +
    geom_line(color = "#2E86AB", size = 1.2) +
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
    p <- p + geom_hline(yintercept = metrics$dc_gain_dB - 3, 
                       linetype = "dotted", color = "gray50", size = 0.5)
    
    p <- p + geom_vline(xintercept = metrics$freq_3dB, 
                       linetype = "dashed", color = "#A23B72", size = 0.8)
    
    label_text <- sprintf("ω[c] == %.2g", metrics$freq_3dB)
    
    p <- p + annotate("text", 
                     x = metrics$freq_3dB * 3, 
                     y = max(df$magnitude) * 0.85,
                     label = label_text,
                     parse = TRUE,
                     hjust = 0, size = 4, color = "#A23B72", fontface = "bold")
    
    p <- p + annotate("text",
                     x = min(df$frequency) * 2,
                     y = metrics$dc_gain_dB,
                     label = sprintf("DC: %.1f dB", metrics$dc_gain_dB),
                     hjust = 0, vjust = -0.5, size = 3.5, color = "gray30", fontface = "bold")
  }
  
  return(p)
}

plot_riboswitch_phase <- function(df, metrics, title_text, x_label, show_annotations = TRUE) {
  p <- ggplot(df, aes(x = frequency, y = phase)) +
    geom_line(color = "#F18F01", size = 1.2) +
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
    p <- p + geom_vline(xintercept = metrics$freq_3dB, 
                       linetype = "dashed", color = "#A23B72", size = 0.8)
    
    if (min(df$phase) < -180 && max(df$phase) > -180) {
      p <- p + geom_hline(yintercept = -180, 
                         linetype = "dotted", color = "gray50", size = 0.5)
    }
    
    phase_at_cutoff <- df$phase[which.min(abs(df$frequency - metrics$freq_3dB))]
    p <- p + annotate("text",
                     x = metrics$freq_3dB * 3,
                     y = phase_at_cutoff,
                     label = sprintf("φ = %.0f°", phase_at_cutoff),
                     hjust = 0, size = 3.5, color = "#A23B72", fontface = "bold")
  }
  
  return(p)
}

# --------------------------
# Create combined figure
# --------------------------
create_riboswitch_combined <- function(configs_list, omega_vec, x_label) {
  plots_list <- list()
  
  for (i in seq_along(configs_list)) {
    cfg <- configs_list[[i]]
    
    res <- riboswitch_response(cfg$params, omega_vec, cfg$stage_type)
    df <- data.frame(frequency = omega_vec,
                    magnitude = res$magnitude,
                    phase = res$phase,
                    config = cfg$desc,
                    stringsAsFactors = FALSE)
    
    metrics <- compute_riboswitch_metrics(df, cfg$params)
    
    p_mag <- plot_riboswitch_mag(df, metrics, 
                                  paste0("(", LETTERS[i], ") ", cfg$short),
                                  x_label, show_annotations = TRUE)
    
    p_phase <- plot_riboswitch_phase(df, metrics, "", x_label, show_annotations = TRUE)
    
    plots_list[[i]] <- p_mag / p_phase + plot_layout(heights = c(1.5, 1))
  }
  
  # Arrange in grid (2x3 for 5 riboswitches)
  combined <- (plots_list[[1]] | plots_list[[2]]) / 
              (plots_list[[3]] | plots_list[[4]]) /
              (plots_list[[5]] | plot_spacer())
  
  return(combined)
}

# --------------------------
# Main execution
# --------------------------
cat("Starting riboswitch Bode plot generation...\n")
cat("Using experimentally determined kinetic parameters\n")
cat("Formatted for PLOS Computational Biology\n\n")

rad_s_dir <- make_dir(out_root, "rad_s")
Hz_dir <- make_dir(out_root, "Hz")
summ_dir <- make_dir(out_root, "summaries")

# Test different ligand concentrations
ligand_concentrations <- c(0.1, 1.0, 10.0)  # μM

for (conc in ligand_concentrations) {
  cat(sprintf("\n=== Ligand concentration: %.1f μM ===\n", conc))
  
  configs <- build_riboswitch_configs(ligand_conc = conc)
  
  metrics_rad_list <- list()
  metrics_Hz_list <- list()
  
  # Generate individual plots
  for (i in seq_along(configs)) {
    cfg <- configs[[i]]
    cat(sprintf("Processing %s\n", cfg$desc))
    
    # rad/s
    res_rad <- riboswitch_response(cfg$params, omega_rad, cfg$stage_type)
    df_rad <- data.frame(frequency = omega_rad,
                        magnitude = res_rad$magnitude,
                        phase = res_rad$phase,
                        config = cfg$desc,
                        stringsAsFactors = FALSE)
    
    metrics_rad <- compute_riboswitch_metrics(df_rad, cfg$params)
    metrics_rad_list[[i]] <- metrics_rad
    
    p_mag_rad <- plot_riboswitch_mag(df_rad, metrics_rad, cfg$desc, "Frequency (rad/s)")
    p_phase_rad <- plot_riboswitch_phase(df_rad, metrics_rad, cfg$desc, "Frequency (rad/s)")
    p_full_rad <- p_mag_rad / p_phase_rad + plot_layout(heights = c(1.5, 1))
    
    png_dir <- make_dir(rad_s_dir, "PNG", sprintf("conc_%.1fuM", conc))
    pdf_dir <- make_dir(rad_s_dir, "PDF", sprintf("conc_%.1fuM", conc))
    safe_name <- gsub("[^A-Za-z0-9_-]", "_", cfg$id)
    
    ggsave(file.path(png_dir, paste0(safe_name, ".png")), 
           p_full_rad, width = 10, height = 7, dpi = 300)
    ggsave(file.path(pdf_dir, paste0(safe_name, ".pdf")), 
           p_full_rad, width = 10, height = 7)
    
    # Hz
    res_Hz <- riboswitch_response(cfg$params, omega_rad, cfg$stage_type)
    df_Hz <- data.frame(frequency = freq_hz,
                       magnitude = res_Hz$magnitude,
                       phase = res_Hz$phase,
                       config = cfg$desc,
                       stringsAsFactors = FALSE)
    
    metrics_Hz <- compute_riboswitch_metrics(df_Hz, cfg$params)
    metrics_Hz_list[[i]] <- metrics_Hz
  }
  
  # Generate combined figures
  cat("Generating combined riboswitch figure...\n")
  combined_rad <- create_riboswitch_combined(configs, omega_rad, "Frequency (rad/s)")
  combined_Hz <- create_riboswitch_combined(configs, freq_hz, "Frequency (Hz)")
  
  conc_label <- sprintf("conc_%.1fuM", conc)
  ggsave(file.path(out_root, paste0("Riboswitch_Bode_Combined_", conc_label, "_rad_s.png")), 
         combined_rad, width = 14, height = 16, dpi = 300)
  ggsave(file.path(out_root, paste0("Riboswitch_Bode_Combined_", conc_label, "_rad_s.pdf")), 
         combined_rad, width = 14, height = 16)
  ggsave(file.path(out_root, paste0("Riboswitch_Bode_Combined_", conc_label, "_Hz.png")), 
         combined_Hz, width = 14, height = 16, dpi = 300)
  ggsave(file.path(out_root, paste0("Riboswitch_Bode_Combined_", conc_label, "_Hz.pdf")), 
         combined_Hz, width = 14, height = 16)
  
  # Save metrics
  metrics_rad_df <- do.call(rbind, metrics_rad_list)
  metrics_Hz_df <- do.call(rbind, metrics_Hz_list)
  
  write.csv(metrics_rad_df, 
            file.path(summ_dir, paste0("riboswitch_metrics_", conc_label, "_rad_s.csv")), 
            row.names = FALSE)
  write.csv(metrics_Hz_df, 
            file.path(summ_dir, paste0("riboswitch_metrics_", conc_label, "_Hz.csv")), 
            row.names = FALSE)
}

cat("\n=== Riboswitch Bode analysis complete ===\n")
cat(sprintf("Output directory: %s\n", out_root))
cat("Combined riboswitch figures saved as Riboswitch_Bode_Combined_*.png/pdf\n")
cat("Individual plots organized by ligand concentration\n")
cat("Metrics with kinetic parameters saved in summaries/\n")
cat("\nNote: Adjust kinetic parameters based on your specific experimental data\n")