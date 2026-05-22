#!/usr/bin/env Rscript
# =====================================================
# SGB*_Result Folder Batch Analysis Pipeline
# Specialized Version - Adapted for SGB*_Result Folder Structure
# Enhanced with evolutionary parameter scenario analysis
# =====================================================

suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org", quiet = TRUE)
    library(optparse)
  }
})

# ----------------------------- Parameter Parsing -----------------------------
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-i", "--input-dir"),
                type = "character",
                help = "Parent directory path containing multiple SGB*_Result folders",
                metavar = "PATH"),
    make_option(c("-o", "--output-root"),
                type = "character",
                default = NULL,
                help = "Output root directory [default: input_dir/Strain_Analysis_Results/]",
                metavar = "PATH"),
    make_option(c("--pattern"),
                type = "character",
                default = "^SGB.*_Result$",
                help = "SGB folder name pattern [default: %default]",
                metavar = "PATTERN"),
    make_option(c("--skip-existing"),
                action = "store_true",
                default = FALSE,
                help = "Skip SGBs with existing analysis results"),
    make_option(c("--debug"),
                action = "store_true",
                default = FALSE,
                help = "Debug mode, verbose output"),
    make_option(c("--genome-size"),
                type = "integer",
                default = 5000000,
                help = "Genome size (bp) for divergence time calculation [default: %default]",
                metavar = "INT"),
    make_option(c("--mutation-rate-low"),
                type = "double",
                default = 1e-6,
                help = "Low mutation rate per generation [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--mutation-rate-medium"),
                type = "double",
                default = 5e-6,
                help = "Medium mutation rate per generation [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--mutation-rate-high"),
                type = "double",
                default = 1e-5,
                help = "High mutation rate per generation [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--generation-time-low"),
                type = "double",
                default = 1,
                help = "Low generation time (hours) [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--generation-time-medium"),
                type = "double",
                default = 2,
                help = "Medium generation time (hours) [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--generation-time-high"),
                type = "double",
                default = 4,
                help = "High generation time (hours) [default: %default]",
                metavar = "FLOAT"),
    make_option(c("--scenarios"),
                type = "character",
                default = "low,medium,high",
                help = "Evolutionary scenarios to run [default: %default]",
                metavar = "SCENARIOS"),
    make_option(c("--custom-params"),
                type = "character",
                default = NULL,
                help = "Custom parameters file (CSV with columns: scenario,mutation_rate,generation_time)",
                metavar = "FILE")
  )
  
  parser <- OptionParser(
    option_list = option_list,
    usage = "%prog [options]",
    description = "Batch processing pipeline for SGB*_Result folder strain evolution analysis with parameter scenarios"
  )
  
  opt <- parse_args(parser)
  
  # Validate required parameters
  if (is.null(opt$`input-dir`)) {
    print_help(parser)
    stop("ERROR: Must specify input directory containing SGB*_Result folders (-i/--input-dir)")
  }
  
  # Set default output directory
  if (is.null(opt$`output-root`)) {
    opt$`output-root` <- file.path(opt$`input-dir`, "Strain_Analysis_Results")
  }
  
  # Parse scenarios
  opt$scenarios <- unlist(strsplit(opt$scenarios, ","))
  
  # Parse custom parameters if provided
  if (!is.null(opt$`custom-params`) && file.exists(opt$`custom-params`)) {
    opt$custom_params <- read.csv(opt$`custom-params`, stringsAsFactors = FALSE)
  } else {
    opt$custom_params <- NULL
  }
  
  return(opt)
}

# ----------------------------- SGB Folder Validation -----------------------------
validate_sgb_folder <- function(sgb_path) {
  folder_name <- basename(sgb_path)
  
  # Extract SGB number (extract SGB4563 from SGB4563_Result)
  sgb_number <- gsub("_Result$", "", folder_name)
  
  # Check required files
  required_files <- c(
    "StrainCells.txt",
    paste0(sgb_number, "_SNPpd.txt")
  )
  
  missing_files <- c()
  for (f in required_files) {
    if (!file.exists(file.path(sgb_path, f))) {
      missing_files <- c(missing_files, f)
    }
  }
  
  if (length(missing_files) > 0) {
    return(list(
      valid = FALSE,
      sgb_name = sgb_number,
      error = paste("Missing files:", paste(missing_files, collapse = ", "))
    ))
  }
  
  # Check if files are empty
  for (f in required_files) {
    file_path <- file.path(sgb_path, f)
    if (file.size(file_path) == 0) {
      return(list(
        valid = FALSE,
        sgb_name = sgb_number,
        error = paste("Empty file:", f)
      ))
    }
  }
  
  return(list(
    valid = TRUE,
    sgb_name = sgb_number,
    sgb_number = sgb_number,
    folder_name = folder_name
  ))
}

# ----------------------------- Basic Statistical Analysis -----------------------------
run_basic_analysis <- function(sgb_name, sgb_output_dir) {
  setwd(sgb_output_dir)
  
  cells <- read.table("StrainCells.txt", header = TRUE, sep = "\t")
  snp_data <- read.table(paste0(sgb_name, "_SNPpd.txt"), 
                        header = TRUE, sep = "\t", check.names = FALSE)
  
  basic_stats <- data.frame(
    Metric = c("SGB_name", "Total_cell_count", "Strain_count", "SNP_site_count"),
    Value = c(sgb_name, nrow(cells), length(unique(cells$Cluster)), nrow(snp_data))
  )
  
  write.table(basic_stats, "basic_statistics.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(list(
    cells = cells,
    snp_data = snp_data,
    stats = basic_stats
  ))
}

# ----------------------------- Genotype Extraction -----------------------------
extract_genotypes <- function(sgb_name, cells, snp_data) {
  # Ensure dplyr package is available
  if (!require("dplyr", quietly = TRUE)) {
    install.packages("dplyr", repos = "https://cloud.r-project.org", quiet = TRUE, quiet = TRUE)
    library(dplyr)
  }
  
  # Get strain list (exclude Mist*)
  strains <- unique(cells[["Cluster"]])
  strains <- strains[!grepl("^Mist", strains)]
  
  if (length(strains) == 0) {
    cat("  Warning: No valid strains found (all start with 'Mist')\n")
    return(NULL)
  }
  
  cat("  Valid strains:", paste(strains, collapse = ", "), "\n")
  
  # Get cell columns
  cell_cols <- grep("Cell", colnames(snp_data), value = TRUE)
  
  # Extract genotypes for each strain
  for (strain in strains) {
    cat("  Processing strain:", strain, "\n")
    
    # Get cells for this strain
    strain_cells <- cells[["Cell"]][cells[["Cluster"]] == strain]
    
    # Find corresponding columns
    strain_cols <- character()
    for (cell in strain_cells) {
      matched <- grep(cell, cell_cols, value = TRUE)
      strain_cols <- c(strain_cols, matched)
    }
    
    if (length(strain_cols) == 0) {
      cat("    Warning: No matching cell columns found\n")
      next
    }
    
    cat("    Found", length(strain_cols), "cells\n")
    
    # Extract data
    strain_snps <- snp_data[, c("contig_list", "location_on_contig_list", 
                               "reference_list", "alter_list", strain_cols)]
    
    # Calculate allele frequencies
    freq_result <- apply(strain_snps[, 5:ncol(strain_snps)], 1, function(row) {
      vals <- as.numeric(row)
      alt_count <- sum(vals == 1, na.rm = TRUE)
      total_count <- sum(vals %in% c(-1, 0, 1), na.rm = TRUE)
      if (total_count > 0) alt_count / total_count else 0
    })
    
    # Dynamic threshold
    mean_freq <- mean(freq_result, na.rm = TRUE)
    cat("    Average frequency:", round(mean_freq, 4), "\n")
    
    threshold <- if (mean_freq < 0.05) 0.1 else if (mean_freq < 0.2) 0.3 else 0.5
    cat("    Using threshold: >", threshold, "\n")
    
    # Create result dataframe
    result <- data.frame(
      contig = strain_snps[, "contig_list"],
      position = strain_snps[, "location_on_contig_list"],
      ref = strain_snps[, "reference_list"],
      alt = strain_snps[, "alter_list"],
      frequency = round(freq_result, 4),
      binary = as.integer(freq_result > threshold)
    )
    
    # Save results
    output_file <- paste0("strain_", strain, "_genotype.tsv")
    write.table(result, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("    Saved to:", output_file, "\n")
    
    # Save binary version
    binary_file <- paste0("strain_", strain, "_binary.tsv")
    write.table(result[, "binary"], binary_file, sep = "\t",
                row.names = FALSE, col.names = FALSE)
    cat("    Binary version:", binary_file, "\n")
  }
  
  return(strains)
}

# ----------------------------- Key SNP Analysis -----------------------------
analyze_key_snps <- function(strains) {
  strain_files <- list.files(pattern = "strain_[0-9]+_genotype[.]tsv$")
  
  if (length(strain_files) == 0) {
    cat("  Warning: No strain genotype files found\n")
    return(NULL)
  }
  
  # Read frequency data
  freq_data <- list()
  for (file in strain_files) {
    strain <- gsub("strain_|_genotype[.]tsv", "", file)
    data <- read.table(file, header = TRUE, sep = "\t")
    freq_data[[strain]] <- data$frequency
  }
  
  # Create frequency matrix
  freq_matrix <- do.call(cbind, freq_data)
  
  # Strain-specific SNPs
  cat("  Identifying strain-specific SNPs...\n")
  for (strain in strains) {
    other_strains <- setdiff(strains, strain)
    
    if (length(other_strains) > 0) {
      strain_vals <- freq_matrix[, strain]
      other_vals <- freq_matrix[, other_strains, drop = FALSE]
      other_max <- apply(other_vals, 1, max)
      specific_idx <- which(strain_vals > 0.6 & other_max < 0.4)
    } else {
      specific_idx <- which(freq_matrix[, strain] > 0.6)
    }
    
    if (length(specific_idx) > 0) {
      sample_data <- read.table(strain_files[1], header = TRUE, sep = "\t")
      specific_snps <- sample_data[specific_idx, 1:4]
      specific_snps$frequency <- sprintf("%.3f", freq_matrix[specific_idx, strain])
      specific_snps$strain <- strain
      
      output_file <- paste0("strain_", strain, "_specific_snps.tsv")
      write.table(specific_snps, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("    Strain", strain, ":", length(specific_idx), "specific SNPs ->", output_file, "\n")
    } else {
      cat("    Strain", strain, ": 0 specific SNPs\n")
    }
  }
  
  # Highly differentiated SNPs
  if (length(strains) > 1) {
    sd_scores <- apply(freq_matrix, 1, sd)
    top_n <- min(20, length(sd_scores))
    top_idx <- order(sd_scores, decreasing = TRUE)[1:top_n]
    
    if (length(top_idx) > 0) {
      sample_data <- read.table(strain_files[1], header = TRUE, sep = "\t")
      diff_snps <- sample_data[top_idx, 1:4]
      diff_snps <- cbind(diff_snps, freq_matrix[top_idx, , drop = FALSE])
      diff_snps$sd_score <- sprintf("%.4f", sd_scores[top_idx])
      
      write.table(diff_snps, "top_differentiated_snps.tsv", 
                  sep = "\t", row.names = FALSE, quote = FALSE)
      cat("    Top", top_n, "differentiated SNPs saved to: top_differentiated_snps.tsv\n")
    }
  }
  
  return(TRUE)
}

# ----------------------------- Phylogenetic Analysis -----------------------------
build_phylogeny <- function() {
  binary_files <- list.files(pattern = "strain_[0-9]+_binary\\.tsv$")
  
  if (length(binary_files) < 2) {
    cat("  Warning: Need at least 2 strains for phylogeny\n")
    return(NULL)
  }
  
  # Ensure ape package is available
  if (!require("ape", quietly = TRUE)) {
    install.packages("ape", repos = "https://cloud.r-project.org", quiet = TRUE, quiet = TRUE)
    library(ape)
  }
  
  strains <- gsub("strain_|_binary\\.tsv", "", binary_files)
  cat("  Building phylogeny for strains:", paste(strains, collapse = ", "), "\n")
  
  # Read binary data
  binary_list <- list()
  for (file in binary_files) {
    strain <- gsub("strain_|_binary\\.tsv", "", file)
    binary_list[[strain]] <- read.table(file, header = FALSE)$V1
  }
  
  # Create matrix
  n_snps <- length(binary_list[[1]])
  binary_matrix <- matrix(0, nrow = length(strains), ncol = n_snps)
  rownames(binary_matrix) <- strains
  
  for (i in 1:length(strains)) {
    binary_matrix[i, ] <- binary_list[[strains[i]]]
  }
  
  # Calculate distance
  dist_matrix <- dist(binary_matrix, method = "manhattan") / n_snps
  write.csv(as.matrix(dist_matrix), "distance_matrix.csv", quote = FALSE)
  cat("    Distance matrix saved to: distance_matrix.csv\n")
  
  # Build NJ tree (if at least 3 strains)
  if (length(strains) >= 3) {
    nj_tree <- nj(as.matrix(dist_matrix))
    
    # Correct negative branch lengths
    if (any(nj_tree$edge.length < 0)) {
      nj_tree$edge.length[nj_tree$edge.length < 0] <- 0.001
    }
    
    # Save tree
    write.tree(nj_tree, "phylogeny_tree.nwk")
    cat("    Phylogenetic tree saved to: phylogeny_tree.nwk\n")
    
    # Plot tree
    png("phylogeny_tree.png", width = 800, height = 600, res = 150)
    plot(nj_tree, main = "Strain Phylogenetic Tree", 
         cex = 1.5, label.offset = 0.02, edge.width = 2)
    add.scale.bar()
    dev.off()
    cat("    Tree plot saved to: phylogeny_tree.png\n")
  } else {
    cat("    Less than 3 strains, skipping tree construction\n")
  }
  
  return(dist_matrix)
}

# ----------------------------- Parameter Scenarios -----------------------------
define_parameter_scenarios <- function(opt) {
  # Define default scenarios
  default_scenarios <- list(
    low = list(
      mutation_rate = opt$`mutation-rate-low`,
      generation_time = opt$`generation-time-low`
    ),
    medium = list(
      mutation_rate = opt$`mutation-rate-medium`,
      generation_time = opt$`generation-time-medium`
    ),
    high = list(
      mutation_rate = opt$`mutation-rate-high`,
      generation_time = opt$`generation-time-high`
    )
  )
  
  # Add custom scenarios if provided
  if (!is.null(opt$custom_params)) {
    custom_scenarios <- list()
    for (i in 1:nrow(opt$custom_params)) {
      scenario_name <- opt$custom_params$scenario[i]
      custom_scenarios[[scenario_name]] <- list(
        mutation_rate = opt$custom_params$mutation_rate[i],
        generation_time = opt$custom_params$generation_time[i]
      )
    }
    return(custom_scenarios)
  }
  
  # Filter to only requested scenarios
  scenarios <- list()
  for (scenario in opt$scenarios) {
    if (scenario %in% names(default_scenarios)) {
      scenarios[[scenario]] <- default_scenarios[[scenario]]
    }
  }
  
  return(scenarios)
}

# ----------------------------- Divergence Time Analysis with Scenarios -----------------------------
estimate_divergence_time_scenarios <- function(sgb_name, scenarios, genome_size) {
  if (!file.exists("distance_matrix.csv")) {
    cat("  Warning: Distance matrix not found, skipping divergence time analysis\n")
    return(NULL)
  }
  
  dist_data <- read.csv("distance_matrix.csv", row.names = 1, check.names = FALSE)
  strains <- rownames(dist_data)
  
  all_results <- list()
  
  # Process each parameter scenario
  for (scenario_name in names(scenarios)) {
    cat("  Calculating divergence times for scenario:", scenario_name, "\n")
    
    params <- scenarios[[scenario_name]]
    mutation_rate <- params$mutation_rate
    generation_time <- params$generation_time
    
    scenario_results <- list()
    
    for (i in 1:(nrow(dist_data) - 1)) {
      for (j in (i + 1):nrow(dist_data)) {
        strain1 <- strains[i]
        strain2 <- strains[j]
        genetic_distance <- dist_data[i, j]
        
        snp_count <- genetic_distance * genome_size
        generations <- snp_count / (2 * mutation_rate * genome_size)
        hours <- generations * generation_time
        days <- hours / 24
        years <- days / 365.25
        
        result <- data.frame(
          Strain1 = strain1,
          Strain2 = strain2,
          Scenario = scenario_name,
          Mutation_Rate = sprintf("%.2e", mutation_rate),
          Generation_Time = sprintf("%.1f", generation_time),
          Genetic_Distance = sprintf("%.6f", genetic_distance),
          Estimated_SNPs = sprintf("%.0f", snp_count),
          Years = sprintf("%.3f", years),
          Days = sprintf("%.0f", days),
          Hours = sprintf("%.0f", hours),
          Generations = sprintf("%.0f", generations)
        )
        
        scenario_results[[paste(strain1, strain2, sep = "_")]] <- result
      }
    }
    
    if (length(scenario_results) > 0) {
      scenario_df <- do.call(rbind, scenario_results)
      all_results[[scenario_name]] <- scenario_df
      
      # Save individual scenario results
      scenario_file <- paste0("divergence_time_", scenario_name, ".tsv")
      write.table(scenario_df, scenario_file, sep = "\t", row.names = FALSE, quote = FALSE)
      cat("    Saved scenario results to:", scenario_file, "\n")
    }
  }
  
  # Combine all scenarios
  if (length(all_results) > 0) {
    combined_results <- do.call(rbind, all_results)
    write.table(combined_results, "divergence_time_all_scenarios.tsv",
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("    All scenario results saved to: divergence_time_all_scenarios.tsv\n")
    
    # Generate summary statistics
    generate_scenario_summary(combined_results)
    
    # Plot comparison
    plot_scenario_comparison(combined_results)
  }
  
  return(combined_results)
}

# ----------------------------- Scenario Summary -----------------------------
generate_scenario_summary <- function(combined_results) {
  if (nrow(combined_results) == 0) return(NULL)
  
  summary_stats <- combined_results %>%
    group_by(Scenario) %>%
    summarise(
      Pairs = n(),
      Mean_Years = mean(as.numeric(Years), na.rm = TRUE),
      Median_Years = median(as.numeric(Years), na.rm = TRUE),
      Min_Years = min(as.numeric(Years), na.rm = TRUE),
      Max_Years = max(as.numeric(Years), na.rm = TRUE),
      SD_Years = sd(as.numeric(Years), na.rm = TRUE),
      Mean_Days = mean(as.numeric(Days), na.rm = TRUE),
      Mean_Generations = mean(as.numeric(Generations), na.rm = TRUE)
    )
  
  write.table(summary_stats, "scenario_summary.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("    Scenario summary saved to: scenario_summary.tsv\n")
  
  return(summary_stats)
}

# ----------------------------- Scenario Comparison Plot -----------------------------
plot_scenario_comparison <- function(combined_results) {
  if (nrow(combined_results) == 0 || length(unique(combined_results$Scenario)) < 2) {
    return(NULL)
  }
  
  if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE)
    library(ggplot2)
  }
  
  # Convert to numeric for plotting
  plot_data <- combined_results
  plot_data$Years_numeric <- as.numeric(plot_data$Years)
  plot_data$Days_numeric <- as.numeric(plot_data$Days)
  plot_data$Generations_numeric <- as.numeric(plot_data$Generations)
  
  # Create unique strain pair identifier
  plot_data$Strain_Pair <- paste(plot_data$Strain1, plot_data$Strain2, sep = "_")
  
  # Plot 1: Divergence time by scenario (years)
  p1 <- ggplot(plot_data, aes(x = Scenario, y = Years_numeric, fill = Scenario)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
    labs(title = "Divergence Time Estimates by Scenario",
         subtitle = "Comparison across different evolutionary parameters",
         y = "Divergence Time (Years)",
         x = "Evolutionary Scenario") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # Plot 2: Facet by strain pair
  if (length(unique(plot_data$Strain_Pair)) <= 10) {  # Only if not too many pairs
    p2 <- ggplot(plot_data, aes(x = Scenario, y = Years_numeric, color = Scenario)) +
      geom_point(size = 3) +
      geom_line(aes(group = Strain_Pair), alpha = 0.3) +
      facet_wrap(~ Strain_Pair, scales = "free_y") +
      labs(title = "Divergence Time Estimates by Strain Pair",
           y = "Divergence Time (Years)",
           x = "Evolutionary Scenario") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("scenario_comparison_by_pair.png", p2, width = 12, height = 8, dpi = 300)
  }
  
  # Plot 3: Parameter sensitivity
  if (!is.null(plot_data$Mutation_Rate) && !is.null(plot_data$Generation_Time)) {
    plot_data$Mutation_Rate_numeric <- as.numeric(plot_data$Mutation_Rate)
    plot_data$Generation_Time_numeric <- as.numeric(plot_data$Generation_Time)
    
    p3 <- ggplot(plot_data, aes(x = Mutation_Rate_numeric, y = Years_numeric,
                               color = factor(Generation_Time_numeric))) +
      geom_point(size = 3, alpha = 0.7) +
      scale_x_log10() +
      scale_y_log10() +
      labs(title = "Parameter Sensitivity Analysis",
           x = "Mutation Rate (log scale)",
           y = "Divergence Time (Years, log scale)",
           color = "Generation Time (hours)") +
      theme_minimal()
    
    ggsave("parameter_sensitivity.png", p3, width = 10, height = 8, dpi = 300)
  }
  
  ggsave("scenario_comparison.png", p1, width = 10, height = 8, dpi = 300)
  cat("    Scenario comparison plots saved\n")
  
  # Save parameter table
  param_table <- unique(plot_data[, c("Scenario", "Mutation_Rate", "Generation_Time")])
  write.table(param_table, "evolutionary_parameters.tsv",
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# ----------------------------- Generate Report with Scenarios -----------------------------
generate_report <- function(sgb_name, scenarios) {
  report_file <- paste0(sgb_name, "_analysis_report.md")
  
  cat("# Strain Evolution Analysis Report\n\n", file = report_file)
  cat("**SGB Name:**", sgb_name, "\n", file = report_file, append = TRUE)
  cat("**Analysis Time:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", 
      file = report_file, append = TRUE)
  
  # Basic statistics
  if (file.exists("basic_statistics.tsv")) {
    cat("## Basic Statistics\n\n", file = report_file, append = TRUE)
    cat("```\n", file = report_file, append = TRUE)
    stats <- readLines("basic_statistics.tsv")
    cat(stats, sep = "\n", file = report_file, append = TRUE)
    cat("\n```\n\n", file = report_file, append = TRUE)
  }
  
  # Evolutionary parameter scenarios
  if (!is.null(scenarios) && length(scenarios) > 0) {
    cat("## Evolutionary Parameter Scenarios\n\n", file = report_file, append = TRUE)
    cat("| Scenario | Mutation Rate | Generation Time (hours) |\n", file = report_file, append = TRUE)
    cat("|----------|---------------|-------------------------|\n", file = report_file, append = TRUE)
    
    for (scenario_name in names(scenarios)) {
      params <- scenarios[[scenario_name]]
      cat(sprintf("| %s | %.2e | %.1f |\n", 
                  scenario_name, params$mutation_rate, params$generation_time),
          file = report_file, append = TRUE)
    }
    cat("\n", file = report_file, append = TRUE)
  }
  
  # Strain-specific SNPs
  specific_files <- list.files(pattern = "strain_[0-9]+_specific_snps\\.tsv$")
  if (length(specific_files) > 0) {
    cat("## Strain-Specific SNPs\n\n", file = report_file, append = TRUE)
    for (file in specific_files) {
      strain <- gsub("strain_|_specific_snps\\.tsv", "", file)
      count <- length(readLines(file)) - 1  # Subtract header
      if (count > 0) {
        cat("- **Strain", strain, ":**", count, "specific SNPs\n", file = report_file, append = TRUE)
      }
    }
    cat("\n", file = report_file, append = TRUE)
  }
  
  # Divergence time summary
  if (file.exists("scenario_summary.tsv")) {
    cat("## Divergence Time Summary\n\n", file = report_file, append = TRUE)
    cat("```\n", file = report_file, append = TRUE)
    summary_data <- readLines("scenario_summary.tsv")
    cat(summary_data, sep = "\n", file = report_file, append = TRUE)
    cat("\n```\n\n", file = report_file, append = TRUE)
  }
  
  # Distance matrix preview
  if (file.exists("distance_matrix.csv")) {
    cat("## Genetic Distance Matrix (Preview)\n\n", file = report_file, append = TRUE)
    cat("```\n", file = report_file, append = TRUE)
    dist_content <- readLines("distance_matrix.csv", n = 6)
    cat(dist_content, sep = "\n", file = report_file, append = TRUE)
    cat("...\n```\n\n", file = report_file, append = TRUE)
  }
  
  # Generated files list
  cat("## Generated Files\n\n", file = report_file, append = TRUE)
  cat("### Data Files\n\n```\n", file = report_file, append = TRUE)
  data_files <- list.files(pattern = "\\.(tsv|csv|txt)$")
  cat(data_files, sep = "\n", file = report_file, append = TRUE)
  cat("\n```\n\n### Tree Files\n\n```\n", file = report_file, append = TRUE)
  tree_files <- list.files(pattern = "\\.(nwk|png)$")
  cat(tree_files, sep = "\n", file = report_file, append = TRUE)
  cat("\n```\n\n### Visualization Files\n\n```\n", file = report_file, append = TRUE)
  plot_files <- list.files(pattern = "\\.(png|pdf|svg)$")
  cat(plot_files, sep = "\n", file = report_file, append = TRUE)
  cat("\n```\n\n", file = report_file, append = TRUE)
  
  # Analysis notes
  cat("## Analysis Notes\n\n", file = report_file, append = TRUE)
  cat("1. Genotype calling: allele frequency > dynamic threshold (0.1-0.5)\n", file = report_file, append = TRUE)
  cat("2. Specific SNP: frequency > 0.6 in one strain and < 0.4 in others\n", file = report_file, append = TRUE)
  cat("3. Distance: normalized Hamming distance\n", file = report_file, append = TRUE)
  cat("4. Phylogeny: Neighbor-Joining method\n", file = report_file, append = TRUE)
  cat("5. Multiple evolutionary scenarios tested with different parameters\n\n", file = report_file, append = TRUE)
  cat("---\n**Analysis Complete**\n", file = report_file, append = TRUE)
  
  # Convert to text format
  txt_report <- paste0(sgb_name, "_analysis_report.txt")
  file.copy(report_file, txt_report, overwrite = TRUE)
  
  return(txt_report)
}

# ----------------------------- Single SGB Analysis -----------------------------
analyze_single_sgb <- function(sgb_path, sgb_name, output_dir, opt) {
  start_time <- Sys.time()
  
  # Create SGB-specific output directory
  sgb_output_dir <- file.path(output_dir, sgb_name)
  
  # Check if to skip existing analysis
  if (opt$`skip-existing` && dir.exists(sgb_output_dir)) {
    result_files <- list.files(sgb_output_dir, pattern = "_analysis_report\\.(md|txt)$")
    if (length(result_files) > 0) {
      cat("Skipping already analyzed SGB:", sgb_name, "\n")
      return(list(
        success = TRUE,
        sgb_name = sgb_name,
        status = "skipped",
        output_dir = sgb_output_dir
      ))
    }
  }
  
  dir.create(sgb_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Setup logging
  log_file <- file.path(sgb_output_dir, paste0(sgb_name, "_analysis.log"))
  if (opt$debug) {
    sink(log_file, append = FALSE, split = TRUE)
  } else {
    sink(log_file, append = FALSE)
  }
  
  cat("========================================\n")
  cat("  Analyzing:", sgb_name, "\n")
  cat("  Time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("  Input:", sgb_path, "\n")
  cat("  Output:", sgb_output_dir, "\n")
  cat("  Genome size:", opt$`genome-size`, "\n")
  cat("  Scenarios:", paste(opt$scenarios, collapse = ", "), "\n")
  cat("========================================\n\n")
  
  tryCatch({
    # 1. Copy input files
    cat("[1/9] Copying input files...\n")
    file.copy(
      from = file.path(sgb_path, "StrainCells.txt"),
      to = file.path(sgb_output_dir, "StrainCells.txt"),
      overwrite = TRUE
    )
    file.copy(
      from = file.path(sgb_path, paste0(sgb_name, "_SNPpd.txt")),
      to = file.path(sgb_output_dir, paste0(sgb_name, "_SNPpd.txt")),
      overwrite = TRUE
    )
    
    # Switch to output directory
    original_wd <- getwd()
    setwd(sgb_output_dir)
    
    # 2. Define parameter scenarios
    cat("[2/9] Defining parameter scenarios...\n")
    scenarios <- define_parameter_scenarios(opt)
    cat("    Defined", length(scenarios), "scenario(s):", paste(names(scenarios), collapse = ", "), "\n")
    
    # 3. Basic analysis
    cat("[3/9] Basic analysis...\n")
    basic_result <- run_basic_analysis(sgb_name, sgb_output_dir)
    
    # 4. Genotype extraction
    cat("[4/9] Extracting genotypes...\n")
    strains <- extract_genotypes(sgb_name, basic_result$cells, basic_result$snp_data)
    
    if (is.null(strains) || length(strains) == 0) {
      cat("  No valid strains to analyze\n")
      setwd(original_wd)
      sink()
      return(list(
        success = FALSE,
        sgb_name = sgb_name,
        error = "No valid strains found"
      ))
    }
    
    # 5. Key SNP analysis
    cat("[5/9] Key SNP analysis...\n")
    analyze_key_snps(strains)
    
    # 6. Phylogenetic analysis
    cat("[6/9] Phylogenetic analysis...\n")
    build_phylogeny()
    
    # 7. Divergence time analysis with scenarios
    cat("[7/9] Divergence time estimation (multiple scenarios)...\n")
    divergence_results <- estimate_divergence_time_scenarios(sgb_name, scenarios, opt$`genome-size`)
    
    # 8. Generate report
    cat("[8/9] Generating report...\n")
    report_file <- generate_report(sgb_name, scenarios)
    
    # 9. Completion
    end_time <- Sys.time()
    elapsed <- difftime(end_time, start_time, units = "secs")
    
    cat("[9/9] Analysis completed!\n")
    cat("  Time elapsed:", round(elapsed, 1), "seconds\n")
    cat("  Output directory:", sgb_output_dir, "\n")
    cat("  Report:", report_file, "\n")
    cat("  Scenarios analyzed:", paste(names(scenarios), collapse = ", "), "\n")
    
    setwd(original_wd)
    sink()
    
    return(list(
      success = TRUE,
      sgb_name = sgb_name,
      output_dir = sgb_output_dir,
      elapsed = elapsed,
      strains = strains,
      scenarios = names(scenarios),
      status = "completed"
    ))
    
  }, error = function(e) {
    setwd(original_wd)
    sink()
    cat("Analysis failed for", sgb_name, ":", e$message, "\n")
    
    return(list(
      success = FALSE,
      sgb_name = sgb_name,
      error = e$message
    ))
  })
}

# ----------------------------- Batch Processing Main Function -----------------------------
main <- function() {
  # Parse arguments
  opt <- parse_arguments()
  
  cat("========================================\n")
  cat("  Multi-SGB Strain Analysis Pipeline\n")
  cat("  Enhanced with Parameter Scenarios\n")
  cat("========================================\n")
  cat("Input directory:", opt$`input-dir`, "\n")
  cat("Output directory:", opt$`output-root`, "\n")
  cat("SGB pattern:", opt$pattern, "\n")
  cat("Genome size:", opt$`genome-size`, "bp\n")
  cat("Scenarios:", paste(opt$scenarios, collapse = ", "), "\n")
  cat("Skip existing:", opt$`skip-existing`, "\n")
  cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("========================================\n\n")
  
  # Find SGB*_Result folders
  all_items <- list.files(opt$`input-dir`, full.names = TRUE)
  sgb_folders <- all_items[dir.exists(all_items) & grepl(opt$pattern, basename(all_items))]
  
  if (length(sgb_folders) == 0) {
    stop("No SGB folders found matching pattern '", opt$pattern, "' in ", opt$`input-dir`)
  }
  
  cat("Found", length(sgb_folders), "SGB folders:\n")
  for (folder in sgb_folders) {
    cat("  -", basename(folder), "\n")
  }
  cat("\n")
  
  # Create output root directory
  dir.create(opt$`output-root`, recursive = TRUE, showWarnings = FALSE)
  
  # Create master log
  master_log <- file.path(opt$`output-root`, "batch_analysis.log")
  sink(master_log, append = FALSE, split = TRUE)
  
  results <- list()
  valid_count <- 0
  
  # Process each SGB
  for (sgb_folder in sgb_folders) {
    folder_name <- basename(sgb_folder)
    
    cat("\n", rep("=", 50), "\n", sep = "")
    cat("Processing:", folder_name, "\n")
    cat(rep("=", 50), "\n\n", sep = "")
    
    # Validate SGB folder
    validation <- validate_sgb_folder(sgb_folder)
    if (!validation$valid) {
      cat("Skipping", folder_name, ":", validation$error, "\n")
      results[[folder_name]] <- list(
        success = FALSE,
        error = validation$error
      )
      next
    }
    
    sgb_name <- validation$sgb_name
    cat("SGB name:", sgb_name, "\n")
    
    # Execute analysis
    result <- analyze_single_sgb(
      sgb_path = sgb_folder,
      sgb_name = sgb_name,
      output_dir = opt$`output-root`,
      opt = opt
    )
    
    results[[sgb_name]] <- result
    
    if (result$success) {
      if (result$status == "skipped") {
        cat("✓ Skipped (already exists):", sgb_name, "\n")
      } else {
        cat("✓ Completed:", sgb_name, "(", round(result$elapsed, 1), "s)\n")
        cat("  Scenarios:", paste(result$scenarios, collapse = ", "), "\n")
        valid_count <- valid_count + 1
      }
    } else {
      cat("✗ Failed:", sgb_name, "(", result$error, ")\n")
    }
  }
  
  # Generate summary report
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("Batch Analysis Summary\n")
  cat(rep("=", 50), "\n\n", sep = "")
  
  successful <- sum(sapply(results, function(x) x$success && x$status != "skipped"))
  skipped <- sum(sapply(results, function(x) x$success && x$status == "skipped"))
  failed <- length(results) - successful - skipped
  
  cat("Total processed:", length(results), "SGB folders\n")
  cat("Successfully analyzed:", successful, "\n")
  cat("Skipped (already exists):", skipped, "\n")
  cat("Failed:", failed, "\n")
  cat("Output directory:", opt$`output-root`, "\n")
  cat("Evolutionary scenarios:", paste(opt$scenarios, collapse = ", "), "\n")
  cat("Completion time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  
  # Save summary results
  summary_data <- data.frame(
    SGB_Folder = names(results),
    SGB_Name = sapply(results, function(x) if (!is.null(x$sgb_name)) x$sgb_name else NA),
    Status = sapply(results, function(x) {
      if (x$success) {
        if (x$status == "skipped") "skipped" else "completed"
      } else {
        "failed"
      }
    }),
    Time_sec = sapply(results, function(x) if (x$success && x$status == "completed") round(x$elapsed, 1) else NA),
    Strains = sapply(results, function(x) if (x$success && x$status == "completed") length(x$strains) else NA),
    Scenarios = sapply(results, function(x) if (x$success && x$status == "completed") paste(x$scenarios, collapse = ";") else NA),
    Output_Dir = sapply(results, function(x) if (x$success) x$output_dir else NA),
    Error = sapply(results, function(x) if (!x$success) x$error else NA)
  )
  
  write.table(summary_data, file.path(opt$`output-root`, "batch_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\nDetailed results saved to: batch_summary.tsv\n")
  
  # Generate batch-level parameter summary
  generate_batch_parameter_summary(opt)
  
  sink()
  
  # Final output
  cat("\n========================================\n")
  cat("Analysis Pipeline Completed!\n")
  cat("========================================\n")
  cat("Summary:\n")
  cat("  Total: ", length(results), "\n")
  cat("  Analyzed: ", successful, "\n")
  cat("  Skipped: ", skipped, "\n")
  cat("  Failed: ", failed, "\n")
  cat("Output: ", opt$`output-root`, "\n")
  cat("Log file: ", master_log, "\n")
  cat("Evolutionary parameters:\n")
  for (scenario in opt$scenarios) {
    if (scenario == "low") {
      cat("  - Low: μ=", opt$`mutation-rate-low`, ", g=", opt$`generation-time-low`, "h\n", sep = "")
    } else if (scenario == "medium") {
      cat("  - Medium: μ=", opt$`mutation-rate-medium`, ", g=", opt$`generation-time-medium`, "h\n", sep = "")
    } else if (scenario == "high") {
      cat("  - High: μ=", opt$`mutation-rate-high`, ", g=", opt$`generation-time-high`, "h\n", sep = "")
    }
  }
  cat("========================================\n")
}

# ----------------------------- Batch Parameter Summary -----------------------------
generate_batch_parameter_summary <- function(opt) {
  param_summary <- data.frame(
    Parameter = c("Genome_size", "Mutation_rate_low", "Mutation_rate_medium", "Mutation_rate_high",
                  "Generation_time_low", "Generation_time_medium", "Generation_time_high",
                  "Scenarios", "Custom_params"),
    Value = c(
      opt$`genome-size`,
      opt$`mutation-rate-low`,
      opt$`mutation-rate-medium`,
      opt$`mutation-rate-high`,
      opt$`generation-time-low`,
      opt$`generation-time-medium`,
      opt$`generation-time-high`,
      paste(opt$scenarios, collapse = ","),
      ifelse(is.null(opt$`custom-params`), "None", opt$`custom-params`)
    )
  )
  
  write.table(param_summary, file.path(opt$`output-root`, "parameter_summary.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Parameter summary saved to: parameter_summary.tsv\n")
}

# ----------------------------- Execute Main Function -----------------------------
if (!interactive()) {
  main()
}