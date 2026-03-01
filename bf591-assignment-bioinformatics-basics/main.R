#!/usr/bin/Rscript
## Author: Neha Rao
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!require("biomaRt", quietly = TRUE)){
  BiocManager::install("biomaRt")
}
library(biomaRt)
library(tidyverse)

#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @examples 
#' `data <- load_expression('data/example_intensity_data_subset.csv')`
load_expression <- function(filepath) {
  # Read the CSV file and return it as a tibble
  data <- read.csv(filepath)
  return(as_tibble(data))
}

#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#' @details This is similar to the filters being implemented in BF528's project 1. 
#' We do not necessarily want to capture all parts of the assay in our analysis, so 
#' filters like this serve to reduce the noise and amount of data to examine.
#'
#' @examples `samples <- filter_15(data_tib)`
#' 
filter_15 <- function(tibble) {
  # Calculate the proportion of values greater than log2(15) for each row
  tibble <- tibble %>%
    rowwise() %>%
    mutate(good_genes = sum(c_across(starts_with("GSM")) > log2(15)) / ncol(tibble) >= 0.15)
  
  # Print the number of probes remaining in filtered dataset
  print(paste0("Probes remaining in filtered dataset: ", sum(tibble$good_genes)))
  
  # Filter the tibble based on the criteria
  filtered_tibble <- tibble %>%
    filter(good_genes) %>%
    select(probe)
  
  return(filtered_tibble)
}


#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
affy_to_hgnc <- function(affy_vector) {
  ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  response <- getBM(
    mart = ensembl,
    attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
    filter = "affy_hg_u133_plus_2",
    values = affy_vector
  )
  return(response)
}



#### ggplot ####

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
reduce_data <- function(expr_tibble, names_ids, good_genes, bad_genes) {
  expr_tibble <- expr_tibble %>%
    # Convert the probe column to character type
    mutate(probe = as.character(probe)) %>%
    # Match probe IDs with HGNC symbols
    left_join(names_ids, by = c("probe" = "affy_hg_u133_plus_2")) %>%
    # Ensure hgnc_symbol is character type
    mutate(hgnc_symbol = as.character(hgnc_symbol)) %>%
    # Categorize genes as 'good,' 'bad,' or NA
    mutate(
      gene_set = case_when(
        !is.na(hgnc_symbol) & hgnc_symbol %in% good_genes ~ "good",
        !is.na(hgnc_symbol) & hgnc_symbol %in% bad_genes ~ "bad",
        TRUE ~ NA_character_
      )
    ) %>%
    # Filter rows with a valid gene_set value (good or bad)
    filter(!is.na(gene_set)) %>%
    # Reorder columns to match the expected order
    select(probe, hgnc_symbol, gene_set, everything())
  
  return(expr_tibble)
}








#' Convert the tibble from wide to long format.
#'
#' @param tibble A tibble of expression data in wide format with information about
#' good and bad genes, gene names, sample names, and expression values.
#'
#' @return A tibble in long format containing the same information.
#'
#' @details This function's primary objective is to reformat the tibble from a 
#' wide format, where there are separate columns for sample name and expression value, 
#' to a long format, where sample names are stored in a column named 'sample' 
#' and their values stored in another column named 'value'.
#'
#' Convert the tibble from wide to long format.
convert_to_long <- function(tibble) {
  # Specify the columns to pivot from wide to long format
  cols_to_pivot <- names(tibble)[grep("^GSM", names(tibble))]
  
  # Convert wide format to long format
  long_tibble <- pivot_longer(tibble, cols = all_of(cols_to_pivot), names_to = "sample", values_to = "value")
  return(long_tibble)
}

#' Plot gene expression values as a boxplot using ggplot2.
#'
#' @param long_format_tibble A tibble in long format containing gene expression data
#' with columns 'hgnc_symbol', 'sample', and 'value'.
#'
#' @return A boxplot visualizing the distribution of gene expression values.
#'
#' @details This function creates a boxplot using ggplot2 to visualize the distribution
#' of gene expression values for the 'good' and 'bad' gene sets.
#'
#' Plot gene expression values as a boxplot using ggplot2.
plot_boxplot <- function(long_format_tibble) {
  # Create a ggplot boxplot
  plot <- ggplot(long_format_tibble, aes(x = gene_set, y = value, fill = gene_set)) +
    geom_boxplot(outlier.shape = NA) +
    labs(
      title = "Boxplot of Gene Expression Levels",
      x = "Gene Set",
      y = "Expression Level (Log2)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 16, face = "bold")
    )
  
  # Display the plot
  print(plot)
}
