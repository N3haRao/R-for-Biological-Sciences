library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(filename, delimiter) {
  # Reading the CSV file as a dataframe
  data <- read.csv(filename, sep = delimiter)
  return(data)
}

#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
  # Calculating the proportion of variance explained by each PC
  pc_variance_explained <- (pca_results$sdev^2) / sum(pca_results$sdev^2)
  return(pc_variance_explained)
}

#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
  # Creating a tibble with PC names, variance explained, and cumulative variance explained
  variance_tibble <- tibble(
    principal_components = paste("PC", 1:length(pca_ve)),
    variance_explained = pca_ve,
    cumulative = cumsum(pca_ve)
  )
  return(variance_tibble)
}


#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
  # Reading the metadata
  metadata_df <- read.csv(metadata)
  
  # Creating a biplot using PC1 and PC2 as x and y, and labelling points with SixSubTypesClassification
  biplot_data <- data.frame(
    PC1 = pca_results$x[, 1],
    PC2 = pca_results$x[, 2],
    Label = metadata_df$SixSubTypesClassification
  )
  
  biplot <- ggplot(biplot_data, aes(x = PC1, y = PC2, label = Label)) +
    geom_point() +
    geom_text(size = 3, vjust = -0.5) +  # Add labels to points
    labs(x = "PC1", y = "PC2") +
    ggtitle("Biplot of PC1 vs. PC2 labeled by SixSubTypesClassification")
  
  return(biplot)
}

#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
  # Filtering the tibble to keep only significant probes based on the FDR threshold
  sig_probes <- diff_exp_tibble$probeid[diff_exp_tibble$padj < fdr_threshold]
  return(sig_probes)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids) {
  # Making sure the sig_ids are character strings
  sig_ids <- as.character(sig_ids)
  
  # Checking if all sig_ids are present in the row names of the intensity matrix
  if (!all(sig_ids %in% rownames(intensity))) {
    stop("Some sig_ids are not found in the intensity matrix.")
  }
  
  # Getting the row indices corresponding to sig_ids
  row_indices <- match(sig_ids, rownames(intensity))
  
  # Sub-setting the intensity matrix based on row indices and convert to a matrix
  de_intensity <- as.matrix(intensity[row_indices, , drop = FALSE])
  
  # Ensuring that the result is a matrix
  if (!is.matrix(de_intensity)) {
    stop("The result is not a matrix.")
  }
  
  return(de_intensity)
}


#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
  # Creating a heatmap using the specified palette
  heatmap(
    de_intensity,
    col = brewer.pal(num_colors, palette),
    main = "Heatmap of Normalized Intensity Values",
    xlab = "Samples",
    ylab = "Probes",
    scale = "none"  # We are not scaling the data for the heatmap
  )
}
