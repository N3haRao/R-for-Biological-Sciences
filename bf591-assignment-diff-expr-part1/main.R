library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')
library('forcats')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile_csv (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(csv_path, metafile_csv, selected_times) {
  # Checking if input files exist and are readable
  if (!file.exists(csv_path) || !file.exists(metafile_csv)) {
    stop("Input files do not exist or are not readable.")
  }
  
  # Debug01: Printing the file paths
  cat("CSV Path: ", csv_path, "\n")
  cat("Metafile Path: ", metafile_csv, "\n")
  
  # Reading counts matrix
  counts <- as.matrix(read.table(csv_path, sep = "\t", header = TRUE, row.names = 1))
  
  # Debug02: Printing the dimensions of the counts matrix
  cat("Counts matrix dimensions: ", dim(counts), "\n")
  
  # Reading sample metadata
  sample_metadata <- read.csv(metafile_csv)
  
  # Debug03: Printing the dimensions of the sample metadata
  cat("Sample metadata dimensions: ", dim(sample_metadata), "\n")
  
  # Ensuring that the sample metadata contains 'samplename' and 'timepoint' columns
  if (!all(c("samplename", "timepoint") %in% colnames(sample_metadata))) {
    stop("Sample metadata must contain 'samplename' and 'timepoint' columns.")
  }
  
  # Converting the timepoint column to a factor
  sample_metadata$timepoint <- factor(sample_metadata$timepoint)
  
  # Setting the reference level for the timepoint factor to 'vP0'
  sample_metadata$timepoint <- relevel(sample_metadata$timepoint, ref = "vP0")
  
  # Filtering the sample metadata to match the selected times
  sample_metadata <- sample_metadata %>%
    filter(timepoint %in% selected_times)
  
  # Debug04: Printing the dimensions of the filtered sample metadata
  cat("Sample metadata dimensions after filtering: ", dim(sample_metadata), "\n")
  
  # Subsetting the counts matrix based on selected timepoints
  selected_samples <- sample_metadata$samplename
  counts_subset <- counts[, selected_samples]
  
  # Debug05: Printing the dimensions of the subsetted counts matrix
  cat("Counts matrix dimensions after subset: ", dim(counts_subset), "\n")
  
  # Additional checks for dimensions
  if (ncol(counts_subset) != nrow(sample_metadata)) {
    stop("Counts matrix and sample metadata dimensions do not match.")
  }
  
  # Creating a SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = counts_subset), colData = DataFrame(sample_metadata))
  
  return(se)
}


#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  # Creating a DESeqDataSet object
  dds <- DESeqDataSet(se, design = design)
  
  # Running DESeq2
  dds <- DESeq(dds)
  
  # Getting the DESeq2 results as a dataframe
  results <- results(dds)
  print(results)
  
  return(list(dds = dds, results = results))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`.
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  # Converting the DESeq2 results to a tibble
  labeled_results <- as_tibble(deseq2_res)
  
  # Adding a column to denote plotting status in volcano plot
  labeled_results$volc_plot_status <- case_when(
    (labeled_results$padj < padj_threshold) & (labeled_results$log2FoldChange > 0) ~ 'UP',
    (labeled_results$padj < padj_threshold) & (labeled_results$log2FoldChange < 0) ~ 'DOWN',
    TRUE ~ 'NS'
  )
  
  # Check the first few rows of labeled_results
  cat("Head of labeled_results:\n")
  print(head(labeled_results))
  
  # Check the data type of the "genes" column in labeled_results
  cat("Data type of 'genes' column in labeled_results:\n")
  print(class(labeled_results$genes))
  
  return(labeled_results)
}



#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  # Creating a ggplot object
  pval_plot <- ggplot(labeled_results, aes(x = padj)) +
    geom_histogram(bins = 30) +
    labs(x = 'Padj', y = 'Count') +
    ggtitle('Histogram of p-values')
  
  return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  # Creating a ggplot object
  log2fc_plot <- ggplot(labeled_results, aes(x = log2FoldChange)) +
    geom_histogram(bins = 30) +
    labs(x = 'Log2FC', y = 'Count') +
    ggtitle('Histogram of log2FC values')
  
  # Filtering to genes that are significant at the given padj threshold
  log2fc_plot <- log2fc_plot + geom_filter(aes(fill = volc_plot_status), data = labeled_results %>% filter(padj < padj_threshold))
  
  return(log2fc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes) {
  # Getting the top ten genes ranked by ascending padj
  top_genes <- labeled_results %>%
    filter(padj < 0.05) %>%
    arrange(padj) %>%
    head(num_genes)
  
  # Creating a DataFrame of normalized counts for the selected genes
  norm_counts_df <- counts(dds_obj, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column(var = 'gene_id') %>%
    filter(gene_id %in% top_genes$`gene_id`)
  
  # Converting the DataFrame to a tibble
  norm_counts_tibble <- as_tibble(norm_counts_df)
  
  # Gathering the data for the scatter plot
  norm_counts_tibble_long <- gather(norm_counts_tibble, key = 'Sample', value = 'Normalized_Counts', -gene_id)
  
  # Creating the scatter plot
  scatter_plot <- ggplot(norm_counts_tibble_long, aes(x = Sample, y = Normalized_Counts)) +
    geom_point() +
    labs(x = 'Sample', y = 'Normalized Counts') +
    ggtitle('Normalized Counts for Top Genes')
  
  return(scatter_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
volcano_plot <- function(labeled_results) {
  # Creating a ggplot object
  volcano_plot <- ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point() +
    labs(x = 'Log2FC', y = '-log10(padj)', color = 'Volcano plot status') +
    ggtitle('Volcano plot')
  
  # Adding horizontal and vertical lines to indicate significance thresholds
  volcano_plot <- volcano_plot + geom_hline(yintercept = -log10(0.1), linetype = 'dashed') +
    geom_vline(xintercept = 0, linetype = 'dashed')
  
  # Labeling the points by significance and fold change
  volcano_plot <- volcano_plot + geom_text_repel(data = labeled_results %>% filter(volc_plot_status == 'UP'), aes(label = Symbol), size = 3) +
    geom_text_repel(data = labeled_results %>% filter(volc_plot_status == 'DOWN'), aes(label = Symbol), size = 3)
  
  return(volcano_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

# Function to generate a named vector ranked by log2FC descending
# Debug code added to diagnose issues
make_ranked_log2fc <- function(results, id2gene_file) {
  # Read the id2gene data
  id2gene <- read.csv(id2gene_file, sep = '\t', header = FALSE, col.names = c("ensg", "genes"))
  
  # Merge the results and id2gene data frames by the 'genes' column
  labeled_results <- merge(results, id2gene, by = 'genes', all.x = TRUE)
  
  # Extract the log2FC values and set the row names to genes
  log2fc_values <- labeled_results$log2fc
  rownames(log2fc_values) <- labeled_results$genes
  
  return(log2fc_values)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # Running fgsea
  fgsea_results <- fgsea(gmt_file_path, rnk_list, minSize = min_size, maxSize = max_size)
  
  # Returning the fgsea results as a tibble
  return(as_tibble(fgsea_results))
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths) {
  # Filtering the fgsea results to the top ten positive and negative NES pathways
  fgsea_results_filtered <- fgsea_results %>%
    arrange(desc(NES)) %>%
    head(num_paths) %>%
    rbind(fgsea_results %>%
            arrange(NES) %>%
            head(num_paths))
  
  # Creating a ggplot object
  fgsea_plot <- ggplot(fgsea_results_filtered, aes(x = SET_NAME, y = NES)) +
    geom_bar(stat = 'identity') +
    labs(x = 'Pathway', y = 'NES') +
    ggtitle('Top ten positive and negative NES pathways')
  
  return(fgsea_plot)
}
