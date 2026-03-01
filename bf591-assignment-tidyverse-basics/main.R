library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data_subset.csv')
read_expression_table <- function(filename) {
  # Read the CSV file
  expression_data <- readr::read_csv(filename)
  
  # Extract subject IDs from the first column (assuming it contains subject IDs)
  subject_ids <- expression_data[, 1]
  
  # Add a "subject_id" column
  expression_data <- expression_data %>%
    dplyr::mutate(subject_id = subject_ids)
  
  return(expression_data)
}




#' Load Metadata from Specified CSV File
#'
#' This function reads the provided CSV file into a dataframe.
#'
#' @param filepath (character) The path to the CSV file.(data/proj_metadata.csv)
#'
#' @return A dataframe containing the loaded metadata.
#'

load_metadata <- function(filepath) {
  #Using function to read the csv file
  metadata <- readr::read_csv(filepath)
  #Returning loaded metadata
  return(metadata)
  }



#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {
  return(stringr::str_replace_all(str, "\\.", "_"))
}


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  data <- data %>%
    dplyr::rename(
      Age = Age_at_diagnosis,
      Subtype = SixSubtypesClassification,
      Batch = normalizationcombatbatch
    ) %>%
    dplyr::select(Sex, Age, TNM_Stage, Tumor_Location, geo_accession, KRAS_Mutation, Subtype, Batch)
  
  return(data)
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
  data <- data %>%
    dplyr::mutate(Stage = paste("stage", TNM_Stage)) %>%
    dplyr::mutate(Stage = factor(Stage))
  
  return(data)
}



#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  mean_age <- data %>%
    dplyr::filter(Sex == sex) %>%
    dplyr::pull(Age) %>%
    mean(na.rm = TRUE)
  
  return(mean_age)
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage. Name the newly created column containing the average, 'mean_avg'
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  age_summary <- data %>%
    dplyr::group_by(Stage) %>%
    dplyr::summarize(mean_avg = mean(Age, na.rm = TRUE))
  
  return(age_summary)
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
  cross_tab <- data %>%
    dplyr::count(Stage, Subtype) %>%
    tidyr::pivot_wider(names_from = Subtype, values_from = n, values_fill = 0)
  
  return(cross_tab)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
summarize_expression <- function(exprs) {
  # Filtering out non-numeric values and convert them to numeric values
  numeric_exprs <- exprs %>%
    select_if(is.numeric) %>%
    mutate_all(as.numeric)
  
  # Calculating mean and variance for (the now numeric0 columns
  avg_expression <- colMeans(numeric_exprs, na.rm = TRUE)
  probe_variability <- apply(numeric_exprs, 2, var, na.rm = TRUE)
  
  # Results: tibble
  summary_data <- tibble(
    mean_exp = avg_expression,
    variance = probe_variability,
    probe = colnames(numeric_exprs)
  )
  
  return(summary_data)
}
