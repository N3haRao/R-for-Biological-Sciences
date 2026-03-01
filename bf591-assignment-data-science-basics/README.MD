# Data Science Basics: Analyzing High-Dimensional Gene Expression Data


## Problem Statement
High dimensional data is typically time and resource-intensive to analyze, interpret and visualize. Gene expression data (including microarray data) often consists of measurements for tens of thousands of known genes for every sample. Learning methods to inspect, reduce and display high dimensional data is necessary for many machine learning and bioinformatics problems.

---

## Learning Objectives
- Understand the fundamental concepts of **Principal Component Analysis (PCA)** and **Hierarchical Clustering** / **Heatmaps**.
- Learn how to **perform PCA in R** and evaluate its outputs.
- Generate basic **clustered heatmaps** to visualize gene expression data.

---

## Skill Set
- **R Functions**: `scale()`, `transpose(t())`, `prcomp()`, `heatmap()`, `ggplot2`
- **Bioinformatics**: Understanding of microarray technology and its application in gene expression analysis.

---

## Background on Microarrays
A **microarray** consists of thousands of specifically designed single-stranded DNA sequences bound to a solid surface (such as glass or nylon). In a typical experiment:
- Extracted RNA or DNA from the sample is labeled with fluorescent dyes.
- These labeled samples are hybridized with the bound probes by flowing across the surface.
- Non-hybridized molecules are washed away.
- A laser excites the attached dye, emitting light, which is detected by a scanner and converted into a digital image.
- Image processing then quantifies gene expression levels from the intensity of the emitted light.

---

## Background on Principal Component Analysis (PCA)
**Principal Component Analysis (PCA)** is an exploratory data analysis technique used to reduce the dimensionality of data while preserving as much variance as possible. PCA creates new, uncorrelated linear combinations of the original variables (called principal components), ordered by the amount of variance they explain. The first principal component (PC1) explains the most variance, the second (PC2) explains the second most, and so on.

PCA is useful in bioinformatics for visualizing and understanding patterns in high-dimensional data, where a few principal components can often capture most of the variance in the dataset, providing a more compact and interpretable representation.

---


## Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391
The example intensity data was taken from the listed publication. In it, the authors proposed the use of gene expression profiling (through microarray technology) to generate a robust and reproducible classifier to identify subtypes of colorectal cancer samples. They identified six subtypes that they demonstrated to be significantly associated with distinct molecular pathways and clinical pathologies. We have provided you a sample expression matrix that was normalized and subjected to batch correction. We will use the example_intensity_data.csv as our data matrix going forward.


---

## Scaling Data with `scale()` in R
When analyzing gene expression data, variables (i.e., gene expression levels) often vary by several orders of magnitude. For example, some genes may be expressed at very low levels (single digits), while others may have counts in the tens of thousands. This discrepancy can make it difficult to visualize or compare data effectively.

To address this, we can standardize or scale the data. The `scale()` function in R centers and scales each column of a data matrix by subtracting the mean and dividing by the standard deviation. This results in Z-scores, where a Z-score of 0 indicates the mean, and a Z-score of 1 corresponds to a value one standard deviation above the mean. This process normalizes the data while retaining the inherent patterns of expression.

To scale the data by gene (rather than by sample), the matrix needs to be transposed using `transpose(t())`, so that scaling is applied within each gene rather than across samples.

### Task: Scaling the Data
1. Read in the `example_intensity_data.csv` as a data frame.
2. Apply `scale()` to center and scale the matrix.

---

## Proportion of Variance Explained by Each Principal Component
Each principal component (PC) in PCA represents a direction of maximal variance in the dataset. It is often helpful to understand how much of the total variance is captured by each principal component. 

### Task: Variance Explained
1. Write a function to compute the proportion of variance explained by each PC.
2. Generate a bar chart showing the variance explained by each principal component.
3. Plot a cumulative sum of variance explained as a scatter and line plot.

**Deliverables**:
- A vector containing the variance explained by each PC.
- A tibble with the PC labels, variance explained, and cumulative variance.
- A bar chart with the variance explained, overlayed with a scatter plot of cumulative variance.

---

## PCA Plotting and Visualization
PCA is often used to visualize the similarity between samples in high-dimensional data. A scatterplot of the first two principal components (PC1 vs. PC2) can reveal clusters of samples, which may correspond to experimental variables such as genotype or disease subtype.

### Task: PCA Plot
1. Plot PC1 vs. PC2 based on the PCA results (`pca_results$x`).
2. Annotate the samples based on their classification (e.g., **SixSubtypesClassification**) and color the points by subtype (C3 vs. C4).

**Deliverables**:
- A scatter plot of PC1 vs. PC2 with sample points colored by their classification.

---

## Hierarchical Clustering and Heatmaps
A **heatmap** is a powerful tool for visualizing the expression of genes across multiple samples. It represents gene expression values as colors, allowing for the identification of patterns and clusters. Before generating a heatmap, hierarchical clustering is often applied to group genes or samples based on similarity.

### Task: Heatmap of Differentially Expressed Genes
We will use hierarchical clustering and heatmaps to visualize the expression of **differentially expressed genes**. 

1. Filter the `differential_expression_results.csv` to find probes with an adjusted p-value < 0.01.
2. Extract the corresponding normalized intensity values from `example_intensity_data.csv`.
3. Create a clustered heatmap using `heatmap()` in R, visualizing the expression of these significant genes.

**Considerations**:
- Use appropriate color palettes for accessibility (e.g., **RColorBrewer**'s colorblind-friendly palettes).

**Deliverables**:
- A list of significant probes (adjusted p-value < 0.01).
- A matrix containing only the significant probes.
- A clustered heatmap of the normalized intensity values for the differentially expressed probes.

---

## References
- [DNA Microarray Technology - NIH](https://www.genome.gov/about-genomics/fact-sheets/DNA-Microarray-Technology)
- Govindarajan, R., et al. "Microarray and its applications." *Journal of Pharmacy & Bioallied Sciences* 4, S310 (2012).
- Trevino, V., et al. "DNA microarrays: A powerful genomic tool for biomedical and clinical research." *Molecular Medicine* 13, 527–541 (2007).
- Shlens, J. "A Tutorial on Principal Component Analysis." [arXiv](https://arxiv.org/abs/1404.1100)
- Lever, J., et al. "Points of Significance: Principal Component Analysis." *Nature Methods* 14, 641–642 (2017).
- Hastie, T., et al. "The elements of statistical learning: Data mining, inference, and prediction." Springer, 2001. [Book PDF](https://hastie.su.domains/ElemStatLearn/printings/ESLII_print12_toc.pdf)
