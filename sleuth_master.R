# Author: Markus Wolfien
# Script for RNA-Seq quantified data analysis with sleuth
# See also https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html

# Preprocessing
# Demultiplexed data with flexbar
# Preprocessed with trimmomatic and pseudoaligned with kallisto (version: 0.43.1)
# Be aware: 
# Please generate bootstraps using "kallisto quant -b"
# The bootstrap is only for estimating the technical variance, it will not improve the original point estimate. The sleuth program, http://pachterlab.github.io/sleuth/ , uses bootstraps to improve differential expression detection. 
# This is why you will get slightly different est_counts in the bootstraps, and the amount of variation gives you an indication of how reliable the initial point estimate is.

# clear workspace
rm(list = ls())

# set working directory
setwd("$PATH-to-Directory")

# Shortcut do save current graphs within the display device as pdf 
dopdf <- function(filename = "dopdf.pdf", pdf.width = 7, pdf.height = 5) {
  dev.copy2pdf(file=filename, width = pdf.width, height = pdf.height)
}

# For use ...  
dopdf()

# To install sleuth start R and first install rhdf5 by typing:
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")

# Then install devtools by typing
install.packages("devtools")

# and install sleuth by typing
devtools::install_github("pachterlab/sleuth")

# afterwards you should be able to open the library
library("sleuth")

#############################################################
######### Read in the data and do the pre formatting ########
#############################################################


# While working in Rstudio - you can disable annoying warning messages about using gui's with this command 
suppressMessages({library("sleuth")}) 


# The first step in a sleuth analysis is to specify where the kallisto results are stored. A variable is created for this purpose with

sample_id <- dir(file.path("kallisto_linux-v0.43.1/hg19_mrna", "results"))

# The result can be displayed by typing

sample_id

#A list of paths to the kallisto results indexed by the sample IDs is collated with

kal_dirs <- file.path("kallisto_linux-v0.43.1/hg19_mrna", "results", sample_id, "")
kal_dirs

# The next step is to load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples:
  
s2c <- read.table(file.path("/kallisto_linux-v0.43.1/hg19_mrna", "meta", "seq_info.txt"), header = TRUE, stringsAsFactors=FALSE)
# s2c <- dplyr::select(s2c, sample, condition)
s2c

# Now the directories must be appended in a new column to the table describing the experiment. 
# This column must be labeled path, otherwise sleuth will report an error. 
# This is to ensure that samples can be associated with kallisto quantifications.

s2c <- dplyr::mutate(s2c, path = kal_dirs)

# It is important to check that the pairings are correct:
  
print(s2c)


#############################################################
######### Perform the main DE detection with sleuth #########
#############################################################

# If you want to include the gene names directly, please go to the next step!


# Next, the “sleuth object” can be constructed. 
# This object will store not only the information about the experiment, but also details of the model to be used for differential testing, and the results. 
# It is prepared and used with four commands that 
# (1) load the kallisto processed data into the object 
# (2) estimate parameters for the sleuth response error measurement (full) model 
# (3) estimate parameters for the sleuth reduced model, and 
# (4) perform differential analysis (testing) using the likelihood ratio test. On a laptop the four steps should take about a few minutes altogether.

# The sleuth object must first be initialized with
so_0 <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# Then the full model is fit with
so_0 <- sleuth_fit(so_0, ~condition, 'full')

# What this has accomplished is to “smooth” the raw kallisto abundance estimates for each sample using a linear model with a parameter that represents the experimental condition (in this case non_res vs. res). 
# To test for transcripts that are differential expressed between the conditions, sleuth performs a second fit to a “reduced” model that presumes abundances are equal in the two conditions. 
# To identify differential expressed transcripts sleuth will then identify transcripts with a significantly better fit with the “full” model.

# The “reduced” model is fit with
so_0 <- sleuth_fit(so_0, ~1, 'reduced')

# and the test is performed with
so_0 <- sleuth_lrt(so_0, 'reduced', 'full')

# In general, sleuth can utilize the likelihood ratio test with any pair of models that are nested, and other walkthroughs illustrate the power of such a framework for accounting for batch effects and more complex experimental designs.
# The models that have been fit can always be examined with the models() function.

models(so_0)

# The results of the test can be examined with

sleuth_table_0 <- sleuth_results(so_0, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_0 <- dplyr::filter(sleuth_table_0, pval <= 0.05) # qval/pval
head(sleuth_significant_0, 20)

###############################################################################################
######### Including gene names into transcript-level analysis       ##########
##############################################################################
# At this point the sleuth object constructed from the kallisto runs has information about the data, the experimental design, the kallisto estimates, the model fit, and the testing. 
# In other words it contains the entire analysis of the data. There is, however, one piece of information that can be useful to add in, but that is optional. 
# In reading the kallisto output sleuth has no information about the genes transcripts are associated with, but this can be added allowing for searching and analysis of significantly differential transcripts by their associated gene names.
# Since the example was constructed with the ENSEMBL human transcriptome, we will add gene names from ENSEMBL using biomaRt:

# installing Biomart
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")
# Open library
library("biomaRt")

# Then collect gene names with
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

write.csv(t2g, file="t2g.csv")


# and add them into the sleuth table with
so <- sleuth_prep(s2c, target_mapping = t2g)

# Repeat the analyses with target gene names
so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) # qval/pval
head(sleuth_significant, 20)

# You can also save the created table into a csv
write.csv(sleuth_significant, file="sleuth_significant.csv")

# This addition of metadata to transcript IDs is very general, and can be used to add in other information.
# The easiest way to view and interact with the results is to generate the sleuth live site that allows for exploratory data analysis:
sleuth_live(so)

# Among the tables and visualizations that can be explored with sleuth live are a number of plots that provide an overview of the experiment. 
# A PCA plot provides a visualization of the samples:
  
plot_pca(so, color_by = 'condition')
plot_pca(so, color_by = 'condition', text_labels = TRUE)


# Various quality control metrics can also be examined. 
# The count distributions for each sample (grouped by condition) can be displayed using the plot_group_density command:
  
plot_group_density(so, use_filtered = TRUE, units = "est_counts", trans = "log", grouping = setdiff(colnames(so$sample_to_covariates), "sample"), offset = 1)


# Convert a sleuth object to a matrix with the condition names.
# A list with an attribute 'data', which contains a matrix of target_ids and transcript expression in which_units
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
head(sleuth_matrix$data) # look at first 5 transcripts, sorted by name

write.csv(sleuth_matrix, file="sleuth_matrix.csv")

# write.xlsx(sleuth_matrix, file="sleuth_matrix.xlsx")


#######################################
######### Search and plot Genes #########
#######################################


# Search for Genes within the gene data set!
geneA <- t2g[grep("^GeneA", t2g$ext_gene),]


t2gmatch <- c("^GeneA","^GeneB","^GeneC","^GeneD","^GeneF","^GeneG","^GeneH","Gene_n")
to_match <- t2g[grep(paste(t2gmatch,collapse="|"), t2g$ext_gene),]
to_match[1:75,]
candidates <- t(sleuth_matrix$data[grep(paste(to_match$target_id[1:75],collapse="|"), rownames(sleuth_matrix$data)),])
candidates
write.csv(candidates, file="candidates.csv")

# Search for ENST ID's within the data set!
sleuth_matrix$data[grep("^ENST00000XXX", rownames(sleuth_matrix$data)),]
rownames(sleuth_matrix$data)[grep("^ENST00000XXX", rownames(sleuth_matrix$data))]

plot_bootstrap(so_0, "ENST00000XXX.X", units = "est_counts", color_by = "condition")


models(so)
tests(so)
so <- sleuth_wt(so, 'condition')

sleuth_genes <- sleuth_gene_table(so, 'condition', test_type ='wt', which_group = 'ext_gene')
head(sleuth_genes) # show info for first 5 genes
sleuth_genes[1:5, 6] # show transcripts for first 5 genes

# Mean-variance relationship
#
# Plot the mean-variance relationship of transcripts as modeled in \code{sleuth}. 
# Each dot represents a transcript. The blue dots represent the transcripts that are used in the shrinkage estimation. The fitted curve represents the smooth fit.
#
# The x-axis represents the mean expression of transcripts pooled across all samples. 
# The y-axis represents the 'biological' variance after the technical variance has been removed.
plot_mean_var(so)
