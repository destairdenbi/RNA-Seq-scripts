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
setwd("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/Plots")

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

sample_id <- dir(file.path("/users/sbi00/mw581/Documents/kallisto_linux-v0.43.1/hg19_mrna", "st_results"))

# The result can be displayed by typing

sample_id

#A list of paths to the kallisto results indexed by the sample IDs is collated with

kal_dirs <- file.path("/users/sbi00/mw581/Documents/kallisto_linux-v0.43.1/hg19_mrna", "st_results", sample_id, "")
kal_dirs

# The next step is to load an auxillary table that describes the experimental design and the relationship between the kallisto directories and the samples:
  
s2c <- read.table(file.path("/users/sbi00/mw581/Documents/kallisto_linux-v0.43.1/hg19_mrna", "meta", "seq_info_2.txt"), header = TRUE, stringsAsFactors=FALSE)
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

# Plotting specific single genes as bar charts per condition
plot_bootstrap(so_0, "ENST00000341259.6", units = "est_counts", color_by = "VEGF") # SH2B3-201
plot_bootstrap(so_0, "ENST00000538307.1", units = "est_counts", color_by = "condition") # SH2B3-202
plot_bootstrap(so_0, "ENST00000550925.2", units = "est_counts", color_by = "condition") # SH2B3-203
plot_bootstrap(so_0, "ENST00000579793.6", units = "est_counts", color_by = "condition") # Notch2-NL

# Plot the most sign. gene as quality control 
plot_bootstrap(so_0, "ENST00000012049.9", units = "est_counts", color_by = "condition") # QPTCL



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
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl",host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

write.csv(t2g, file="/media/sdc/mw581/Projects/Steinhoff/t2g.csv")


# and add them into the sleuth table with
so <- sleuth_prep(s2c, target_mapping = t2g)

# Repeat the analyses with target gene names
so <- sleuth_fit(so, ~condition, 'full')

so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, pval <= 0.05) # qval/pval
head(sleuth_significant, 20)

# You can also save the created table into a csv
write.csv(sleuth_significant, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant.csv")

# This addition of metadata to transcript IDs is very general, and can be used to add in other information.
# The easiest way to view and interact with the results is to generate the sleuth live site that allows for exploratory data analysis:
sleuth_live(so)

# Among the tables and visualizations that can be explored with sleuth live are a number of plots that provide an overview of the experiment. 
# A PCA plot provides a visualization of the samples:
  
plot_pca(so, color_by = 'condition')
plot_pca(so, color_by = 'condition', text_labels = TRUE)
plot_pca(so, color_by = 'lanes', text_labels = TRUE)
plot_pca(so, color_by = 'SH2B3', text_labels = TRUE)
plot_pca(so, color_by = 'VEGF', text_labels = TRUE)
plot_pca(so, color_by = 'EPO', text_labels = TRUE)
plot_pca(so, color_by = 'Vitronectin', text_labels = TRUE)
plot_pca(so, color_by = 'Throm', text_labels = TRUE)
plot_pca(so, color_by = 'CD34_EPC', text_labels = TRUE)


# Various quality control metrics can also be examined. 
# The count distributions for each sample (grouped by condition) can be displayed using the plot_group_density command:
  
plot_group_density(so, use_filtered = TRUE, units = "est_counts", trans = "log", grouping = setdiff(colnames(so$sample_to_covariates), "sample"), offset = 1)


# Convert a sleuth object to a matrix with the condition names.
# A list with an attribute 'data', which contains a matrix of target_ids and transcript expression in which_units
sleuth_matrix <- sleuth_to_matrix(so, 'obs_norm', 'tpm')
head(sleuth_matrix$data) # look at first 5 transcripts, sorted by name

write.csv(sleuth_matrix, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_matrix.csv")

# write.xlsx(sleuth_matrix, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_matrix.xlsx")


#######################################
######### Search and plot Genes #########
#######################################


# Search for Genes within the gene data set!
t2g[grep("^NFKB", t2g$ext_gene),]   # NFKb hits
sh2b3 <- t2g[grep("^SH2B3", t2g$ext_gene),]
t2g[grep("^TP53", t2g$ext_gene),]   # p53 hits
Vegfb <- t2g[grep("^VEGFB", t2g$ext_gene),]   # VEGF hits
t2g[grep("^MT", t2g$ext_gene),]   # MT hits
cd34 <- t2g[grep("^CD34", t2g$ext_gene),]   # cd34 hits
Prom1 <- t2g[grep("^PROM1", t2g$ext_gene),]   # cd133 hits
Notch2 <- t2g[grep("^NOTCH2", t2g$ext_gene),]   # notch hits
kit <- t2g[grep("^KIT", t2g$ext_gene),]   # c-kit hits
gata4 <- t2g[grep("^GATA4", t2g$ext_gene),]   # gata4 hits
atxn1 <- t2g[grep("ATXN1", t2g$ext_gene),]   # sca1 hits

t2gmatch <- c("^SH2B3","^VEGFB","^CD34","^PROM1","^NOTCH2","^KIT","^GATA4","ATXN1")
to_match <- t2g[grep(paste(t2gmatch,collapse="|"), t2g$ext_gene),]
to_match[1:75,]
candidates <- t(sleuth_matrix$data[grep(paste(to_match$target_id[1:75],collapse="|"), rownames(sleuth_matrix$data)),])
candidates
write.csv(candidates, file="/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/correlation/candidates.csv")

# Search for ENST ID's within the data set!
sleuth_matrix$data[grep("^ENST00000341259", rownames(sleuth_matrix$data)),]
rownames(sleuth_matrix$data)[grep("^ENST00000341259", rownames(sleuth_matrix$data))]

plot_bootstrap(so_0, "ENST00000229239.9", units = "est_counts", color_by = "VEGF") # SH2B3-201
plot_bootstrap(so_0, "ENST00000309422.6", units = "est_counts", color_by = "CD34_EPC") # VEGF
plot_bootstrap(so_0, "ENST00000376838.5", units = "est_counts", color_by = "SH2B3") # mTOR



models(so)
tests(so)
so <- sleuth_wt(so, 'conditionres')

sleuth_genes <- sleuth_gene_table(so, 'conditionres', test_type ='wt', which_group = 'ext_gene')
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


##############################################
######### Further DE detection ################
##############################################

install.packages("cowplot")
library("cowplot")


# Removal of outliers can greatly improve results, and is sometimes warranted due to botched sample prep, problems with sequencing or occasionally book-keeping errors/sample mix-ups. 
# However removing outliers can also accidentally, if not intentionally, become a form of data fishing.
# Therefore before removing outliers it is prudent to try to understand why a sample might be an outlier. 
# To do so it is helpful to examine the PCA loadings, i.e. the primary genes whose linear combinations define the principal components:

plot_loadings(so, pc_input = 1)
plot_loadings(so, pc_input = 2)
plot_loadings(so, pc_input = 3)



# The first gene driving PC1 is

plot_bootstrap(so_0, 'ENST00000251595.10', color_by = 'condition')
plot_bootstrap(so_0, 'ENST00000397806.1', color_by = 'condition')


# The outlier can be removed with

sample_to_condition <- dplyr::filter(sample_to_condition, sample != 'SRR1654638')

# Once the sample has been removed it is useful to re-examine the data at a high level. 
# The first two principal components shown in the PCA plot now reveal a better behaved experiment and the separation between conditions is evident:
  
so <- sleuth_prep(sample_to_condition, target_mapping = ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)


# this is your DF1    
DF1 <- data.frame("A" = rep(c("A","B"), 18),
                  "B" = rep(c("C","D","E"), 12),
                  "NUM"= rep(rnorm(36,10,1)),
                  "TEST" = rep(NA,36))

#this is a DF2 i created, with unique A, B, VAL
DF2 <- data.frame("A" = rep(c("A","B"),3),
                  "B" = rep(c("C","D","E"),2),
                  "VAL" = rep(1:6))

# and this is the answer of what i assume you want      
tmp <- merge(DF1,DF2, by=c("A","B"), all.x=TRUE, all.y=FALSE)
DF1[4] <- tmp[5]

############################################
######### Subgroup analysis       ##########
############################################
print(s2c)

# Subgroup 1
s2c_1 <- s2c[c(3,5,9,11),]
so_1 <- sleuth_prep(s2c_1, target_mapping = t2g)
so_1 <- sleuth_fit(so_1, ~condition, 'full')
so_1 <- sleuth_fit(so_1, ~1, 'reduced')
so_1 <- sleuth_lrt(so_1, 'reduced', 'full')
sleuth_table_1 <- sleuth_results(so_1, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_1 <- dplyr::filter(sleuth_table_1, pval <= 0.05) # qval/pval
head(sleuth_significant_1, 20)
so_1 <- sleuth_wt(so_1, 'conditionres')
sleuth_genes_1 <- sleuth_gene_table(so_1, 'conditionres', test_type ='wt', which_group = 'ext_gene')
head(sleuth_genes_1) # show info for first 5 genes
sleuth_genes_1[1:5, 6] # show transcripts for first 5 genes
write.csv(sleuth_significant_1, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant_1.csv")
sleuth_live(so_1)


# Subgroup 2 -- appears not to work ... more than 2 samples are needed -- this two groups are added to group 3 within the analysis of (joint) group 5
s2c_2 <- s2c[c(1,4),]
print(s2c_2)
so_2 <- sleuth_prep(s2c_2, target_mapping = t2g)
so_2 <- sleuth_fit(so_2, ~condition, 'full')
so_2 <- sleuth_fit(so_2, ~1, 'reduced')
so_2 <- sleuth_lrt(so_2, 'reduced', 'full')
sleuth_table_2 <- sleuth_results(so_2, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_2 <- dplyr::filter(sleuth_table_2, pval <= 0.05) # qval/pval
head(sleuth_significant_2, 20)
write.csv(sleuth_significant_2, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant_2.csv")
sleuth_live(so_2)

# Subgroup 3
s2c_3 <- s2c[c(6,8,10,14),]
print(s2c_3)
so_3 <- sleuth_prep(s2c_3, target_mapping = t2g)
so_3 <- sleuth_fit(so_3, ~condition, 'full')
so_3 <- sleuth_fit(so_3, ~1, 'reduced')
so_3 <- sleuth_lrt(so_3, 'reduced', 'full')
sleuth_table_3 <- sleuth_results(so_3, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_3 <- dplyr::filter(sleuth_table_3, pval <= 0.05) # qval/pval
head(sleuth_significant_3, 20)
so_3 <- sleuth_wt(so_3, 'conditionres')
sleuth_genes_3 <- sleuth_gene_table(so_3, 'conditionres', test_type ='wt', which_group = 'ext_gene')
head(sleuth_genes_3) # show info for first 5 genes
sleuth_genes_3[1:5, 6] # show transcripts for first 5 genes
write.csv(sleuth_significant_3, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant_3.csv")
sleuth_live(so_3)

# Subgroup 4
s2c_4 <- s2c[c(2,7,12,13),]
print(s2c_4)
so_4 <- sleuth_prep(s2c_4, target_mapping = t2g)
so_4 <- sleuth_fit(so_4, ~condition, 'full')
so_4 <- sleuth_fit(so_4, ~1, 'reduced')
so_4 <- sleuth_lrt(so_4, 'reduced', 'full')
sleuth_table_4 <- sleuth_results(so_4, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_4 <- dplyr::filter(sleuth_table_4, pval <= 0.05) # qval/pval
head(sleuth_significant_4, 20)
so_4 <- sleuth_wt(so_4, 'conditionres')
sleuth_genes_4 <- sleuth_gene_table(so_4, 'conditionres', test_type ='wt', which_group = 'ext_gene')
head(sleuth_genes_4) # show info for first 5 genes
sleuth_genes_4[1:5, 6] # show transcripts for first 5 genes
write.csv(sleuth_significant_4, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant_4.csv")
sleuth_live(so_4)

# Subgroup 5 --- joint group of 2 and 3
s2c_5 <- s2c[c(1,4,6,8,10,14),]
print(s2c_5)
so_5 <- sleuth_prep(s2c_5, target_mapping = t2g)
so_5 <- sleuth_fit(so_5, ~condition, 'full')
so_5 <- sleuth_fit(so_5, ~1, 'reduced')
so_5 <- sleuth_lrt(so_5, 'reduced', 'full')
sleuth_table_5 <- sleuth_results(so_5, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant_5 <- dplyr::filter(sleuth_table_5, pval <= 0.05) # qval/pval
head(sleuth_significant_5, 20)
so_5 <- sleuth_wt(so_5, 'conditionres')
sleuth_genes_5 <- sleuth_gene_table(so_5, 'conditionres', test_type ='wt', which_group = 'ext_gene')
head(sleuth_genes_5) # show info for first 5 genes
sleuth_genes_5[1:5, 6] # show transcripts for first 5 genes
write.csv(sleuth_significant_5, file="/media/sdc/mw581/Projects/Steinhoff/sleuth_significant_5.csv")
sleuth_live(so_5)

##############################################
######### Machine learning stuff ##############
##############################################

## load libraries
library(caret) # for Naive Bayes
library(klaR) # libraries needed by caret
library(adabag) # for adaBoost
library(tree)
library(readxl) # read Excel files
library(fastAdaboost)
library(randomForest)
library(e1071) # package for SVMs (svmLinear2)
library(kernlab) # package for SVMs
library(plyr)
library(mlbench)

library(tsne) # unsupervised ML algorithm
library(scatterplot3d) # for visualizing super cool 3d scatterplots
library(rgl) # for interactive 3d plots
library(rsm) # for curve 3d fitting
library(reshape2)
library(ggplot2)

library(plot3D)# library for 4d graphs
library(plot3Drgl)#library for rotating 4d graphs

# select which features to omit from machine learning 
omit_features <- 
  c(
    "RTC...",
    NA
  )

#remove unwanted features from data
data_select <- sleuth_matrix
data_select2 <- t(data_select$data) # transposing the matrix for the needed input!


# use this if you want to apply the omiited ones to your data >> data_select <- data[ , -which(names(data) %in% omit_features)] <<

# replace NANs with zeros. If there might be any ... 
data_select2[is.na(data_select2)] <- 0

# remove features with nearly zero variance
nzv_cols <- nearZeroVar(data_select2)
if(length(nzv_cols) > 0) data_select2 <- data_select2[,-nearZeroVar(data_select2)]

# Reordering Matrix to add meta data
s2c_n <- s2c[c("8","1", "2", "14", "10", "6", "12", "4", "7", "9", "11", "13", "5", "3"),]

data <- cbind(s2c_n, data_select2)

# [OPTIONAL] next step is to remove highly correlated features -----
# skipping this step for now
# calculate correlation matrix
correlationMatrix <- cor(data_select2[,1:10000])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
# print indexes of highly correlated attributes
print(highlyCorrelated)

tmp <- cor(data_select2[,1:100000])
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

# data.new <- data_select2[,1:10000]0[, apply(tmp,2,function(x) all(x<=0.99))]
library(directlabels)
colors = rainbow(length(unique(data$Res)))
names(colors) = unique(data$condition)
col.list=c("red","grey")

ecb = function(x,y){ plot(x, color = colors[data$Res], pch=16, type="p")}
tsne_2 = tsne(data[,48000:50000], epoch_callback = ecb, perplexity=30, k=2 , max_iter = 130, epoch = 5, initial_dims = 10, whiten = TRUE, min_cost =0.1)
par(new=TRUE)
plot(tsne_2, xaxt='n', ann=FALSE, yaxt='n', frame.plot=FALSE, type = "p", pch=19, col=col.list[data$Res])
par(new=TRUE)
plot(tsne_2,t='n', xaxt='n', ann=FALSE, yaxt='n', frame.plot=FALSE); text(tsne_2, col=colors[data$Res])
# Epoch: Iteration #150 error is: 0.495513394050989 with k=2 perplex=30, 10 dims, 48000-50000

ecb_2_3d = function(x,y,z){ scatterplot3d(x, color = col.list[data$Res], pch=16, type="p")}
tsne2_3d = tsne(data[,13:141209], epoch_callback = ecb_2_3d, perplexity=30, k=3 , max_iter = 350, epoch = 50, initial_dims = 10, whiten = FALSE, min_cost =0.1)

plotrgl() # transforms the last plot of plot3D into fancy rotating version



## Data fitting for 3d tsne  --------
DATA <- as.data.frame(cbind(1:14, tsne2_3d[,1],tsne2_3d[,2],tsne2_3d[,3]))
x_wt <- DATA$V1
y_disp <- DATA$V2
z_mpg <- DATA$V3
fit <- lm(z_mpg ~ poly(x_wt, y_disp, degree = 3), data = DATA)
image(fit, y_disp ~ x_wt)
contour(fit, y_disp ~ x_wt)
persp(fit, y_disp ~ x_wt, zlab = "z_mpg")
# Use rsm package to create surface model.
SurfMod <- contour(fit, y_disp ~ x_wt)
# extract list values from rsm Surface Model 
Xvals <- SurfMod$`x_wt ~ y_disp`[1]
Yvals <- SurfMod$`x_wt ~ y_disp`[2]
Zvals <- SurfMod$`x_wt ~ y_disp`[3]
# Construct matrix with col and row names 
SurfMatrix <- Zvals$z
colnames(SurfMatrix) <- Yvals$y
rownames(SurfMatrix) <- Xvals$x

# Convert matrix to data frame
SurfDF <- melt(SurfMatrix)

gg <- ggplot(data = SurfDF) + 
  geom_tile(data = SurfDF, aes(Var1, Var2,z = value, fill = value)) +
  stat_contour(data = SurfDF, aes(Var1, Var2, z = value, color = ..level..)) +
  scale_colour_gradient(low = "green", high = "red") 


direct.label(gg, "bottom.pieces")
par(new=TRUE)
# adds the single patient
plot(tsne2_3d,t='n', xaxt='n', ann=FALSE, yaxt='n', frame.plot=FALSE); text(tsne2_3d, col=colors[data$Res])
# adds the class

plot(tsne2_3d, xaxt='n', ann=FALSE, yaxt='n', frame.plot=FALSE, type = "p", pch=19, col=col.list[data$Res])
legend("topright", "(x,y)", c("Responder","Non-responder"), pch = 19, col=c("red","grey"), box.lty = 0, cex = 0.85, title = "True class", inset = 0)


########################################################
######## ROC - curves ---------------
########################################################

library(caret)
library(mlbench)
data(Sonar)
ctrl <- trainControl(method="cv", 
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)
rfFit <- train(Class ~ ., data=Sonar, 
               method="rf", preProc=c("center", "scale"), 
               trControl=ctrl)
library(pROC)
# Select a parameter setting
selectedIndices <- rfFit$pred$mtry == 2
# Plot:
plot.roc(rfFit$pred$obs[selectedIndices],
         rfFit$pred$M[selectedIndices])

#########################################################
# 4D tSNE #################
# Add small dots on basal plane and on the depth plane
scatter3D_fancy <- function(x, y, z,..., colvar = z)
{
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 3, add = TRUE, colkey = FALSE)
    
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 3, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
            colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75)) 
}

ecb_2_4d_f = function(x,y,z,w){ scatter3D_fancy(x[,1], x[,2], x[,3],colvar = x[,4],bty = "u", ticktype = "detailed",col.grid = "white",col.panel ="grey95", pch = 20,cex = 2, phi = 0)}
tsne2_4d_f = tsne(data[,30:1000], epoch_callback = ecb_2_4d_f, perplexity=10, k=4 , max_iter = 2000, epoch = 50, initial_dims = 450, whiten = TRUE, min_cost =0)

col.list=c("red","grey")
scatter3D(x = tsne2_4d_f[,1], y = tsne2_4d_f[,2], z = tsne2_4d_f[,3], add = TRUE, pt.cex = 2.5, cex = 2.5, colkey = FALSE, col = col.list[data$lanes])
# adds the single patient
text3D(x = tsne2_4d_f[,1], y = tsne2_4d_f[,2], z = tsne2_4d_f[,3], labels = data$sample, add = TRUE, colkey = FALSE, cex = 0.6)
legend("right", "(x,y)", c("Responder","Non-responder"), pch = 21, col=c("red","grey"), box.lty = 0, cex = 0.9, title = "True class", inset = 0, pt.lwd = 1.5, pt.cex = 1.5)

########################################################
######### Co-expression analyses ##################
###################################################

library(WGCNA)
library(cluster)
library(ProCoNA)
library(gplots)
library(RColorBrewer)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
ALLOW_WGCNA_THREADS=24

# Take a quick look at what is in the data set:
dim(sleuth_matrix);
head(sleuth_matrix);

datExpr02 = as.data.frame(t(sleuth_matrix$data));

names(datExpr02) = rownames(sleuth_matrix$data);

gsg2 = goodSamplesGenes(datExpr02);
gsg2$allOK

# if not 

if (!gsg2$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg2$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr02)[!gsg2$goodGenes], collapse = ", ")));
  if (sum(!gsg2$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr02)[!gsg2$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr02 = datExpr02[gsg2$goodSamples, gsg2$goodGenes]
}



Tree = hclust(dist(datExpr02), method = "ward.D2");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(Tree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

abline(h = 60000, col = "red");

clust = cutreeStatic(Tree, cutHeight = 2000000, minSize =5)
table(clust)
keepSamples = (clust==1)
datExpr2 = datExpr02[keepSamples,]
nGenes = ncol(datExpr2)
nSamples = nrow(datExpr2)

collectGarbage();
Tree2 = hclust(dist(datExpr02), method = "complete");
plot(Tree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)



## Automatic construction of the gene network and identifcation of modules
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 4, to=40, by=4))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr02, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.75,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

datExpr2<-sapply(datExpr02, as.numeric)
datExpr2[datExpr2 == 0] <- 0.000001


ADJ1=abs(cor(datExpr2[,1:20000],use="p"))^10

# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1

hierADJ=hclust(as.dist(dissADJ), method="average" )
colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=50))
# Plot the resulting clustering tree together with the true color assignment
sizeGrWindow(10,5);
plotDendroAndColors(hierADJ, colors = colorStaticADJ, dendroLabels = FALSE, hang = 0.03,
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )


branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )


colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ,
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))
# Plot results of all module detection methods together:
sizeGrWindow(10,5)
plotDendroAndColors(dendro = hierADJ,
                    colors=data.frame(colorStaticADJ,
                                      colorDynamicADJ, colorDynamicHybridADJ),
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")


dissTOM = 1-TOMsimilarityFromExpr(datExpr2[,1:20000], power = 10);
dissTOM2 = 1-TOMsimilarityFromExpr(datExpr2[,1:60000], power = 10);

# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average");
# The reader should vary the height cut-off parameter h1
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.8, minSize=30))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.99,
                                                    deepSplit=2, pamRespectsDendro = FALSE))
# Now we plot the results
sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=data.frame(colorStaticTOM,
                                      colorDynamicTOM, colorDynamicHybridTOM),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")




k=softConnectivity(datE=datExpr2,power=7)
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
#########################################
#datExpr=datExpr2[, rank(-k,ties.method="first" )<=14000]
#gene.names= names(datExpr)
#plotNetworkHeatmap(datExpr, plotGenes = gene.names,
#                   networkType="signed", useTOM=FALSE,
#                   power=7, main="signed correlations")

#ADJ2=abs(cor(datExpr,use="p"))^10
#dissADJ=1-ADJ2
#hierADJ=hclust(as.dist(dissADJ), method="average" )
#sizeGrWindow(10,5);
#plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03,
#                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )


#dissTOM=TOMdist(ADJ1)
#dissTOM2=TOMdist(ADJ2)

#hierADJ=hclust(as.dist(dissADJ), method="average" )
#plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03,
#                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )

#cmd1=cmdscale(as.dist(dissTOM2),2)
#collectGarbage()

####################################### ?

net = blockwiseModules(datExpr2, power = 10,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "testTOM",
                       verbose = 3)

table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[30]], mergedColors[net$blockGenes[[30]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net$dendrograms[[30]], mergedColors[net$blockGenes[[30]]],
                    "Module colors",abHeight = 0.99,
                    dendroLabels = FALSE )

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

net2 = blockwiseModules(datExpr2, power = 10,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       getTOMs = NULL,
                       saveTOMs = FALSE,
                       maxBlockSize = 10000,
                       saveTOMFileBase = "blockwiseTOM",
                       deepSplit = 2,
                       detectCutHeight = 0.995,
                       verbose = 2)


table(net2$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors2 = labels2colors(net2$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net2$dendrograms[[2]], mergedColors2[net2$blockGenes[[2]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(net2$dendrograms[[2]], mergedColors2[net2$blockGenes[[2]]],
                    "Module colors",abHeight = 0.99,
                    dendroLabels = FALSE )

moduleLabels2 = net2$colors
moduleColors2 = labels2colors(net2$colors)
MEs2 = net2$MEs;
geneTree2 = net2$dendrograms;

##################
##### code chunks from -- https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-blockwise.R
############

bwLabels = matchLabels(net2$colors, moduleLabels2, pThreshold = 1e-7);
bwColors = labels2colors(bwLabels)

# Here a more flexible way of plotting several trees and colors on one page
sizeGrWindow(12,6)
pdf(file = "BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6);
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4);
nBlocks = length(net2$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks)
  plotDendroAndColors(net2$dendrograms[[block]], moduleColors2[net2$blockGenes[[block]]],
                      "Module colors", 
                      main = paste("Gene dendrogram and module colors in block", block), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      setLayout = FALSE)
dev.off();


## Quantifying module–trait associations - 1. Entering traits 
library(readxl) # read Excel files
traitData <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/cardiac_data_machine_learning/Perfect_table_trait_final_clin_res170510_1t.xlsx", sheet = "Tabelle1") # all pre operative time point
traitData <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/cardiac_data_machine_learning/Perfect_table_trait3_final_clin_res170510_1t.xlsx", sheet = "Tabelle1") # selected pre operative time point
traitData <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/cardiac_data_machine_learning/Perfect_table_trait20_final_clin_res170510_1t.xlsx", sheet = "Tabelle1") # selected pre operative time point from Machine learning analysis of Steinhoff et al. 2017a
dim(traitData)
names(traitData)

Samples = rownames(datExpr02);
traitRows = match(Samples, traitData$RTC_Nr);
datTraits = traitData[traitRows, -1];
datTraits[is.na(datTraits)] <- 0
# rownames(datTraits) = traitData[traitRows, 1];

nzv_cols <- nearZeroVar(datTraits)
if(length(nzv_cols) > 0) datTraits <- datTraits[, -nearZeroVar(datTraits)]

nums <- sapply(datTraits, is.numeric)

datTraits2 <- datTraits[ , nums]

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr2), method = "complete")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits2, signed = FALSE);
sizeGrWindow(6,10)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, autoColorHeight = TRUE, colorHeight = 0.4, cex.rowText = 0.2,cex.colorLabels = 0.7, addTextGuide = TRUE,
                    groupLabels = names(datTraits2),
                    main = "Sample dendrogram and trait heatmap")

plotDendroAndColors(Tree2, traitColors, autoColorHeight = TRUE, colorHeight = 0.4, cex.rowText = 0.1,cex.colorLabels = 0.7, addTextGuide = TRUE,
                    groupLabels = names(datTraits2),
                    main = "Sample dendrogram and trait heatmap") # Tree from above measurement

dopdf(filename = "Dendro_and_traits_20.pdf")

########################################################
# Quantifying module–trait associations - 2. Integrating expression and traits
# Define numbers of genes and samples
nGenes = ncol(datExpr2);
nSamples = nrow(datExpr2);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr2, moduleColors2)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(20,13.5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits2),
               cex.lab.x = 0.7,
               cex.lab.y = 0.4,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.15,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


dopdf(filename = "module_trait_relationship_20.pdf", pdf.width = 20, pdf.height = 13.5)


#######################################


##################################################
################# Tried gene sub grouping ########
##################################################

library(sna) # this is needed for closeness
#library(MASS)
#library(class)
library(cluster)
#library(sma) # different from sna! this is needed for plot.mat below
library(impute) # needed for imputing missing value before principal component analysis
#library(splines) # for the spline predictor to estimate the number of clusters
library(Hmisc)    # probably you won’t need this
library(bio3d)

thresholds1= c(seq(.1,.5, by=.1), seq(.55,.9, by=.05) )

RdichotTable=pickHardThreshold(datExpr2, cutvector= thresholds1)[[2]]

tau= 0.74

corhelp=cor(datExpr02[,1:20000],use="pairwise.complete.obs")


AdjMat1 = I(abs(corhelp)>tau); diag(AdjMat1)=0

Degree = apply(AdjMat1,2,sum)
# Let’s create a scale free topology plot.
# The black curve corresponds to scale free topology and
# the red curve corresponds to truncated scale free topology.
par(mfrow=c(1,1))
scaleFreePlot(Degree, truncated1=TRUE); 

# ===================================================
# Module Detection/TOMk Plots
# ===================================================

# Restrict the analysis to the genes with high connectivity
DegCut = 1000
DegreeRank = rank(-Degree)
rest1 = DegreeRank <= DegCut
sum(rest1)
# Comment: this is different from DegCut due to ties in the ranks...

# Remove any isolated nodes, they affect the correlation between different GTOM measures
conn.comp1 = component.dist(AdjMat1[rest1,rest1])   # need the SNA library
conn.comp1$csize
rest1[rest1] = (conn.comp1$membership==1);

summary(Degree[rest1])

datSummary = datExpr02[,1:20000]

AdjMat1rest = AdjMat1[rest1,rest1]
datSummaryrest=datSummary[rest1,]
rm(DegreeRank,AdjMat1); collect_garbage();

# Compute the distance measures
dissCor1 = 1 - abs(corhelp[rest1,rest1])
dissCor6 = 1 - corhelp[rest1,rest1]^6
dissGTOM0 = 1 - AdjMat1rest; diag(dissGTOM0) = 0
# Compute the generalized TOM based dissimilarity 
dissGTOM1 <- TOMdist(AdjMat1rest,TOMType = "unsigned", verbose = 1)
dissGTOM2 <- TOMdist(AdjMat1rest,TOMType = "unsigned", verbose = 2)
dissGTOM3 <- TOMdist(AdjMat1rest,TOMType = "unsigned", verbose = 3)

# This creates a hierarchical clustering using the TOM matrix 
hierGTOM0 = hclust(as.dist(dissGTOM0),method="average")
hierGTOM1 = hclust(as.dist(dissGTOM1),method="average")
hierGTOM2 = hclust(as.dist(dissGTOM2),method="average")
hierGTOM3 = hclust(as.dist(dissGTOM3),method="average")

par(mfrow=c(1,4))
plot(hierGTOM0, main="Adjacency Matrix: GTOM0", labels=F)
plot(hierGTOM1, main="Standard TOM Measure: GTOM1", labels=F)
plot(hierGTOM2, main="New TOM Measure: GTOM2", labels=F)
plot(hierGTOM3, main="New TOM Measure: GTOM3", labels=F)

# This suggests a height cut-off of h1 for the GTOM1 but one should check the robustness of this choice

colorGTOM1=as.character(labels2colors(hierGTOM1,h1=.95,minsize1=40))
colorGTOM2=as.character(modulecolor2(hierGTOM2,h1=.4))

par(mfrow=c(3,3),mar=c(2,2))
plot(hierGTOM0, main="Adjacency Matrix: GTOM0", labels=F)
plot(hierGTOM1, main="Standard TOM Measure: GTOM1", labels=F)
plot(hierGTOM2, main="New TOM Measure: GTOM2", labels=F)
hclustplot(hierGTOM0,colorGTOM1, title1="")
hclustplot(hierGTOM1,colorGTOM1, title1="Colored by GTOM1 modules")
hclustplot(hierGTOM2,colorGTOM1, title1="")
hclustplot(hierGTOM0,colorGTOM2, title1="")
hclustplot(hierGTOM1,colorGTOM2, title1="Colored by GTOM2 modules")
hclustplot(hierGTOM2,colorGTOM2, title1="")


##############################################################
########## Correlation analyses --------------
##############################################################

# Correlations with significance levels
library(Hmisc)
library(readxl) # read Excel files
library(gplots)

corr <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/correlation.xlsx", sheet = "seq_0") # full dataset
corr <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/correlation.xlsx", sheet = "seq_10") # full dataset
corr <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/correlation.xlsx", sheet = "red") # full dataset
corr <- read_excel("/users/sbi00/mw581/CS/home/Desktop/MMRP/Perfect/Sequencing/correlation.xlsx", sheet = "red_full") # full dataset



# calculate correlation matrix
nums <- sapply(corr, is.numeric)
corr_num <- corr[ , nums]
# corr[,2:108]<-sapply(corr[,2:108], as.numeric)
# replace NANs with zeros. If there might be any ... 
corr_num[is.na(corr_num)] <- 0
correlationMatrix <- cor(corr_num)
# summarize the correlation matrix
print(correlationMatrix)
correlationMatrix[is.na(correlationMatrix)] <- 0
# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
# print indexes of highly correlated attributes
print(highlyCorrelated)

tmp <- cor(corr_num)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0

heatmap(correlationMatrix, margins = c(10, 10))
dopdf(filename = "cor_heatmap.pdf", pdf.width = 16, pdf.height = 11)

corr_high <- corr_num[ , highlyCorrelated]
correlationMatrix_high <- cor(corr_high)
heatmap(correlationMatrix_high, margins = c(10, 10))
dopdf(filename = "high_86_cor_heatmap.pdf", pdf.width = 16, pdf.height = 11)


highlyCorrelated2 <- findCorrelation(correlationMatrix, cutoff=0.80)
# print indexes of highly correlated attributes
print(highlyCorrelated2)
corr_high2 <- corr_num[ , highlyCorrelated2]
correlationMatrix_high2 <- cor(corr_high2)
heatmap.2(correlationMatrix_high2, margins = c(10, 10))
heatmap.2(correlationMatrix_high2, margins = c(13, 13),density.info=c("none"),trace=c("none"), cexRow = 0.25 ,cexCol = 0.25)
dopdf(filename = "high_80_cor_heatmap_legend_10_seq_mean_fill.pdf", pdf.width = 13, pdf.height = 13)


rcorr <- rcorr(as.matrix(corr_num), type="spearman") # type can be pearson or spearman
highlyCorrelated3 <- findCorrelation(rcorr$r, cutoff=0.75)
# print indexes of highly correlated attributes
print(highlyCorrelated3)
corr_high3 <- corr_num[ , highlyCorrelated3]
correlationMatrix_high3 <- cor(corr_high3)
heatmap.2(correlationMatrix, margins = c(12, 12),density.info=c("none"),trace=c("none"), cexRow = 1 ,cexCol = 1)
dopdf(filename = "cor_heatmap_legend_red_sh2b3_combined.pdf", pdf.width = 13, pdf.height = 13)

# First Correlogram Example
library(corrgram)
corrgram(corr_high, order=TRUE, lower.panel=panel.shade,
         upper.panel=panel.pie,
         main="het")

corrgram(corr_high, order=TRUE, lower.panel=panel.shade,
         upper.panel=NULL, font.labels = 0.2, cex.labels = 0.3,
         main="Perfect Correlation Analysis")

corrgram(corr_high2, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pie,
        font.labels = 0.2, cex.labels = 0.4,
         main="Perfect Correlation Analysis (cor>0.86)")

dopdf(filename = "high_86_cor_heatmap_alternative.pdf", pdf.width = 16, pdf.height = 11)

corrgram(rcorr$r, order=TRUE, lower.panel=panel.shade,
         upper.panel=NULL, font.labels = 0.2, cex.labels = 0.4,
         main="Perfect Correlation Analysis")


# data imputation for correlation?
library(Hmisc)
library(missForest)

# impute with mean value example
NA2mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
corr_num[] <- lapply(corr_num, NA2mean)

