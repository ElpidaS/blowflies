# Differential Gene Expression in Crysomya

'''
We have the tsv files (6) of two types of individuals (female and male). 
3 reaptes each
embryos; 0-4h

We will use the 
'''
# load libraries
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)

##### 01 Input the transcripts.tsv files into R ####
# you can put the Directory that YOU have saved the .tsv files
setwd("C:\\Users\\Elpida\\Desktop\\scripts\\Heterozygocity_Rpart")

# Set the directory path
directory <- "C:\\Users\\Elpida\\Desktop\\scripts\\Heterozygocity_Rpart"

# Get a list of all .tsv files in the directory
tsv_files <- list.files(directory, pattern = "*.tsv", full.names = TRUE)


# Loop through each .tsv file and import the data as a separate table
for (file in tsv_files) {
  file_name <- basename(file)
  sample_name <- sub("^.*?_([^/]+)\\.transcript.tsv$", "\\1", file_name)
  
  x <- read.table(file, header = TRUE, sep = "\t")
  # delete all columns expect est_counts(4th) 
  #x <- x[, -c(2, 3, 5)]
  
  # Assign the name of the sample to the table
  assign(sample_name, x)
  
}


############### 02 for loop to rename the columns of each table #############

# Filter variables starting with "male" or "female" from all variable names 
# in the enviroment (get it with ls() )
filtered_vars <- grep("^male|^female", ls(), value = TRUE)
# get a list that contains the names of the files
file_name <- basename(filtered_vars)
sample_name <- sub("^.*?_([^/]+)\\.transcript.tsv$", "\\1", file_name)

# Create a list of the filtered variables
table_list <- mget(filtered_vars)
file_names <- c("f75", "f76", "f87","m98","m05","m06")
n <- 1

# this loop is possible because we allready know the sequence of variables in
# table list

for (i in table_list) {
  # make a c() with the names of the wanted column names
  col_names <- c("target_id", paste("sample", file_names[n], sep = "_"))
  # change the names of the columns 
  colnames(i) <- col_names
  # assign as the name of each file in table_list as the corresponding one
  # from file_names
  assign(sample_name[n],i)
  n <- n + 1
}

# output -> each sample file will have two columns 1st gene id, 2nd count (named after the sample)


################# 03 combine tables into one ####################

# List of tables to merge
table_list <- list(females_SRR5519675, females_SRR5519676, females_SRR5519687, males_SRR5519698,males_SRR5519705, males_SRR5519706)

# Merge tables based on the "target_id" column
combined_table <- Reduce(function(x, y) merge(x, y, by = "target_id", all = TRUE), table_list)

# asigns the values from the "gene" column of the 
# "all_counts" data frame as the row names of the data frame.
rownames(combined_table) <- combined_table$target_id

# delete the first columns as we have the info about genes as the rows' names
combined_table <- combined_table[, -c(1)]


# we have our main table;
# row names = gene names
# one column for each sample

#################### 04 Create Metadata ########################

# Create fake sample metadata
coldata <- data.frame(
  sex = rep(c("female", "male"), each =3),
  stringsAsFactors = FALSE
)
# add the names of the samples as row
row.names(coldata) <- colnames(combined_table)


#### matrix transformation + make all values numeric

my_matrix <- as.numeric(as.matrix(combined_table))

# turn values to intergers
combined_table <- round(combined_table)

#### 05 Prepapre for the PCA ######

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = combined_table, colData = coldata, design = ~ sex)

# Perform DESeq2 analysis and estimate size factors
dds <- DESeq(dds)
sizeFactors(dds) <- median(estimateSizeFactors(dds))

# Perform regularized log transformation (rlog)
rld <- rlog(dds)

# Extract PCA coordinates
pcaData <- plotPCA(rld, intgroup = "sex", returnData = TRUE)

# Create PCA plot
pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = sex)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), size = 3, hjust = 1, vjust = 2) + 
  theme_minimal()

# Print the PCA plot
print(pcaPlot)

#####################
###### SPM ##########
#####################

for (file in tsv_files) {
  file_name <- basename(file)
  sample_name <- sub("^.*?_([^/]+)\\.transcript.tsv$", "\\1", file_name)
  
  x <- read.table(file, header = TRUE, sep = "\t")
  # delete all columns expect est_counts(4th) 
  x <- x[, -c(2, 3, 4)]
  
  # Assign the name of the sample to the table
  assign(paste(sample_name,"TPM", sep = "_"), x)
  
}
###############
### 3 vs 3  ###
###############

all_TPM <- data.frame(females_SRR5519675_TPM$target_id,
             females_SRR5519675_TPM$tpm,
             females_SRR5519676_TPM$tpm,
             females_SRR5519687_TPM$tpm,
             males_SRR5519698_TPM$tpm,
             males_SRR5519705_TPM$tpm,
             males_SRR5519706_TPM$tpm)

rownames(all_TPM) <- females_SRR5519675_TPM$target_id
all_TPM <- all_TPM [, -c(1)]

mean_females <- rowMeans(all_TPM[, 1:3])
mean_males <- rowMeans(all_TPM[, 4:6])
TPM_means <- data.frame(mean_females,mean_males)

# filter TPMs ?
# use the rob's suggestion ?


# calculate SPM values  
TPM_SPMs <- TPM_means
TPM_SPMs$SPM <- (TPM_SPMs$mean_females^2) / ((TPM_SPMs$mean_males^2) + (TPM_SPMs$mean_females^2))

# sex bias
TPM_SPMs$sexbias <- 'unbiased'
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM > 0.7] <- 'Female-biased')
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM < 0.3] <- 'Male-biased')

# ggplot
ggplot(TPM_SPMs, aes(x = SPM, fill = sexbias)) +
  geom_histogram(bins = 30, colour = 'black', position = 'identity') +
  geom_vline(xintercept = c(0.3, 0.7), linetype = 'dotted', colour = 'black') +
  theme_classic() +
  scale_fill_manual("", breaks = c('Male-biased', 'Female-biased', 'unbiased'),
                    labels = c('male-biased', 'female-biased', 'unbiased'),
                    values = c('blue3', 'red3', 'grey40')) +
  labs(x = 'Female SPM', y = 'Frequency') +
  theme(plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

##############
#### ZGA ? ###
##############
ZGA_TPM <- data.frame(females_SRR5519675_TPM$target_id,
                      females_SRR5519687_TPM$tpm,
                      males_SRR5519706_TPM$tpm)

rownames(ZGA_TPM) <- females_SRR5519675_TPM$target_id
ZGA_TPM <- ZGA_TPM [, -c(1)]

# calculate SPM values  
TPM_SPMs <- ZGA_TPM
TPM_SPMs$SPM <- (TPM_SPMs$females_SRR5519687_TPM.tpm^2) / ((TPM_SPMs$males_SRR5519706_TPM.tpm^2) + (TPM_SPMs$females_SRR5519687_TPM.tpm^2))

# sex bias
TPM_SPMs$sexbias <- 'unbiased'
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM > 0.7] <- 'Female-biased')
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM < 0.3] <- 'Male-biased')

# ggplot
ggplot(TPM_SPMs, aes(x = SPM, fill = sexbias)) +
  geom_histogram(bins = 30, colour = 'black', position = 'identity') +
  geom_vline(xintercept = c(0.3, 0.7), linetype = 'dotted', colour = 'black') +
  theme_classic() +
  scale_fill_manual("", breaks = c('Male-biased', 'Female-biased', 'unbiased'),
                    labels = c('male-biased', 'female-biased', 'unbiased'),
                    values = c('blue3', 'red3', 'grey40')) +
  labs(x = 'Female SPM', y = 'Frequency') +
  theme(plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))

###############
### 2 vs 2  ###
# MAternal ? ##
###############

mat_TPM <- data.frame(females_SRR5519675_TPM$target_id,
                      females_SRR5519675_TPM$tpm,
                      females_SRR5519676_TPM$tpm,
                      males_SRR5519698_TPM$tpm,
                      males_SRR5519705_TPM$tpm)

rownames(mat_TPM) <- females_SRR5519675_TPM$target_id
mat_TPM <- mat_TPM [, -c(1)]

mean_females <- rowMeans(mat_TPM[, 1:2])
mean_males <- rowMeans(mat_TPM[, 3:4])
TPM_means <- data.frame(mean_females,mean_males)

# filter TPMs ?
# use the rob's suggestion ?


# calculate SPM values  
TPM_SPMs <- TPM_means
TPM_SPMs$SPM <- (TPM_SPMs$mean_females^2) / ((TPM_SPMs$mean_males^2) + (TPM_SPMs$mean_females^2))

# sex bias
TPM_SPMs$sexbias <- 'unbiased'
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM > 0.7] <- 'Female-biased')
TPM_SPMs <- within(TPM_SPMs, sexbias[SPM < 0.3] <- 'Male-biased')

# ggplot
ggplot(TPM_SPMs, aes(x = SPM, fill = sexbias)) +
  geom_histogram(bins = 30, colour = 'black', position = 'identity') +
  geom_vline(xintercept = c(0.3, 0.7), linetype = 'dotted', colour = 'black') +
  theme_classic() +
  scale_fill_manual("", breaks = c('Male-biased', 'Female-biased', 'unbiased'),
                    labels = c('male-biased', 'female-biased', 'unbiased'),
                    values = c('blue3', 'red3', 'grey40')) +
  labs(x = 'Female SPM', y = 'Frequency') +
  theme(plot.title = element_text(size = 14, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
