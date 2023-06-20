#### SPMs  ####
### filtered, to acount for possible different TPM distribution in each sample

# load libraries
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(data.table)

setwd("C:\\Users\\Elpida\\Desktop\\scripts\\Heterozygocity_Rpart")

# Set the directory path
directory <- "C:\\Users\\Elpida\\Desktop\\scripts\\Heterozygocity_Rpart"

# Get a list of all .tsv files in the directory
tsv_files <- list.files(directory, pattern = "*.tsv", full.names = TRUE)

# 1% filter
for (file in tsv_files) {
  file_name <- basename(file)
  sample_name <- sub("^.*?_([^/]+)\\.transcript.tsv$", "\\1", file_name)
  
  x <- read.table(file, header = TRUE, sep = "\t")
  # delete all columns expect est_counts(4th) 
  x <- x[, -c(2, 3, 4)]
  
  # Calculate the bottom 1% threshold
  bottom_threshold <- quantile(x$tpm, probs = 0.01)

  # Subset the data based on the condition TPM > bottom_threshold
  x<- x[x$tpm > bottom_threshold, ]

  # Assign the name of the sample to the table
  assign(paste(sample_name,"TPM_01", sep = "_"), x)
  
}


###############
### 3 vs 3  ###
###############
# after the filtering each sample has a different number of genes
# when we merge we should only include the ones that exist in all 
# of them (all = FALSE)

# Merge data frames based on row names
data_frames <- list(females_SRR5519675_TPM_01, 
                    females_SRR5519676_TPM_01, 
                    females_SRR5519687_TPM_01, 
                    males_SRR5519698_TPM_01, 
                    males_SRR5519705_TPM_01, 
                    males_SRR5519706_TPM_01)

all_TPM <- Reduce(function(x, y) merge(x, y, by = "target_id", all = FALSE), data_frames)

colnames(all_TPM)[2] <- "f75"
colnames(all_TPM)[3] <- "f76"
colnames(all_TPM)[4] <- "f87"
colnames(all_TPM)[5] <- "m98"
colnames(all_TPM)[6] <- "m05"
colnames(all_TPM)[7] <- "m06"

rownames(all_TPM) <- all_TPM$target_id
all_TPM <- all_TPM [, -c(1)]

mean_females <- rowMeans(all_TPM[, 1:3])
mean_males <- rowMeans(all_TPM[, 4:6])
TPM_means <- data.frame(mean_females,mean_males)


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

###################
#### 2 seperate ###
###################
# Merge data frames based on row names
data_frames <- list(females_SRR5519687_TPM_01, 
                    males_SRR5519706_TPM_01)

all_TPM_2 <- Reduce(function(x, y) merge(x, y, by = "target_id", all = FALSE), data_frames)

colnames(all_TPM_2)[2] <- "f87"
colnames(all_TPM_2)[3] <- "m06"

rownames(all_TPM_2) <- all_TPM_2$target_id
all_TPM_2 <- all_TPM_2 [, -c(1)]


# calculate SPM values  
TPM_SPMs <- all_TPM_2
TPM_SPMs$SPM <- (TPM_SPMs$f87^2) / ((TPM_SPMs$f87^2) + (TPM_SPMs$m06^2))

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
###############
# Merge data frames based on row names
data_frames <- list(females_SRR5519675_TPM_01, 
                    females_SRR5519676_TPM_01, 
                    males_SRR5519698_TPM_01, 
                    males_SRR5519705_TPM_01)

TPM_2vs2 <- Reduce(function(x, y) merge(x, y, by = "target_id", all = FALSE), data_frames)

colnames(TPM_2vs2)[2] <- "f75"
colnames(TPM_2vs2)[3] <- "f76"
colnames(TPM_2vs2)[4] <- "m98"
colnames(TPM_2vs2)[5] <- "m05"


rownames(TPM_2vs2) <- TPM_2vs2$target_id
TPM_2vs2 <- TPM_2vs2 [, -c(1)]

mean_females <- rowMeans(TPM_2vs2[, 1:2])
mean_males <- rowMeans(TPM_2vs2[, 3:4])
TPM_means <- data.frame(mean_females,mean_males)


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

########### MODELS of SEX DETERMINATION ################
## 01 Functional
# there should be one(/or more) female-limited transcripts

## 02 Non Functional
# there should be one(/or more) male biased transcripts 

### model 01 ####
# SPM = 1 --- zero transcripts
Model_01_transcripts_f_lim <- TPM_SPMs[TPM_SPMs$SPM == 1, , drop = FALSE] 
# SPM > 0.7 --- 1691
Model_01_transcripts_f_biased <- TPM_SPMs[TPM_SPMs$SPM > 0.7, , drop = FALSE]
# SPM > 0.9 --- 254
Model_01_transcripts_f_biased09 <-TPM_SPMs[TPM_SPMs$SPM > 0.9, , drop = FALSE]

# SPM > 0.98 --- 67
# use that one because 67 is a manageable number and >0.98 is pretty close to 1
Model_01_transcripts_f_biased098 <-TPM_SPMs[TPM_SPMs$SPM > 0.98, , drop = FALSE]
# Write the row names to the text file
writeLines(rownames(Model_01_transcripts_f_biased098), con = "model01_98.txt")

# SPM > 0.99 --- 0
Model_01_transcripts_f_biased099 <-TPM_SPMs[TPM_SPMs$SPM >= 0.99, , drop = FALSE]


#### model 02 ####
# male double (1.75-2.25) more than female (expression)
# in SPM -> 0.16 < SPM < 0.24 --- 817
Model_02_transcripts <- TPM_SPMs[TPM_SPMs$SPM > 0.16 & TPM_SPMs$SPM < 0.24, , drop = FALSE] 
# for double (x2) is SPM == 0.2 -> 0.2 <SPM< 0.21 ---90
Model_02_transcripts_2x <- TPM_SPMs[TPM_SPMs$SPM > 0.2 & TPM_SPMs$SPM < 0.21, , drop = FALSE] 
writeLines(rownames(Model_02_transcripts_2x), con = "model02_x2.txt")
