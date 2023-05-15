# Original date (Rob) 28.05.20
# Edit (Elpida) 15.5.23
# blobtools bestsum output: filtering out contaminant contigs

setwd('C:/Users/Elpida/Desktop') 
library(data.table) 

# read in blobtools table, name columns appropriately
# in our case is either the AF_* or TF_* (should change manually)
contigs.bestsum <- read.table('TF_blobplot.blobDB.table.txt', stringsAsFactors = FALSE)
colnames(contigs.bestsum) <- c('name', 'length', 'GC', 'N' , 'bam0', 'phylum.t.6%s', 'phylum.s.7%s', 'phylum.c.8')
head(contigs.bestsum)
sum(contigs.bestsum$length)

# filter out athropoda contigs (for keeping)
contigs.athropoda <- contigs.bestsum[(contigs.bestsum$`phylum.t.6%s` == 'Arthropoda'),]
# tottall length of all the contings and their mean length
sum(contigs.athropoda$length)
mean(contigs.athropoda$length)

# filter out no hit contigs (also for keeping)
# we keep no hits because we wave reasons beliving that they also belong to athropoda
contigs.no.hit <- contigs.bestsum[(contigs.bestsum$`phylum.t.6%s` == 'no-hit'),]
sum(contigs.no.hit$length)
mean(contigs.no.hit$length)

# combine Arthropoda + no hits
contigs.athropoda.no.hit <- rbind(contigs.athropoda,contigs.no.hit)
nrow(contigs.athropoda.no.hit)
sum(contigs.athropoda.no.hit$length)

# only include contigs greater than 220mb in length and have more than 2 coverage (bam0 > 2)
# According to Andere et al. 2020 expected genome size is 426 Mb
# 220 and higher allow us to get a genome (contigs) of 430 Mb in total
scop.contigs.keep <- contigs.athropoda.no.hit[(contigs.athropoda.no.hit$length > 220) & (contigs.athropoda.no.hit$bam0 > 2),]

head(scop.contigs.keep)
nrow(scop.contigs.keep)
sum(scop.contigs.keep$length)

contig_keep_list <- scop.contigs.keep[, 1]
contig_keep_col <- data.frame(matrix(unlist(contig_keep_list), nrow=length(contig_keep_list), byrow=T))
head(contig_keep_col)
nrow(contig_keep_col)

write.table(contig_keep_col, file = 'contig_keep_list_TF.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)


