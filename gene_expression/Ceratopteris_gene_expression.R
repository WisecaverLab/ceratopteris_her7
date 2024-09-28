library(stringr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
rm(list = ls())

# Transcript to genes lookup table
t2gfile <- 'Desktop/gene_expression_jen_local/10_gene_exp/transcripts2genes.txt'
tx2gene <- read.table(t2gfile, sep="\t", header=TRUE)

######### male v herm overview ######################################################################
# Path to kallisto abundance files
filesfile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_abundance_files_male_v_herm.txt'
files <- scan(filesfile, what = 'list')
file_names <- str_split(files, '/', simplify = TRUE)
file_names <- file_names[,ncol(file_names)-1]
names(files) <- file_names

# sample conditions
samplefile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_sample_conditions_male_v_herm.txt'
sampleTable <- read.table(samplefile, sep="\t", header=TRUE)

# import kallisto counts with tximport
btxi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")

# create deseq dataset
bddsTxi <- DESeqDataSetFromTximport(btxi, colData = sampleTable, design = ~ condition)
bdds <- DESeq(bddsTxi)
#res <- results(bdds, contrast=c("condition","Male","Herm"), alpha = 0.1, lfcThreshold = 0.58)
res <- results(bdds, contrast=c("condition","Male","Herm"), alpha = 0.1, lfcThreshold = 1)
summary(res)

#res
#shrink log fold changes
bresLFC <- lfcShrink(bdds, contrast=c("condition","Male","Herm"), res = res, type="ashr")
summary(bresLFC)
write.table(as.data.frame(bresLFC), file="Desktop/gene_expression_jen_local/10_gene_exp/burow_male_v_herm_de_results.txt", quote = FALSE, sep = "\t")

# plot DE results
#ggmaplot(bresLFC, fdr = 0.1, fc = 1.5, size = 0.1, top = 0, legend = "top", label.select = c("Ceric.05G026200", 'Ceric.29G061800', 'Ceric.1Z203700', 'Ceric.1Z024300', 'Ceric.07G096400'), label.rectangle = TRUE)
ggmaplot(bresLFC, fdr = 0.1, fc = 2, size = 0.1, top = 0, legend = "top", label.select = c("Ceric.05G026200", 'Ceric.29G061800', 'Ceric.1Z203700', 'Ceric.1Z024300', 'Ceric.07G096400'), label.rectangle = TRUE)
#parta <- ggmaplot(bresLFC, main = "male v WT hermaphrodite", fdr = 0.1, fc = 1.5, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
parta <- ggmaplot(bresLFC, main = "male v WT hermaphrodite", fdr = 0.1, fc = 2, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
parta <- ggpar(parta, ylim = c(-12, 12))
#parta <- ggdraw(parta) + draw_label("4868", x = 0.9, y = 0.8, color = "#B31B21", size = 7, fontface = "bold") + draw_label("4885", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
parta <- ggdraw(parta) + draw_label("3464", x = 0.9, y = 0.8, color = "#B31B21", size = 7, fontface = "bold") + draw_label("3087", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
parta

######### her19 v herm overview ######################################################################
# Path to kallisto abundance files
filesfile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_abundance_files_her19_v_herm.txt'
files <- scan(filesfile, what = 'list')
file_names <- str_split(files, '/', simplify = TRUE)
file_names <- file_names[,ncol(file_names)-1]
names(files) <- file_names

# sample conditions
samplefile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_sample_conditions_her19_v_herm.txt'
sampleTable <- read.table(samplefile, sep="\t", header=TRUE)

# import kallisto counts with tximport
btxi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")

# create deseq dataset
bddsTxi <- DESeqDataSetFromTximport(btxi, colData = sampleTable, design = ~ condition)
bdds <- DESeq(bddsTxi)
#res <- results(bdds, contrast=c("condition","her19","Herm"), alpha = 0.1, lfcThreshold = 0.58)
res <- results(bdds, contrast=c("condition","her19","Herm"), alpha = 0.1, lfcThreshold = 1)
summary(res)

#shrink log fold changes
bresLFC <- lfcShrink(bdds, contrast=c("condition","her19","Herm"), res = res, type="ashr")
summary(bresLFC)
write.table(as.data.frame(bresLFC), file="Desktop/gene_expression_jen_local/10_gene_exp/burow_her19_v_herm_de_results.txt", quote = FALSE, sep = "\t")

# plot DE results
ggmaplot(bresLFC, fdr = 0.1, fc = 1.5, size = 0.1, top = 0, legend = "top")
#partb <- ggmaplot(bresLFC, main = "her19 v WT hermaphrodite", fdr = 0.1, fc = 1.5, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
partb <- ggmaplot(bresLFC, main = "her19 v WT hermaphrodite", fdr = 0.1, fc = 2, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
partb <- ggpar(partb, ylim = c(-12, 12))
#partb <- ggdraw(partb) + draw_label("209", x = 0.9, y = 0.8, color = "#B31B21", size = 7, fontface = "bold") + draw_label("1982", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
partb <- ggdraw(partb) + draw_label("44", x = 0.9, y = 0.8, color = "#B31B21", size = 7, fontface = "bold") + draw_label("1422", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
partb

######### male v her19 overview ######################################################################
# Path to kallisto abundance files
filesfile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_abundance_files_male_v_her19.txt'
files <- scan(filesfile, what = 'list')
file_names <- str_split(files, '/', simplify = TRUE)
file_names <- file_names[,ncol(file_names)-1]
names(files) <- file_names

# sample conditions
samplefile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_sample_conditions_male_v_her19.txt'
sampleTable <- read.table(samplefile, sep="\t", header=TRUE)

# import kallisto counts with tximport
btxi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")

# create deseq dataset
bddsTxi <- DESeqDataSetFromTximport(btxi, colData = sampleTable, design = ~ condition)
bdds <- DESeq(bddsTxi)
#res <- results(bdds, contrast=c("condition","Male","her19"), alpha = 0.1, lfcThreshold = 0.58)
res <- results(bdds, contrast=c("condition","Male","her19"), alpha = 0.1, lfcThreshold = 1)
summary(res)

#shrink log fold changes
bresLFC <- lfcShrink(bdds, contrast=c("condition","Male","her19"), res = res, type="ashr")
summary(bresLFC)
write.table(as.data.frame(bresLFC), file="Desktop/gene_expression_jen_local/10_gene_exp/burow_male_v_her19_de_results.txt", quote = FALSE, sep = "\t")

# plot DE results
#ggmaplot(bresLFC, fdr = 0.1, fc = 1.5, size = 0.1, top = 0, legend = "top", label.select = c("Ceric.05G026200", 'Ceric.29G061800', 'Ceric.1Z203700', 'Ceric.1Z024300', 'Ceric.07G096400'), label.rectangle = TRUE)
#partc <- ggmaplot(bresLFC, main = "male v her19", fdr = 0.1, fc = 1.5, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
partc <- ggmaplot(bresLFC, main = "male v her19", fdr = 0.1, fc = 2, size = 0.1, alpha = 0.5, top = 0, legend = "none", ggtheme = ggplot2::theme_classic(base_size = 6, base_line_size = 0.25))
partc <- ggpar(partc, ylim = c(-12, 12))
#partc <- ggdraw(partc) + draw_label("6222", x = 0.9, y = 0.9, color = "#B31B21", size = 7, fontface = "bold") + draw_label("5568", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
partc <- ggdraw(partc) + draw_label("4727", x = 0.9, y = 0.9, color = "#B31B21", size = 7, fontface = "bold") + draw_label("3754", x = 0.9, y = 0.3, color = "#1465AC", size = 7, fontface = "bold")
partc

########## Combined dataset for boxplots ######################################################################
filesfile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_abundance_files.txt'
files <- scan(filesfile, what = 'list')
file_names <- str_split(files, '/', simplify = TRUE)
file_names <- file_names[,ncol(file_names)-1]
names(files) <- file_names

samplefile <- 'Desktop/gene_expression_jen_local/10_gene_exp/burow_sample_conditions.txt'
sampleTable <- read.table(samplefile, sep="\t", header=TRUE)

ctxi <- tximport(files, type = "kallisto", tx2gene = tx2gene, countsFromAbundance = "no")
cddsTxi <- DESeqDataSetFromTximport(ctxi, colData = sampleTable, design = ~ condition)

rlogd <- rlog(cddsTxi)
normalized_rlog_matrix <- assay(rlogd)
write.table(as.data.frame(normalized_rlog_matrix), file="Desktop/gene_expression_jen_local/10_gene_exp/burow_combined_kallisto_normalizedcounts_rlog_matrix.txt", quote = FALSE, sep = "\t" )




cddsTxi$condition <- relevel(cddsTxi$condition, ref = "Herm")
cdds <- DESeq(cddsTxi)

# PCA plot of sample variance
rld <- rlog(cdds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic() +
  xlim(-100, 100) +
  ylim(-100, 100)
ggsave2('Desktop/gene_expression_jen_local/PCA_supplemental_figure.pdf', device = "pdf", width = 10, height = 8, units = "cm")

my_comparisons <- list( c("her19", "Male"), c("Herm", "Male"))

d <- plotCounts(cdds, gene='Ceric.29G061800', intgroup="condition", returnData=TRUE)
#ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + theme_cowplot(font_size = 7) +labs(y = "count", x = NULL, title = 'HER7', subtitle = 'Ceric.29G061800')
her7 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + theme_cowplot(font_size = 7) +labs(y = "count", x = NULL, title = 'HER7', subtitle = 'Ceric.29G061800')
#her7 <- her7 + stat_compare_means(comparisons = my_comparisons, size = 2) 
her7

d <- plotCounts(cdds, gene='Ceric.1Z203700', intgroup="condition", returnData=TRUE)
dwf4_1 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + theme_cowplot(font_size = 7) +labs(y = "count", x = NULL, title = 'CYP90B1-1', subtitle = 'Ceric.1Z203700')
dwf4_1 <- dwf4_1 + stat_compare_means(comparisons = my_comparisons, size = 2) 

d <- plotCounts(cdds, gene='Ceric.1Z024300', intgroup="condition", returnData=TRUE)
dwf4_2 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + theme_cowplot(font_size = 7) +labs(y = "count", x = NULL, title = 'CYP90B1-2', subtitle = 'Ceric.1Z024300')
dwf4_2 <- dwf4_2 + stat_compare_means(comparisons = my_comparisons, size = 2) 

d <- plotCounts(cdds, gene='Ceric.07G096400', intgroup="condition", returnData=TRUE)
cpd <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + theme_cowplot(font_size = 7) +labs(y = "count", x = NULL, title = 'CYP90A1', subtitle = 'Ceric.07G096400')
#cpd <- cpd + stat_compare_means(comparisons = my_comparisons, size = 2) 

d <- plotCounts(cdds, gene='Ceric.05G026200', intgroup="condition", returnData=TRUE)
rot3 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'CYP90C1', subtitle = 'Ceric.05G026200') + theme_cowplot(font_size = 7) 
rot3 <- rot3 + stat_compare_means(comparisons = my_comparisons, size = 2) 

d <- plotCounts(cdds, gene='Ceric.11G061100', intgroup="condition", returnData=TRUE)
aba1 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'ABA1', subtitle = 'Ceric.11G061100') + theme_cowplot(font_size = 7) 
aba1 <- aba1 + stat_compare_means(comparisons = my_comparisons, size = 2) 
aba1

d <- plotCounts(cdds, gene='Ceric.34G051400', intgroup="condition", returnData=TRUE)
ga3 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'CYP701/GA3', subtitle = 'Ceric.34G051400') + theme_cowplot(font_size = 7) 
ga3 <- ga3 + stat_compare_means(comparisons = my_comparisons, size = 2) 
ga3

d <- plotCounts(cdds, gene='Ceric.01G052700', intgroup="condition", returnData=TRUE)
kao <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'CYP88A3/KAO', subtitle = 'Ceric.01G052700') + theme_cowplot(font_size = 7) 
kao <- ga3 + stat_compare_means(comparisons = my_comparisons, size = 2) 
kao

d <- plotCounts(cdds, gene='Ceric.21G024600', intgroup="condition", returnData=TRUE)
ga4h <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'GA3OX/GA4H', subtitle = 'Ceric.21G024600') + theme_cowplot(font_size = 7) 
ga4h <- ga3 + stat_compare_means(comparisons = my_comparisons, size = 2) 
ga4h

d <- plotCounts(cdds, gene='Ceric.21G054300', intgroup="condition", returnData=TRUE)
her7l <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'HER7-like', subtitle = 'Ceric.21G054300') + theme_cowplot(font_size = 7) 
her7l <- ga3 + stat_compare_means(comparisons = my_comparisons, size = 2) 
her7l

d <- plotCounts(cdds, gene='Ceric.34G043600', intgroup="condition", returnData=TRUE)
gaia1 <- ggplot(d, aes(x=condition, y=count)) + geom_boxplot() + scale_y_log10() + labs(y = "count", x = NULL, title = 'GAIA1', subtitle = 'Ceric.34G043600') + theme_cowplot(font_size = 7) 
gaia1 <- ga3 + stat_compare_means(comparisons = my_comparisons, size = 2) 
gaia1


# Put it all together with cowplot
maplots <- plot_grid(parta, partb, partc, labels = 'AUTO', ncol = 1, label_size = 10)
maplots
boxplots <- plot_grid(her7, dwf4_1, cpd, NULL, dwf4_2, rot3, labels = c('D', 'F', 'H', 'E', 'G', 'I'), ncol = 3, label_size = 10)
boxplots
plot_grid(maplots, boxplots, labels = c('', ''), rel_widths = c(1, 2.5))
ggsave2('Desktop/gene_expression_jen_local/Figure4.pdf', device = "pdf", width = 16.5, height = 9, units = "cm")

