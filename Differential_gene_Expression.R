#Script for running the differential genes expression of the clover samples using edgeR###

library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(dplyr)
library(pheatmap)
library(tibble)
library(tidyr)
library(stringr)
library(viridis)
library(pheatmap)
library(mixOmics)
library(VennDiagram)
library(patchwork)
library(UpSetR)
library(ComplexUpset)
library(topGO)
library(rrvgo)

#Read complete data table

d <- read.delim ("completed.expression.txt")
rownames(d) <- d$GENE
d <- d[,-1]
d <- d[ , order(names(d))]
colnames(d)

#Importing the table with sampleInfo 
sampleInfo <- read.csv("design.csv", sep = ";")

sampleInfo$Genotype <- as.factor(sampleInfo$Genotype)
sampleInfo$Condition <- as.factor(sampleInfo$Condition)

Group <- factor(paste(sampleInfo$Genotype, sampleInfo$Condition, sep="_"))
sampleInfo <- cbind(sampleInfo,Group=Group)
sampleInfo


#Creating the DGEList
DGEList <- DGEList(d, remove.zeros = TRUE, group=Group)

DGEList$samples$Condition <- sampleInfo$Condition
DGEList$samples$Genotype <- sampleInfo$Genotype
DGEList$samples$group <- sampleInfo$Group

#Calculating presudocounts 
pseudoCounts <- log2(DGEList$counts + 1)

#Histogram
hist(pseudoCounts[ ,"S10_F_S1_1_1_To"], main = "", xlab = "pseudoCounts")

#Boxplot for pseudo-counts
par(mar = c(12,2,1,2))
boxplot(pseudoCounts, col = "gray", las = 3, cex.names = 1)


#MDS for presudo-counts------------------------------ 

colConditions <- brewer.pal(8, "Dark2")
colConditions <- colConditions[match(sampleInfo$Condition,
                                     levels(sampleInfo$Condition))]

#plotMDS(pseudoCounts, pch = pchGenotypes, col = colConditions)

plotMDS(pseudoCounts, col = colConditions)

legend("topright", lwd = 2, col = brewer.pal(8, "Dark2")[1:8], 
       legend = levels(sampleInfo$Condition))


#PlotMDS for S10 samples
colnames(pseudoCounts)
pseudoCounts_S10 <- pseudoCounts[,c(1:48)]

plotMDS(pseudoCounts_S10, col = colConditions)


#PlotMDS for Tienshan  samples
pseudoCounts_Ti <- pseudoCounts[,c(49:72)]

plotMDS(pseudoCounts_Ti, col = colConditions)

#PlotMDS for X88  samples
pseudoCounts_X88 <- pseudoCounts[,c(97:120)]

plotMDS(pseudoCounts_X88, col = colConditions)


#PlotMDS for progenitor  samples
pseudoCounts_Progenitor <- pseudoCounts[,c(73:96)]

plotMDS(pseudoCounts_Progenitor, col = colConditions)


#Filtering using the build in function in edgeR###

keep <- filterByExpr(DGEList, group = Group)
DGEList <- DGEList[keep, , keep.lib.sizes=FALSE]
table(keep)


#Calculate normalizing factors

DGEList <- calcNormFactors(DGEList, method="TMM")
DGEList$samples


#Extracting normalized counts

normCounts <- edgeR::cpm(DGEList)
pseudoNormCounts <- edgeR::cpm(DGEList, log = TRUE, prior.count = 1)

###########################################

#Creating the design matrix

design.matrix <- model.matrix(~0+Group)

rownames(design.matrix) <- colnames(DGEList)
colnames(design.matrix) <- levels(sampleInfo$Group)
design.matrix


#Creating contrasts 

my.contrasts <- makeContrasts(
  S10_To_A2 = S10_To_A2-S10_To_A1,
  S10_To_A3 = S10_To_A3-S10_To_A1,
  S10_To_A4 = S10_To_A4-S10_To_A1,
  S10_To_A1_F1 = S10_To_F1-S10_To_A1,
  S10_Tp_A2 = S10_Tp_A2-S10_Tp_A1,
  S10_Tp_A3 = S10_Tp_A3-S10_Tp_A1,
  S10_Tp_A4 = S10_Tp_A4-S10_Tp_A1,
  S10_Tp_A1_F1 = S10_Tp_F1-S10_Tp_A1,
  S10_To_F2 = S10_To_F2-S10_To_F1,
  S10_To_F3 = S10_To_F3-S10_To_F1,
  S10_To_F4 = S10_To_F4-S10_To_F1,
  S10_Tp_F2 = S10_Tp_F2-S10_Tp_F1,
  S10_Tp_F3 = S10_Tp_F3-S10_Tp_F1,
  S10_Tp_F4 = S10_Tp_F4-S10_Tp_F1,
  Tienshan_To_F2 = Tienshan_To_F2-Tienshan_To_F1,
  Tienshan_To_F3 = Tienshan_To_F3-Tienshan_To_F1,
  Tienshan_To_F4 = Tienshan_To_F4-Tienshan_To_F1,
  Tienshan_Tp_F2 = Tienshan_Tp_F2-Tienshan_Tp_F1,
  Tienshan_Tp_F3 = Tienshan_Tp_F3-Tienshan_Tp_F1,
  Tienshan_Tp_F4 = Tienshan_Tp_F4-Tienshan_Tp_F1,
  TO_F2 = TO_F2-TO_F1,
  TO_F3 = TO_F3-TO_F1,
  TO_F4 = TO_F4-TO_F1,
  TP_F2 = TP_F2-TP_F1,
  TP_F3 = TP_F3-TP_F1,
  TP_F4 = TP_F4-TP_F1,
  X88_To_F2 = X88_To_F2-X88_To_F1,
  X88_To_F3 = X88_To_F3-X88_To_F1,
  X88_To_F4 = X88_To_F4-X88_To_F1,
  X88_Tp_F2 = X88_Tp_F2-X88_Tp_F1,
  X88_Tp_F3 = X88_Tp_F3-X88_Tp_F1,
  X88_Tp_F4 = X88_Tp_F4-X88_Tp_F1,
  levels=design.matrix)

my.contrasts

#Test for DGE using the f-test (glmQLFTest)
#Estimate the dispersion and fit the model 

DGEList <- estimateDisp(DGEList, design.matrix)
fit_glmQL <- glmQLFit(DGEList, design.matrix)


colnames(my.contrasts)

DE_Genes_GLM <- NULL
DGE_Number <- NULL
for (i in (1:32)) {
  glmqlf_test <- glmQLFTest(fit_glmQL, contrast=my.contrasts[,i])
  
  #print(topTags(glmqlf_test))
  
  res_glmqlf<- topTags(glmqlf_test, n = nrow(DGEList$counts))
  selected_glmqlf <- res_glmqlf$table$FDR < 0.001 & abs(res_glmqlf$table$logFC) > 2
  selected_glmqlf <- res_glmqlf$table[selected_glmqlf, ]
  selected_glmqlf <- rownames_to_column(selected_glmqlf) %>% rename(gene_ID = rowname)
  selected_glmqlf$Group <- colnames(my.contrasts)[i]
  
  for (j in 1:nrow(selected_glmqlf))
    if (selected_glmqlf$logFC[j] > 0){
      selected_glmqlf$UD[j] <- "1"
      
    } else {
      selected_glmqlf$UD[j] <- "-1"
      
    }
  
  DE_Genes_GLM <- rbind(DE_Genes_GLM, selected_glmqlf)
  
  #
  DGE_Number$Group[i] <- colnames(my.contrasts)[i]
  DGE_Number$Number[i] <- nrow(selected_glmqlf)
  DGE_Number$Upregulated[i] <- nrow(filter(selected_glmqlf, logFC >= 2))
  DGE_Number$Downregulated[i] <- nrow(filter(selected_glmqlf, logFC <= 2))
}  

DE_Genes_GLM <- as.data.frame(DE_Genes_GLM)
DGE_Number <- as.data.frame(DGE_Number)

write.xlsx(DE_Genes_GLM, "DE_Genes_GLMqlf.xlsx")


#Select only the unregulated genes####### 

DE_Genes_GLM <- readxl::read_xlsx("DE_Genes_GLMqlf.xlsx")
DE_Genes_GLM_upregulated <- filter(DE_Genes_GLM, UD == 1)
#Selecting the genes upregulated in the ambient samples 

DE_Genes_GLM_upregulated_ambient<- filter(DE_Genes_GLM_upregulated, Group %in% c("S10_To_A2", "S10_Tp_A2", "S10_To_A3", "S10_Tp_A3", "S10_To_A4", "S10_Tp_A4", "S10_To_A1_F1", "S10_Tp_A1_F1"))


#Eliminating the ambient unregulated genes from the rest of the samples 
DE_Genes_GLM_upregulated_no_ambient <- anti_join(DE_Genes_GLM_upregulated, DE_Genes_GLM_upregulated_ambient, by="gene_ID")
Groups <- unique(DE_Genes_GLM_upregulated_no_ambient$Group)


#Adding a column for To, Tp, private or common 

DE_GLM_Up_no_amb_separated <- NULL

for (i in c(1:3, 7:9, 13:15, 19:21)){
  To <- filter(DE_Genes_GLM_upregulated_no_ambient,  Group %in% Groups[i])
  Tp <- filter(DE_Genes_GLM_upregulated_no_ambient,  Group %in% Groups[i+3])
  
  for (j in 1:nrow(To)){
    if (To$gene_ID[j] %in% Tp$gene_ID){
      To$Separation[j] <- "Common"
      
    } else {
      To$Separation[j] <- "To_private"
    }
  }
  
  for (n in 1:nrow(Tp)){
    if (Tp$gene_ID[n] %in% To$gene_ID){
      Tp$Separation[n] <- "Common"
      
    } else {
      Tp$Separation[n] <- "Tp_private"
    }
  }
  DE_GLM_Up_no_amb_separated <- rbind(DE_GLM_Up_no_amb_separated, To, Tp)
  
}


#Adding functional annotations######

Annotatios_all <- read.delim("To_Annotations.txt", sep ="\t")

names(Annotatios_all)[names(Annotatios_all) == "Gene"] <- "gene_ID"

#intersect with the differentially expressed genes
DE_GLM_Up_no_amb_separated_Annotations <- inner_join(DE_GLM_Up_no_amb_separated, Annotatios_all, by="gene_ID")


#####Generate the PCA plots####

resPCA <- pca(t(logcpm), ncomp = 6)
#resPCA <- pca(t(logcpm[, c(25:120)]), ncomp = 6)

colnames(logcpm)
plot(resPCA)

plotIndiv(resPCA, group = sampleInfo$Group)
plotIndiv(resPCA)

#write.table(sampleInfo, "sampleInfo_PCA.txt", sep="\t", quote = F)

sampleInfo <- read.delim("sampleInfo_PCA.txt", sep="\t")
sampleInfo$Genome <- as.factor(sampleInfo$Genome)

plot1 <- plotIndiv(resPCA, group = sampleInfo$Genotype1, col.per.group = colorBlindGrey8[1:5], pch = sampleInfo$Genome,
                   legend = T, legend.title = 'Samples', title = 'PCA Hylite results - All Samples',
                   style = 'ggplot2', size.xlabel = 15, size.ylabel = 15, ind.names = F, cex = 5, point.lwd = 1.5,
                   size.legend = 15, size.legend.title = 15)


#Could not change the order of the samples in the legend using just mixomics, had to use ggplot
# extract the data.frame used for plot
df <- plot1$df
labels_x <- plot1$graph$labels[[1]]
labels_y <- plot1$graph$labels[[2]]
colnames(df)

#classes <- unique(sampleInfo$Genotype1)
#df$group <- factor(df$group, levels = c("T. repens S10", "T. repens Ti", "T. repens 880",  "T. occidentale", "T. pallescens"))
df$group <- factor(df$group, levels = c( "T. occidentale", "T. pallescens", "T. repens 880", "T. repens Ti","T. repens S10"))
df$genome <- factor(df$pch.legend, levels = c("To", "Tp", "TrTo", "TrTp"))

#cols <- c("#A3CC7A","#54990F", "#1f78b4", "#a6cee3","#E6B8BF")
cols <- c("#a6cee3","#E6B8BF", "#1f78b4", "#54990F", "#A3CC7A")
names(cols) <- levels(df$group)

labels=c(substitute(paste(italic("T. occidentale"))),
         substitute(paste(italic("T. pallescens"))),
         "Tp?To 880",
         substitute(paste(italic("T. repens"), "Ti")),
         substitute(paste(italic("T. repens"), "S10")))


p1 <- ggplot(df, aes(x, y, col = group, shape=genome, fill=genome)) + 
  theme_bw() +
  geom_point(size=11, alpha=1, stroke = 2) +
  scale_shape_manual(values=c(21, 22, 24, 23), labels = c("To", "Tp", bquote('Tr'['To']), bquote('Tr'['Tp'])), name="Genome")+
  scale_fill_manual(values=c("black", "black", "lightgrey", "lightgrey"), labels = c("To", "Tp", bquote('Tr'['To']), bquote('Tr'['Tp'])), name="Genome")+
  scale_colour_manual(values = cols, labels=labels)+
  labs(x=labels_x, y=labels_y, color = "Genotype")+
  theme_classic()+
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(size=25),
        legend.key.size = unit(1.3, 'cm'), legend.text = element_text(size=25), legend.title = element_text(size=25),
        strip.text = element_text(size=25))
p1


#Plot PCA for S10 samples############

resPCA_S10 <- pca(t(logcpm[, c(1:48)]), ncomp = 6)
plot(resPCA_S10)

plotIndiv(resPCA_S10, comp = c(2,3))

sampleInfo$Genotype2 <- as.factor(paste(sampleInfo$Genotype1, sampleInfo$Condition, sep = " "))
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "Tp")] <- "Tp"
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "To")] <- "To"
sampleInfo$Subgenome <- as.factor(sampleInfo$Subgenome)

plot2 <- plotIndiv(resPCA_S10, group = sampleInfo$Condition[c(1:48)], pch=sampleInfo$Subgenome[c(1:48)],
                   col.per.group = colorBlindGrey8[c(7, 8,6, 5, 1:4)], legend = T, legend.title = 'Condition', legend.title.pch = 'Subgenome',
                   title = 'PCA plot - S10 Samples', style = 'ggplot2',size.xlabel = 20, size.ylabel = 20, cex = 5, point.lwd = 1.5,
                   size.legend = 20, size.legend.title = 20)

#Re-make the figure in ggplot

df <- plot2$df
labels_x <- plot2$graph$labels[[1]]
labels_y <- plot2$graph$labels[[2]]
colnames(df)

df$group <- factor(df$group, levels = c("A1", "A2", "A3",  "A4", "F1", "F2", "F3", "F4"),
                   labels = c("A0", "A1", "A2",  "A3", "F0", "F1", "F2", "F3"))

cols <- c("#CFE6B8","#A3CC7A","#78B33E", "#54990F",  "#B8DEE6","#7ABECC", "#3B9AB2", "#286879")


cols <- c("#c2eb9c","#94d05c","#78b33e","#1e6a00","#A8B6CC","#7Ca1CC","#3D65A5","#001c6b")
names(cols) <- levels(df$group)

p <- ggplot(df, aes(x, y, fill= group, shape=as.factor(pch.legend))) + 
  geom_point(size=11, alpha=1, stroke = 1) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=c(24,23), labels = c(bquote('Tr'['To']), bquote('Tr'['Tp'])), name="Subgenome")+
  #scale_colour_manual(values = cols)+
  labs(x=labels_x, y=labels_y, fill = "Condition")+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_classic()+
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(size=25),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25), legend.title = element_text(size=25),
        strip.text = element_text(size=25))
p


##PCA Plot for Ti samples##########################

colnames(logcpm)

resPCA_Ti <- pca(t(logcpm[, c(49:72)]), ncomp = 6)
plot(resPCA_Ti)

plotIndiv(resPCA_Ti, comp = c(2,3))

sampleInfo$Genotype2 <- as.factor(paste(sampleInfo$Genotype1, sampleInfo$Condition, sep = " "))
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "Tp")] <- "Tp"
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "To")] <- "To"
sampleInfo$Subgenome <- as.factor(sampleInfo$Subgenome)

plot2 <- plotIndiv(resPCA_Ti, group = sampleInfo$Condition[c(49:72)], pch=sampleInfo$Subgenome[c(49:72)],
                   col.per.group = colorBlindGrey8[c(7, 8,6, 5)], legend = T, legend.title = 'Condition', legend.title.pch = 'Subgenome',
                   title = 'PCA plot - S10 Samples', style = 'ggplot2',size.xlabel = 20, size.ylabel = 20, cex = 5, point.lwd = 1.5,
                   size.legend = 20, size.legend.title = 20)

#Re-make the figure in ggplot

df <- plot2$df
labels_x <- plot2$graph$labels[[1]]
labels_y <- plot2$graph$labels[[2]]
colnames(df)

df$group <- factor(df$group, levels = c("F1", "F2", "F3", "F4"),
                   labels = c("F0", "F1", "F2", "F3"))

cols <- c("#A8B6CC","#7Ca1CC","#3D65A5","#001c6b")
names(cols) <- levels(df$group)

p <- ggplot(df, aes(x, y, fill= group, shape=as.factor(pch.legend))) + 
  geom_point(size=11, alpha=1, stroke = 1) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=c(24,23), labels = c(bquote('Tr'['To']), bquote('Tr'['Tp'])), name="Subgenome")+
  #scale_colour_manual(values = cols)+
  labs(x=labels_x, y=labels_y, fill = "Condition")+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_classic()+
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(size=25),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25), legend.title = element_text(size=25),
        strip.text = element_text(size=25))
p


##PCA Plot for X88 samples############################### 
colnames(logcpm)

resPCA_X88 <- pca(t(logcpm[, c(97:120)]), ncomp = 6)
plot(resPCA_X88)

plotindiv(resPCA_X88, comp = c(2,3))

sampleInfo$Genotype2 <- as.factor(paste(sampleInfo$Genotype1, sampleInfo$CondiX88on, sep = " "))
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "Tp")] <- "Tp"
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "To")] <- "To"
sampleInfo$Subgenome <- as.factor(sampleInfo$Subgenome)

plot2 <- ploX88ndiv(resPCA_X88, group = sampleInfo$Condition[c(97:120)], pch=sampleInfo$Subgenome[c(97:120)],
                   col.per.group = colorBlindGrey8[c(7, 8,6, 5)], legend = T, legend.title = 'Condition', legend.title.pch = 'Subgenome',
                   title = 'PCA plot - S10 Samples', style = 'ggplot2',size.xlabel = 20, size.ylabel = 20, cex = 5, point.lwd = 1.5,
                   size.legend = 20, size.legend.title = 20)

#Re-make the figure in ggplot

df <- plot2$df
labels_x <- plot2$graph$labels[[1]]
labels_y <- plot2$graph$labels[[2]]
colnames(df)

df$group <- factor(df$group, levels = c("F1", "F2", "F3", "F4"),
                   labels = c("F0", "F1", "F2", "F3"))

cols <- c("#A8B6CC","#7Ca1CC","#3D65A5","#001c6b")
names(cols) <- levels(df$group)

p <- ggplot(df, aes(x, y, fill= group, shape=as.factor(pch.legend))) + 
  geom_point(size=11, alpha=1, stroke = 1) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=c(24,23), labels = c(bquote('Tr'['To']), bquote('Tr'['Tp'])), name="Subgenome")+
  #scale_colour_manual(values = cols)+
  labs(x=labels_x, y=labels_y, fill = "Condition")+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_classic()+
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(size=25),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25), legend.title = element_text(size=25),
        strip.text = element_text(size=25))
p


##PCA Plot for Progenitors##########
colnames(logcpm)

resPCA_Prog <- pca(t(logcpm[, c(73:96)]), ncomp = 6)
plot(resPCA_Prog)

plotIndiv(resPCA_Prog, comp = c(2,3))

sampleInfo$Genotype2 <- as.factor(paste(sampleInfo$Genotype1, sampleInfo$Condition, sep = " "))
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "Tp")] <- "Tp"
sampleInfo$Subgenome[str_which(sampleInfo$Genotype, "To")] <- "To"
sampleInfo$Subgenome <- as.factor(sampleInfo$Subgenome)

plot2 <- plotIndiv(resPCA_Prog, group = sampleInfo$Condition[c(73:96)], pch=sampleInfo$Genome[c(73:96)],
                   col.per.group = colorBlindGrey8[c(7, 8,6, 5)], legend = T, legend.title = 'Condition', legend.title.pch = 'Subgenome',
                   title = 'PCA plot - S10 Samples', style = 'ggplot2',size.xlabel = 20, size.ylabel = 20, cex = 5, point.lwd = 1.5,
                   size.legend = 20, size.legend.title = 20)

#Re-make the figure in ggplot

df <- plot2$df
labels_x <- plot2$graph$labels[[1]]
labels_y <- plot2$graph$labels[[2]]
colnames(df)

df$group <- factor(df$group, levels = c("F1", "F2", "F3", "F4"),
                   labels = c("F0", "F1", "F2", "F3"))

cols <- c("#A8B6CC","#7Ca1CC","#3D65A5","#001c6b")
names(cols) <- levels(df$group)

p <- ggplot(df, aes(x, y, fill= group, shape=as.factor(pch.legend))) + 
  geom_point(size=11, alpha=1, stroke = 1) +
  scale_fill_manual(values=cols) +
  scale_shape_manual(values=c(24,23), labels = c("To", "Tp"), name="Genome")+
  #scale_colour_manual(values = cols)+
  labs(x=labels_x, y=labels_y, fill = "Condition")+
  guides(fill = guide_legend(override.aes = list(shape = 21)))+
  theme_classic()+
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(size=25),
        legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25), legend.title = element_text(size=25),
        strip.text = element_text(size=25))
p


#Plotting the figure with the number of upregulated and downregulated genes in each condition for all the genotypes###

data <- read.delim("Data_Up_Down_Genes.txt")

p1<-ggplot(data) +
  aes(x = Condition1, fill = Subgenome3, weight = value) +
  geom_bar() +
  labs(
    x=" ",
    y = "Number of differentially expressed transcripts"
  ) + 
  theme_classic()+
  theme(axis.text.x=element_text(hjust=1,vjust=0.5, size = 30), axis.title.y = element_text(size=25, face = "bold"), axis.text.y = element_text(size=22),
        legend.title = element_text(size=22), legend.text = element_text(size=20), strip.text = element_text(size=30, face="bold"),
        panel.spacing = unit(0.2, "lines"))+
  facet_grid(cols = vars(reorder(genotype, Order)),  scales = "free", space = "free")+
  geom_text(aes(label = gsub("-", "", ..count..)), stat = "count", vjust = ifelse(data$value >= 0, 1, -0.25), colour = "black", size=9)+
  scale_fill_manual(name = "Expression change", labels = c("To Downregulated", "To Upregulated", "Tp Downregulated", "Tp Upregulated"),
                    values = brewer.pal(4,"Paired"))+ scale_y_continuous(breaks = pretty(data$value), labels = abs(pretty(data$value)))
p1


####Plotting the percentage of genes in each of the TO private, TP private and common groups

Group_separation <- read.delim("Group_separation_Fig.txt")

Group_separation$type <- factor(Group_separation$type, levels = c("To_Private", "Common", "Tp_Private"))

p1<- ggplot(data=Group_separation, aes(y=percent, x=Condition, fill= type)) +
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values=c("darkolivegreen3", "slategray2", "forestgreen"),
                    labels = c(bquote('Specific Tr'['To']), bquote('Common Tr'['To']*'Tr'['Tp']), bquote('Specific Tr'['Tp'])))+
  geom_text(aes(label=percent1), color="black", size=9, position = position_stack(vjust = 0.5))+
  ggtitle(" ")+
  theme_classic() + xlab(" ") + ylab("Percent upregulated transcripts") +
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(25),
        legend.title = element_blank(), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25),
        strip.text = element_text(size=25,face="bold"), panel.spacing = unit(0.2, "lines"))+
  facet_grid(cols=vars(reorder(Genotype1, Order)),  scales = "free", space = "free")
p1




###Upset plots####
###Generate the Upset plots comparing the upregulated genes between the TO and its derived subgenomes

test <- NULL
test[["S10_To"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "S10_To"),]$gene_ID)
test[["Ti_To"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "Tienshan_To"),]$gene_ID)
test[["X88_To"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "X88_To"),]$gene_ID)
test[["TO"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "TO"),]$gene_ID)

#Complex Upset
set_size(5, 3.6)

data <- fromList(test)
p1 <- ComplexUpset::upset(data, c("TO","X88_To", "Ti_To", "S10_To"), sort_sets=FALSE, stripes = "white",guides='over',
                          height_ratio=0.25, width_ratio = 0.3,
                          labeller =function(labels) {c("")},
                          themes=upset_default_themes(element_blank()),
                          base_annotations=list(
                            'Intersection size'=(
                              intersection_size(text=list(size=7),width=0.8,
                                                bar_number_threshold=1  # show all numbers on top of bars 
                              ))),
                          matrix=(
                            intersection_matrix(geom=geom_point(shape=21, size=12, stroke=0.30),
                                                segment=geom_segment(
                                                  linetype='solid', size=1.5, color="#666666",
                                                ))
                            + scale_color_manual(
                              values=c('S10_To'="#A3CC7A", 'Ti_To'="#54990F", 'X88_To'= "#1f78b4", 'TO'="#a6cee3"),
                              labels=c(substitute(paste(italic("T. repens"), "S10")),
                                       substitute(paste(italic("T. repens"), "Ti")),
                                       "Tp?To 880",
                                       substitute(paste(italic("T. occidentale")))),
                              guide=guide_legend(override.aes=list(shape='circle'))
                            )
                          ),
                          queries=list(
                            upset_query(set='S10_To', fill="#A3CC7A"),
                            upset_query(set='Ti_To', fill="#54990F"),
                            upset_query(set='X88_To', fill="#1f78b4"),
                            upset_query(set='TO', fill="#a6cee3")
                          )
) & theme(legend.margin=margin(t=3.5, r=0, b=3.5, l=2.5, 'cm'), legend.text = element_text(size=25),
          axis.text = element_text(size=20), axis.title =element_text(size=25))
p1


svg("Figures/Upsetplot_TO_samples.svg", width = 16, height = 10)
p1
dev.off()


###Generate the Upset plots comparing the upregulated genes between the TP and its derived subgenomes

test <- NULL
test[["S10_Tp"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "S10_Tp"),]$gene_ID)
test[["Ti_Tp"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "Tienshan_Tp"),]$gene_ID)
test[["X88_Tp"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "X88_Tp"),]$gene_ID)
test[["TP"]] <- unique(DE_Genes[str_which(DE_Genes$Group, "TP"),]$gene_ID)

#upsetR
data<- fromList(test)
upset(data, order.by = "freq", sets=c("S10_Tp","Ti_Tp","X88_Tp","TP"), keep.order=T,	nsets=4, text.scale = c(3, 3, 3, 3, 2, 3),  sets.bar.color=c("#bd0026","#0c2c84", "#005a32", "#F2AD00"), point.size = 4)

data <- fromList(test)
p1 <- ComplexUpset::upset(data,c("TP","X88_Tp", "Ti_Tp", "S10_Tp"), sort_sets=FALSE, stripes = "white",guides='over',
                          height_ratio=0.25, width_ratio = 0.3,
                          labeller =function(labels) {c("")},
                          themes=upset_default_themes(element_blank()),
                          base_annotations=list(
                            'Intersection size'=(
                              intersection_size(text=list(size=7),width=0.8,
                                                bar_number_threshold=1  # show all numbers on top of bars 
                              ))),
                          matrix=(
                            intersection_matrix(geom=geom_point(shape=21, size=12, stroke=0.30),
                                                segment=geom_segment(
                                                  linetype='solid', size=1.5, color="#666666",
                                                ))
                            + scale_color_manual(
                              values=c('S10_Tp'="#A3CC7A", 'Ti_Tp'="#54990F", 'X88_Tp'= "#1f78b4", 'TP'="#a6cee3"),
                              labels=c(substitute(paste(italic("T. repens"), "S10")),
                                       substitute(paste(italic("T. repens"), "Ti")),
                                       "Tp?To 880",
                                       substitute(paste(italic("T. pallescens")))),
                              guide=guide_legend(override.aes=list(shape='circle'))
                            )
                          ),
                          queries=list(
                            upset_query(set='S10_Tp', fill="#A3CC7A"),
                            upset_query(set='Ti_Tp', fill="#54990F"),
                            upset_query(set='X88_Tp', fill="#1f78b4"),
                            upset_query(set='TP', fill="#a6cee3")
                          )
) & theme(legend.margin=margin(t=3.5, r=0, b=3.5, l=2.5, 'cm'), legend.text = element_text(size=25),
          axis.text = element_text(size=25), axis.title =element_text(size=25))
p1


svg("Figures/Upsetplot_TP_samples.svg", width = 16, height = 10)
p1
dev.off()



#####Separating the downregulated genes###

DE_Genes_GLM <- readxl::read_xlsx("DE_Genes_GLMqlf.xlsx")
DE_Genes_GLM_down <- filter(DE_Genes_GLM, UD == -1)

#Selecting the genes downregulated in the ambient samples 

DE_Genes_GLM_down_ambient<- filter(DE_Genes_GLM_down, Group %in% c("S10_To_A2", "S10_Tp_A2", "S10_To_A3", "S10_Tp_A3", "S10_To_A4", "S10_Tp_A4", "S10_To_A1_F1", "S10_Tp_A1_F1"))

#Eliminating the ambient unregulated genes from the rest of the samples 
DE_Genes_GLM_down_no_ambient <- anti_join(DE_Genes_GLM_down, DE_Genes_GLM_down_ambient, by="gene_ID")
Groups <- unique(DE_Genes_GLM_down_no_ambient$Group)

Down_Genes_Separated <- NULL

for (i in c(1:3, 7:9, 13:15, 19:21)){
  To <- filter(DE_Genes_GLM_down_no_ambient,  Group %in% Groups[i])
  Tp <- filter(DE_Genes_GLM_down_no_ambient,  Group %in% Groups[i+3])
  To$Separation <- "nimic"
  Tp$Separation <- "nimic"
  
  for (j in 1:nrow(To)){
    if (To$gene_ID[j] %in% Tp$gene_ID){
      To$Separation[j] <- "Common"
      
    } else {
      To$Separation[j] <- "To_private"
    }
  }
  
  for (n in 1:nrow(Tp)){
    if (Tp$gene_ID[n] %in% To$gene_ID){
      Tp$Separation[n] <- "Common"
      
    } else {
      Tp$Separation[n] <- "Tp_private"
    }
  }
  Down_Genes_Separated <- rbind(Down_Genes_Separated, To, Tp)
  
}

Down_Genes_separated <- as.data.frame(Down_Genes_separated)


#plot the number of genes in the TO private, TP private and the common group

Down_Genes_Separated$Group_separation <- paste0(Down_Genes_Separated$Group, "_", Down_Genes_Separated$Separation)

Down_separation <- as.data.frame(table(Down_Genes_Separated$Group_separation))

Down_separation$type <- "empty"
Down_separation[str_which(Down_separation$Var1, "To_private"),]$type <- "To_private"
Down_separation[str_which(Down_separation$Var1, "Tp_private"),]$type <- "Tp_private"
Down_separation[str_which(Down_separation$Var1, "Common"),]$type <- "Common"

Down_separation$Condition <- "empty"
Down_separation[str_which(Down_separation$Var1, "F1"),]$Condition <- "F0"
Down_separation[str_which(Down_separation$Var1, "F2"),]$Condition <- "F1"
Down_separation[str_which(Down_separation$Var1, "F3"),]$Condition <- "F2"
Down_separation[str_which(Down_separation$Var1, "F4"),]$Condition <- "F3"
Down_separation$Condition <- as.factor(Down_separation$Condition)


Down_separation$Genotype <- "empty"
Down_separation[str_which(Down_separation$Var1, "S10"),]$Genotype <- "T. repens S10"
Down_separation[str_which(Down_separation$Var1, "Ti"),]$Genotype <- "T. repens Ti"
Down_separation[str_which(Down_separation$Var1, "X88"),]$Genotype <- "Tp?To 880"
Down_separation[str_which(Down_separation$Var1, "TO|TP"),]$Genotype <- "Progenitors"

Down_separation$type <- factor(Down_separation$type, levels = c("To_private", "Common", "Tp_private"))
Down_separation$Genotype <- factor(Down_separation$Genotype, levels = c("T. repens S10", "T. repens Ti", "Tp?To 880","Progenitors"))


Down_separation$Condition1 <-  paste0(Down_separation$Genotype, "_", Down_separation$Condition)

Down_separation <-  Down_separation[-c(7,9,11,19,21,23,31,33,35,43,45,47),]

Down_separation <-Down_separation %>% 
  group_by(Condition1) %>%
  mutate(per = 100 *Freq/sum(Freq)) %>%
  mutate(per1=paste0(round(100 *Freq/sum(Freq),1),'%'))%>%ungroup


p1<- ggplot(data=Down_separation, aes(y=per, x=Condition, fill= type)) +
  geom_bar(stat="identity", position = "stack")+
  scale_fill_manual(values=c("darkolivegreen3", "slategray2", "forestgreen"),
                    labels = c(bquote('Specific Tr'['To']), bquote('Common Tr'['To']*'Tr'['Tp']), bquote('Specific Tr'['Tp'])))+
  geom_text(aes(label=per1), color="black", size=9, position = position_stack(vjust = 0.5))+
  ggtitle(" ")+
  theme_classic() + xlab(" ") + ylab("Percent downregulated transcripts") +
  theme(axis.text.x=element_text(size=25),axis.title=element_text(size=25,face="bold"), axis.text.y = element_text(25),
        legend.title = element_blank(), legend.key.size = unit(1, 'cm'), legend.text = element_text(size=25),
        strip.text = element_text(size=25,face="bold"), panel.spacing = unit(0.2, "lines"))+
  facet_grid(cols=vars(Genotype),  scales = "free", space = "free")
p1

png("Figures/Downregulated_separated_groups.png", width = 8000, height = 4500, res = 400)
p1
dev.off()


###Generate the Upset plots comparing the downregulated genes between the TO and its derived subgenomes

test <- NULL
test[["S10_To"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "S10_To"),]$gene_ID)
test[["Ti_To"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "Tienshan_To"),]$gene_ID)
test[["X88_To"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "X88_To"),]$gene_ID)
test[["TO"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "TO"),]$gene_ID)

#upsetR
data<- fromList(test)
UpSetR::upset(data, order.by = "freq", sets=c("S10_To","Ti_To","X88_To","TO"), keep.order=T,	nsets=4, text.scale = c(3, 3, 3, 3, 2, 3),  sets.bar.color=c("#bd0026","#0c2c84", "#005a32", "#F2AD00"), point.size = 4)

#Complex Upset
set_size(5, 3.6)

data <- fromList(test)
p1 <- ComplexUpset::upset(data, c("TO","X88_To", "Ti_To", "S10_To"), sort_sets=FALSE, stripes = "white",guides='over',
                          height_ratio=0.25, width_ratio = 0.3,
                          labeller =function(labels) {c("")},
                          themes=upset_default_themes(element_blank()),
                          base_annotations=list(
                            'Intersection size'=(
                              intersection_size(text=list(size=7),width=0.8,
                                                bar_number_threshold=1  # show all numbers on top of bars 
                              ))),
                          matrix=(
                            intersection_matrix(geom=geom_point(shape=21, size=12, stroke=0.30),
                                                segment=geom_segment(
                                                  linetype='solid', size=1.5, color="#666666",
                                                ))
                            + scale_color_manual(
                              values=c('S10_To'="#A3CC7A", 'Ti_To'="#54990F", 'X88_To'= "#1f78b4", 'TO'="#a6cee3"),
                              labels=c(substitute(paste(italic("T. repens"), "S10")),
                                       substitute(paste(italic("T. repens"), "Ti")),
                                       "Tp?To 880",
                                       substitute(paste(italic("T. occidentale")))),
                              guide=guide_legend(override.aes=list(shape='circle'))
                            )
                          ),
                          queries=list(
                            upset_query(set='S10_To', fill="#A3CC7A"),
                            upset_query(set='Ti_To', fill="#54990F"),
                            upset_query(set='X88_To', fill="#1f78b4"),
                            upset_query(set='TO', fill="#a6cee3")
                          )
) & theme(legend.margin=margin(t=3.5, r=0, b=3.5, l=2.5, 'cm'), legend.text = element_text(size=25),
          axis.text = element_text(size=20), axis.title =element_text(size=25))
p1

svg("Figures/Upsetplot_TO_samples_Downregulated.svg", width = 16, height = 10)
p1
dev.off()


###Generate the Upset plots comparing the downregulated genes between the TP and its derived subgenomes

test <- NULL
test[["S10_Tp"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "S10_Tp"),]$gene_ID)
test[["Ti_Tp"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "Tienshan_Tp"),]$gene_ID)
test[["X88_Tp"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "X88_Tp"),]$gene_ID)
test[["TP"]] <- unique(Down_Genes_Separated[str_which(Down_Genes_Separated$Group, "TP"),]$gene_ID)

#upsetR
data<- fromList(test)
UpSetR::upset(data, order.by = "freq", sets=c("S10_Tp","Ti_Tp","X88_Tp","TP"), keep.order=T,	nsets=4, text.scale = c(3, 3, 3, 3, 2, 3),  sets.bar.color=c("#bd0026","#0c2c84", "#005a32", "#F2AD00"), point.size = 4)

data <- fromList(test)
p1 <- ComplexUpset::upset(data,c("TP","X88_Tp", "Ti_Tp", "S10_Tp"), sort_sets=FALSE, stripes = "white",guides='over',
                          height_ratio=0.25, width_ratio = 0.3,
                          labeller =function(labels) {c("")},
                          themes=upset_default_themes(element_blank()),
                          base_annotations=list(
                            'Intersection size'=(
                              intersection_size(text=list(size=7),width=0.8,
                                                bar_number_threshold=1  # show all numbers on top of bars 
                              ))),
                          matrix=(
                            intersection_matrix(geom=geom_point(shape=21, size=12, stroke=0.30),
                                                segment=geom_segment(
                                                  linetype='solid', size=1.5, color="#666666",
                                                ))
                            + scale_color_manual(
                              values=c('S10_Tp'="#A3CC7A", 'Ti_Tp'="#54990F", 'X88_Tp'= "#1f78b4", 'TP'="#a6cee3"),
                              labels=c(substitute(paste(italic("T. repens"), "S10")),
                                       substitute(paste(italic("T. repens"), "Ti")),
                                       "Tp?To 880",
                                       substitute(paste(italic("T. pallescens")))),
                              guide=guide_legend(override.aes=list(shape='circle'))
                            )
                          ),
                          queries=list(
                            upset_query(set='S10_Tp', fill="#A3CC7A"),
                            upset_query(set='Ti_Tp', fill="#54990F"),
                            upset_query(set='X88_Tp', fill="#1f78b4"),
                            upset_query(set='TP', fill="#a6cee3")
                          )
) & theme(legend.margin=margin(t=3.5, r=0, b=3.5, l=2.5, 'cm'), legend.text = element_text(size=25),
          axis.text = element_text(size=25), axis.title =element_text(size=25))
p1

svg("Figures/Upsetplot_TP_samples_Downregulated.svg", width = 16, height = 10)
p1
dev.off()



#GO enrichment##########################
library(topGO)
library(rrvgo)

geneID2GO <- readMappings("GOMAP_terms2.txt")
#names(geneID2GO)

all.genes <- names(geneID2GO)
all.genes <- all.genes[-1]

Gene_enrichment <- function(gene_list, group){
  
  geneList <- ifelse(all.genes %in% unique(gene_list), 1, 0)
  names(geneList) <- all.genes
  
  table(geneList)
  
  # Create topGOdata object
  GOdata <- new("topGOdata",
                ontology = "BP", # use biological process ontology
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO)
  
  
  resultFisher_weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  results <- GenTable(GOdata, Fisher = resultFisher_weight01, topNodes = 200, numChar = 100) %>% 
    mutate(Fisher = gsub("<", "", Fisher)) %>%
    filter(as.numeric(Fisher)<=0.05)
  
  results <- results[results$Annotated >1,]
  
  writexl::write_xlsx(results, paste0( group, ".xlsx"))
  return(results)
}

Gene_enrichment_fig <- function(enrichment_result, name){
  simMatrix <- calculateSimMatrix(enrichment_result$GO.ID,
                                  orgdb="org.At.tair.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(as.numeric(enrichment_result$Fisher)), enrichment_result$GO.ID)
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.At.tair.db")
  
  pdf(paste0(name, ".pdf"), width = 12, height =8)
  treemapPlot(reducedTerms)
  dev.off()
  
  #7 and 7 for the module figure
  
  return(treemapPlot(reducedTerms))
}


#GO enrichment for all the upregulated genes###

results <- Gene_enrichment(unique(Up_genes$gene_ID), "DE_Up_Genes")
Gene_enrichment_fig(results, "DE_Up_Genes")


#GO enrichment for all the downregulated genes###

results <- Gene_enrichment(unique(Down_genes$gene_ID), "GO_Down_Genes")
Gene_enrichment_fig(results, "GO_Down_Genes")

