library(pheatmap)
library(stringr)
library(tidyr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggplotify)
library(tidyverse)
library(mixOmics)
library(hrbrthemes)
library(viridis)
library(UpSetR)
library(corrplot)

#Script for normalizing the count for the clover TO and TP subgenomes together###

#Read the complete data table

d <- read.delim ("../completed.expression.txt")
rownames(d) <- d$GENE
d <- d[,-1]
d <- d[ , order(names(d))]
colnames(d)

dim(d)

#reorder
d_WC <- d[, c(1:72, 97:120)]
colnames(d_WC)

samples <- str_sub(colnames(d_WC), 1, -4)
samples <- unique(samples)

#add up the counts for the two subgenomes in the white clovers
d_WC_sum <- NULL
for (i in samples){
  print(i)
  d_WC_sum[[i]] <- d[, str_which(colnames(d), paste0(i, "_To"))] + d[, str_which(colnames(d), paste0(i, "_Tp"))]
}

d_WC_sum_df <- data.frame(d_WC_sum)
rownames(d_WC_sum_df) <- rownames(d_WC)
colnames(d_WC_sum_df)

#bind the progentitor samples
d_WC_sum_df_prog <- cbind(d_WC_sum_df, d[, c(73:96)])

colnames(d_WC_sum_df_prog)


colnames(d_WC_sum_df_prog)
groups <- c(rep("S10_A1", 3), rep("S10_A2", 3), rep("S10_A3", 3), rep("S10_A4", 3), rep("S10_F1", 3)
            , rep("S10_F2", 3), rep("S10_F3", 3), rep("S10_F4", 3), rep("Ti_F1", 3), rep("Ti_F2", 3),
            rep("Ti_F3", 3), rep("Ti_F4", 3), rep("X88_F1", 3), rep("X88_F2", 3), rep("X88_F3", 3)
            , rep("X88_F4", 3), rep("TO_F1", 3), rep("TO_F2", 3), rep("TO_F3", 3), rep("TO_F4", 3),
            rep("TP_F1", 3), rep("TP_F2", 3), rep("TP_F3", 3), rep("TP_F4", 3))


#normalise data using TMM
{
  ## TMM ( edgeR)
  library(edgeR)
  # make expression object
  ncol(d_WC_sum_df_prog)
  colnames(d_WC_sum_df_prog)
  #d.m.for.norm <- d.m[,c(9:16,65:88,91:118)]
  y <- DGEList(counts=d_WC_sum_df_prog)
  keep <- filterByExpr(y, group = groups)
  y <- y[keep, , keep.lib.sizes=FALSE]
  table(keep)
  y <- calcNormFactors(y)
  y$samples
  #plotMDS(y)
  d.m.norm.1 <- edgeR::cpm(y)
  d.m.norm <- as.data.frame(d.m.norm.1) 
  head(d_WC_sum_df_prog)
  colnames(d_WC_sum_df_prog)
  head(d.m.norm)
  colnames(d.m.norm)
}

test <- t(log(1+d.m.norm))
resPCA <- pca(test, ncomp = 6)
plotIndiv(resPCA)



###Multiply with the normalization factors

y$samples

norm.factors <- y$samples$norm.factors

norm.facs.sub <- c(rep(norm.factors[1:48],each=2), 2*norm.factors[49:72])

length(norm.facs.sub)
ncol(d)
colnames(d)
d.reordered <- d[,c(1:72,97:120,73:96)]
colnames(d.reordered)

#devide the inital values by the normalization factors
d.norm <- data.frame(t(t(d.reordered)/norm.facs.sub))
d.norm <- d.norm[rownames(d.m.norm),]

colnames(d.norm)
summary(d.norm)

#Calculate the mean normalized counts for the replicates####

d.norm1 <-  as.data.frame(t(d.norm))
rownames(d.norm1)
groups_WC <- paste(str_sub(rownames(d.norm1)[1:96], 1, -6), rep(c("To", "Tp"), 48), sep="_")
groups_prog <- str_sub(rownames(d.norm1)[97:120], 1, -3)
Groups <- c(groups_WC, groups_prog)


d.norm1$Groups <- Groups

d.norm2 <- d.norm1 %>% 
  group_by(Groups) %>%
  summarise(across(everything(), mean))


d.norm3 <-  as.data.frame(t(d.norm2))
colnames(d.norm3) <- d.norm3[1,]
d.norm3 <- d.norm3[-1,]

colnames(d.norm3)

d.norm3[,] = apply(d.norm3[,], 2, function(x) as.numeric(as.character(x)))

d.norm3_1 <- as.data.frame(t(d.norm3))


### Make table with log subgenome ratios for all samples###############

colnames(d.norm)

i=1
ToTp <- NULL
for (i in seq(1, 95, by=2)){
  ToTp[[i]] <- log((d.norm[,i]+1)/(d.norm[,i+1]+1))
}

ToTp_Ratios<- as.data.frame(do.call(cbind, ToTp))
head(ToTp_Ratios)
colnames(ToTp_Ratios) <- colnames(d.norm[seq(1, 95, by=2)])
rownames(ToTp_Ratios) <- rownames(d.norm)


ToTp.prog <- NULL
for (i in 97:108){
  ToTp.prog[[i]] <- log((d.norm[,i]+1)/(d.norm[,i+12]+1))
}

prog.ratios <- as.data.frame(do.call(cbind, ToTp.prog))
colnames(prog.ratios) <- colnames(d.norm[,97:108])
rownames(prog.ratios) <- rownames(d.norm)

head(prog.ratios)


all.ratios_split <- cbind(ToTp_Ratios, prog.ratios)
all.ratios_split.cor <- cor(all.ratios_split)
heatmap(all.ratios_split.cor)

colnames(all.ratios_split)


ratios_split1 <- as.data.frame(t(all.ratios_split))


groups <- c(rep("S10_A1", 3), rep("S10_A2", 3), rep("S10_A3", 3), rep("S10_A4", 3), rep("S10_F1", 3)
            , rep("S10_F2", 3), rep("S10_F3", 3), rep("S10_F4", 3), rep("Ti_F1", 3), rep("Ti_F2", 3),
            rep("Ti_F3", 3), rep("Ti_F4", 3), rep("X88_F1", 3), rep("X88_F2", 3), rep("X88_F3", 3)
            , rep("X88_F4", 3), rep("TP_TO_F1", 3), rep("TP_TO_F2", 3), rep("TP_TO_F3", 3), rep("TP_TO_F4", 3))


ratios_split1$Groups <- groups

#Calculate the mean ratios for the three replicates
ratios_split2 <- ratios_split1 %>% 
  group_by(Groups) %>%
  summarise(across(everything(), mean))

ratios_split3 <-  as.data.frame(t(ratios_split2))
colnames(ratios_split3) <- ratios_split3[1,]
ratios_split3 <- ratios_split3[-1,]
ratios_split3[,] = apply(ratios_split3[,], 2, function(x) as.numeric(as.character(x)))


###Filter the genes that were unregulated 

ratios_split3_DE <- ratios_split3[unique(DE_upregulated_genes$gene_ID),]

#DE Ratios correlations figure############

Ratios_DE_cor <- cor(ratios_split3_DE[,c(1:20)], method = "spearman")

Ratios_DE_cor <- Ratios_DE_cor[c(1:8, 13:20, 9:12),c(1:8, 13:20, 9:12)]


corrplot(Ratios_DE_cor, type = "upper",  method="ellipse", tl.col = "black", tl.srt = 45, addCoef.col = 'black', col.lim = c(0,1), col = COL2('RdBu'))


####Calculate the variance in log ratios across the frost conditions############

variance <- NULL

colnames(ratios_split3)

var.func.log.ratios <- function(x) { var( x[1:4])}
variance$S10_A <- apply(ratios_split3, 1, var.func.log.ratios )

var.func.log.ratios <- function(x) { var( x[5:8])}
variance$S10_F <- apply(ratios_split3, 1, var.func.log.ratios )

var.func.log.ratios <- function(x) { var( x[9:12])}
variance$prog <- apply(ratios_split3, 1, var.func.log.ratios )

var.func.log.ratios <- function(x) { var( x[13:16])}
variance$Ti <- apply(ratios_split3, 1, var.func.log.ratios )

var.func.log.ratios <- function(x) { var( x[17:20])}
variance$X88 <- apply(ratios_split3, 1, var.func.log.ratios )

variance <- as.data.frame(variance)

as.data.frame(colMeans(variance))

as.data.frame(colMeans(variance_DE))

sd_cutoff = 1.5

var.func.log.ratios <- function(x) { var( x[5:20])}
variance$All <- apply(ratios_split3, 1, var.func.log.ratios )

for (i in colnames(variance)[1:6]){
  # set total log ratio variance cutoff based on 
  hist(log(as.numeric(1+variance[,i])))
  mean.log.var <- mean(log(1+variance[,i]))
  mean.log.var
  sd.log.var <- sd(log(1+variance[,i]), na.rm=TRUE)
  sd.log.var
  var.cutoff<- mean.log.var + sd_cutoff*sd.log.var
  var.cutoff
  abline(v=var.cutoff, col="red")
  
  var.outliers <- subset(variance, log(X88) > var.cutoff)
  var.outliers.filter <- log(variance[,i]) > var.cutoff
  variance[,paste0(i, "outliers")] <- var.outliers.filter
  head(subset(variance, paste0(i, "outliers") ==TRUE))
}

nrow(var.outliers.X88)
head(var.outliers.X88)

#Variance for the DE genes

variance_DE <- variance[unique(DE_Gene_F1$gene_ID), ]

#Generate smooth plots for variance outliers and the other genes#####

variance_out_S10 <-  variance[variance$S10_Foutliers=="TRUE",]
variance_out_Ti <-  variance[variance$Tioutliers=="TRUE",]
variance_out_880 <-  variance[variance$X88outliers=="TRUE",]
variance_out_prog <-  variance[variance$progoutliers=="TRUE",]

variance_out_S10_DE <-  variance_DE[variance_DE$S10_Foutliers=="TRUE",]
variance_out_Ti_DE <-  variance_DE[variance_DE$Tioutliers=="TRUE",]
variance_out_880_DE <-  variance_DE[variance_DE$X88outliers=="TRUE",]
variance_out_prog_DE <-  variance_DE[variance_DE$progoutliers=="TRUE",]


####Smooth plot for all genes vs outliers####

# set colors
colors <- viridis(n=30,begin=0, end=1)
palette <- colorRampPalette(colors)

colnames(ratios_split3)

#outliers
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_S10), c(5:8)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_Ti), c(9:12)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_prog), c(13:16)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_880), c(17:20)]

variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_S10_DE), c(5:8)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_Ti_DE), c(9:12)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_prog_DE), c(13:16)]
variance_out_ratios <- ratios_split3[rownames(ratios_split3) %in% rownames(variance_out_880_DE), c(17:20)]


colnames(variance_out_ratios) <-  c("F0", "F1", "F2", "F3")

variance_out_ratios_long <- pivot_longer(variance_out_ratios, cols = c(1:4))

variance_out_ratios_long$Condition1 <- "empty"
variance_out_ratios_long[str_which(variance_out_ratios_long$name, "F0"),]$Condition1 <- "1"
variance_out_ratios_long[str_which(variance_out_ratios_long$name, "F1"),]$Condition1 <- "2"
variance_out_ratios_long[str_which(variance_out_ratios_long$name, "F2"),]$Condition1 <- "3"
variance_out_ratios_long[str_which(variance_out_ratios_long$name, "F3"),]$Condition1 <- "4"

lr.plot.rand.out <- runif(length(variance_out_ratios_long$name), min=-0.2, max=0.2)
variance_out_ratios_long$Condition2 <- as.numeric(variance_out_ratios_long$Condition1)+lr.plot.rand.out

png("Figures/Smooth_Ratios_Outliers_880_DE.png", res = 600, width = 4800, height =3600)
smoothScatter(variance_out_ratios_long[, c(4,2)], colramp=palette, ylim=c(-8, 8), nrpoints = 0,
              ylab = bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"), xlab=" ")
dev.off()

###All genes
#5:8 - S10 F
#9:12 - Ti
#13:16 - Prog
#17:20 - 880

#All
all_genes_test <- pivot_longer(ratios_split3[, c(17:20)], cols=c(1:4))

#DE
all_genes_test <- pivot_longer(ratios_split3_DE[, c(5:8)], cols=c(1:4))

all_genes_test$Condition1 <- "nimic"
all_genes_test[str_which(all_genes_test$name, "F1"),]$Condition1 <- "1"
all_genes_test[str_which(all_genes_test$name, "F2"),]$Condition1 <- "2"
all_genes_test[str_which(all_genes_test$name, "F3"),]$Condition1 <- "3"
all_genes_test[str_which(all_genes_test$name, "F4"),]$Condition1 <- "4"

lr.plot.rand.out <- runif(length(all_genes_test$name), min=-0.2, max=0.2)
all_genes_test$Condition2 <- as.numeric(all_genes_test$Condition1)+lr.plot.rand.out


png("Figures/Smooth_Ratios_Allgenes_S10_DE.png", res = 600, width = 4800, height =3600)
smoothScatter(all_genes_test[, c(4,2)], colramp=palette, ylim=c(-8, 8), nrpoints = 0, 
              ylab = bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"), xlab=" ")
dev.off()


###Generate the heatmap with the genes from the raffinose pathway###

d.norm4 <-d.norm3[, c(1, 3, 5, 7, 2, 4, 6,8, 9, 11, 13, 15, 10, 12, 14, 16, 17, 19, 21, 23, 18, 20, 22, 24, 25:28, 29:32, 33, 35, 37, 39, 34, 36, 38, 40)]
colnames(d.norm4)

log_d.norm4 <- log2(1+ d.norm4)

colnames(log_d.norm4)

#reorder
colnames(log_d.norm4)
log_d.norm4 <- log_d.norm4[, c(1:24,33:40, 25:32)]

####All the galactios synthases in the clover genome#####

GolS <- c("jg35812.t1", "jg44037.t1", "jg70060.t1")

#####Raffinose-syntase

Raf <- c("jg37708.t1", "jg58284.t1", "jg73690.t1")

#####Stachyose synthase###

Sta <- c("jg16432.t1")


test <- pheatmap(log_d.norm4[c(GolS, Raf, Sta), ], cluster_rows = FALSE, cluster_cols = FALSE,
                 fontsize = 22, angle_col = c("45"), color =viridis(12, 0.9), fontsize_row=15, fontsize_col=10,  cellheight = 20,
                 cellwidth = 20, gaps_col = c(seq(0,40, by=4)), gaps_row = seq(1, 6,  by=1)) 
test

pdf("Galactinol_figure/Galactinol_Raf.pdf", width = 29, height = 18)
test
dev.off()


#####Plot the expression of the genes involved in the sucrose metabolism####
#Keep only the expression in F1

colnames(log_d.norm4)
log_d.norm4_F1 <- log_d.norm4[, c(2,6,10,14,18,22,26,30,34,38)]
colnames(log_d.norm4_F1)
colnames(log_d.norm4_F1) <- c("S10_A1_To","S10_A1_Tp","S10_F1_To","S10_F1_Tp","Ti_F1_To",
                              "Ti_F1_Tp","880_F1_To","880_F1_Tp","TO_F1","TP_F1")

DE_Genes_F1 <- DE_Genes[str_which(DE_Genes$Group1, "F1"),]

#Beta-amylase

beta_amylase <- c("jg33586.t1", "jg37743.t1", "jg41600.t1", "jg52568.t1", "jg52571.t1")


#UDP-glucose 4-epimerase

UDP_glucose_epimerase <- c("jg9824.t1", "jg11991.t1","jg11286.t1", "jg29679.t1", "jg32526.t1", "jg42033.t1", "jg71684.t1")


#Alpha-galactosidase

alpha_galactosidase <- c("jg5642.t1", "jg32292.t1", "jg31107.t1", "jg30764.t1", "jg30852.t1")

#Sucrose-phosphate synthase

Sucrose_phosphate_synthase <- c("jg6737.t1", "jg28877.t1", "jg50351.t1")

#Alkaline-invertase

alkaline_invertase <- c("jg7039.t1","jg7040.t1","jg6030.t1","jg38547.t1","jg50395.t1","jg50397.t1","jg58089.t1","jg74325.t1")

#Glucan phosphorylase

glucan_phosphorylase <- c("jg15758.t1", "jg45889.t1", "jg36733.t1", "jg70827.t1","jg73474.t1")


#UDP-Glucose pyrophosphorylase

UDP_glucose_pyrophosphorylase <- c("jg23560.t1", "jg39651.t1")


####Maltose phosphorylase, amylomaltase, 4-alpha-glucanotransferase DPE2  

maltose_phosphorylase <- c("jg28340.t1")

#Sucrose phosphatase

sucrose_phosphatase <- c("jg2759.t1","jg3651.t1", "jg17542.t1")


#####Hexokinase

hexokinase <- c("jg241.t1", "jg53055.t1", "jg22511.t1", "jg29652.t1", "jg38456.t1", "jg30099.t1", "jg70786.t1")


#Glucose phosphase isomerase

glucose_phosphase_isomerase <- c("jg9980.t1", "jg9983.t1", "jg29583.t1", "jg52180.t1", "jg71064.t1")

#Fructokinase

fructokinase <- "jg44614.t1"


##Galactinol synthase 

galactinol <- c( "jg35812.t1", "jg44037.t1", "jg70060.t1")


#Raffinose synthase 

raffinose <- c("jg37708.t1", "jg58284.t1","jg73690.t1")

#Stachyose synthase 

stachyose <- c("jg16432.t1")


#####Heatmap all the genes sucrose metabolism###

genes <- c(beta_amylase, UDP_glucose_epimerase, alpha_galactosidase, Sucrose_phosphate_synthase, alkaline_invertase, glucan_phosphorylase, UDP_glucose_pyrophosphorylase, maltose_phosphorylase, sucrose_phosphatase,
           hexokinase, glucose_phosphase_isomerase, fructokinase, galactinol, raffinose, stachyose)


test <- pheatmap(log_d.norm4_F1[genes, ], cluster_rows = FALSE, cluster_cols = FALSE,
                 fontsize = 22, angle_col = c("45"), color =viridis(12, 0.9), fontsize_row=15, fontsize_col=10,  cellheight = 20,
                 cellwidth = 20, gaps_col = c(2), gaps_row = c(5, 12, 17, 20, 28, 33, 35, 36, 39, 46, 51, 52, 55, 58))


pdf("Sucrose/Heatmap_sucrose.pdf", width = 21, height = 29)
test
dev.off()



