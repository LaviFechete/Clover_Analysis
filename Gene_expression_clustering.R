
#This script uses the normalized gene counts obtained from "Subgenome_Normalization.R"

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



######Identify the genes that have a very low or very high subgenome ratio######

sd_cutoff = 1.5

colnames(ratios_split3)[1:20]

for (i in colnames(ratios_split3)[1:20]){
  #i="Ti_F3"
  hist(ratios_split3[,i])
  mean.log.ratio <- mean(ratios_split3[,i])
  mean.log.ratio
  sd.log.ratio <- sd(ratios_split3[,i], na.rm=TRUE)
  sd.log.ratio
  ratio.cutoff <- mean.log.ratio + sd_cutoff*sd.log.ratio
  ratio.cutoff
  abline(v=ratio.cutoff, col="red")
  abline(v=-ratio.cutoff, col="red")
  
  ratios_split3.subset <- subset(ratios_split3, ratios_split3[,i] > ratio.cutoff)
  ratios_split3.subset.filter <- ratios_split3[,i] > ratio.cutoff
  ratios_split3[,paste0(i, "outliers_High")] <- ratios_split3.subset.filter
  head(subset(ratios_split3, paste0(i, "outliers") ==TRUE))
  
  ratios_split3.subset <- subset(ratios_split3, ratios_split3[,i] < -(ratio.cutoff))
  ratios_split3.subset.filter <- ratios_split3[,i] < -(ratio.cutoff)
  ratios_split3[,paste0(i, "outliers_Low")] <- ratios_split3.subset.filter
  head(subset(ratios_split3, paste0(i, "outliers") ==TRUE))
  
}

colnames(ratios_split3)


####From the differentially expressed genes select the outlier genes for
#F1, F2 and F3 condition#######

ratios_split3_DE <- ratios_split3[unique(DE_Gene_F1$gene_ID),]

#High ratio
ratios_split3_DE$S10_outliers_High <- "FALSE"
ratios_split3_DE[ratios_split3_DE$S10_F2outliers_High=="TRUE"&
                   ratios_split3_DE$S10_F3outliers_High=="TRUE"& ratios_split3_DE$S10_F4outliers_High=="TRUE",]$S10_outliers_High <- "TRUE"

ratios_split3_DE$Ti_outliers_High <- "FALSE"
ratios_split3_DE[ratios_split3_DE$Ti_F2outliers_High=="TRUE"&
                   ratios_split3_DE$Ti_F3outliers_High=="TRUE"& ratios_split3_DE$Ti_F4outliers_High=="TRUE",]$Ti_outliers_High <- "TRUE"

ratios_split3_DE$X88_outliers_high <- "FALSE"
ratios_split3_DE[ratios_split3_DE$X88_F2outliers_High=="TRUE"&
                   ratios_split3_DE$X88_F3outliers_High=="TRUE"& ratios_split3_DE$X88_F4outliers_High=="TRUE",]$X88_outliers_high <- "TRUE"

table(ratios_split3_DE$S10_outliers_High)
table(ratios_split3_DE$Ti_outliers_High)
table(ratios_split3_DE$X88_outliers_high)

####Low  ratio###

ratios_split3_DE$S10_outliers_Low <- "FALSE"
ratios_split3_DE[ratios_split3_DE$S10_F2outliers_Low=="TRUE"&
                   ratios_split3_DE$S10_F3outliers_Low=="TRUE"& ratios_split3_DE$S10_F4outliers_Low=="TRUE",]$S10_outliers_Low <- "TRUE"

ratios_split3_DE$Ti_outliers_Low <- "FALSE"
ratios_split3_DE[ratios_split3_DE$Ti_F2outliers_Low=="TRUE"&
                   ratios_split3_DE$Ti_F3outliers_Low=="TRUE"& ratios_split3_DE$Ti_F4outliers_Low=="TRUE",]$Ti_outliers_Low <- "TRUE"

ratios_split3_DE$X88_outliers_Low <- "FALSE"
ratios_split3_DE[ratios_split3_DE$X88_F2outliers_Low=="TRUE"&
                   ratios_split3_DE$X88_F3outliers_Low=="TRUE"& ratios_split3_DE$X88_F4outliers_Low=="TRUE",]$X88_outliers_Low <- "TRUE"

table(ratios_split3_DE$S10_outliers_Low)
table(ratios_split3_DE$Ti_outliers_Low)
table(ratios_split3_DE$X88_outliers_Low)



#Generate the hierarchial clustering and plotting in Figure 3 b and c 

c("#a6cee3","#E6B8BF", "#1f78b4", "#54990F", "#A3CC7A")

######
ratio_outliers_X880_low <- filter(ratios_split3_DE, X880_outliers_Low == "TRUE")

test1 <- unique((ratio_outliers_X880_low[,c(9:12, 17:20)]))
pheatmap(cor(test1, method = "spearman"))

gene_dist <- dist(test1)
gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("X880_", "8800 ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "8800"),]$genotype <- "8800"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"

freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("8800 F0", "8800 F1", "8800 F2", "8800 F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"),
                          labels = c("F0 8800", "F1 8800", "F2 8800", "F3 8800", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p1 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#1f78b4", "#a6cee3"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("8800 Low ratios")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+ 
  geom_text(data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -0.5,
  size=5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 3))
p1


###########

ratio_outliers_Ti_low <- filter(ratios_split3_DE, Ti_outliers_Low == "TRUE")

test1 <- unique((ratio_outliers_Ti_low[,c(9:16)]))

deviance_Ti_out_low <- deviance_test_DE[rownames(ratio_outliers_Ti_low), "Ti_Tp_prog"]
mean(deviance_Ti_out_low)

plotDensities(deviance_Ti_out_low)

gene_dist <- dist(test1)
gene_hclust <- stats::hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("Ti_", "Ti ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "Ti"),]$genotype <- "Ti"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"


freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("Ti F0", "Ti F1", "Ti F2", "Ti F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"),
                          labels = c("F0 Ti", "F1 Ti", "F2 TI", "F3 Ti", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p2 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#a6cee3","#54990F"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("Ti Low ratios")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+ 
  geom_text(data    = dat_text,
          mapping = aes(x = -Inf, y = -Inf, label = label),
          hjust   = -0.1,
          vjust   = -0.5,
          size=5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 2))
p2

#####################################

ratio_outliers_S10_low <- filter(ratios_split3_DE, S10_outliers_Low == "TRUE")

test1 <- unique((ratio_outliers_S10_low[,c(5:12)]))

gene_dist <- dist(test1)
gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("S10_", "S10 ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "S10"),]$genotype <- "S10"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"

freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("S10 F0", "S10 F1", "S10 F2", "S10 F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"),
                          labels = c("F0 S10", "F1 S10", "F2 S10", "F3 S10", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p3 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#a6cee3", "#A3CC7A"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("S10 Low ratios")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+ 
  geom_text(data    = dat_text,
          mapping = aes(x = -Inf, y = -Inf, label = label),
          hjust   = -0.1,
          vjust   = -0.5,
          size=5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 3))
p3

############

ratio_outliers_X880_high <- filter(ratios_split3_DE, X880_outliers_high == "TRUE")

test1 <- unique((ratio_outliers_X880_high[,c(9:12, 17:20)]))

gene_dist <- dist(test1)
gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("X880_", "8800 ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "8800"),]$genotype <- "8800"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"

freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("8800 F0", "8800 F1", "8800 F2", "8800 F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"),
                          labels = c("F0 8800", "F1 8800", "F2 8800", "F3 8800", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p4 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#1f78b4", "#a6cee3"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("8800 High ratios")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+
  geom_text(data    = dat_text,
          mapping = aes(x = -Inf, y = -Inf, label = label),
          hjust   = -0.1,
          vjust   = -0.5,
          size=5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 3))
p4

###########

ratio_outliers_Ti_high <- filter(ratios_split3_DE, Ti_outliers_High == "TRUE")

test1 <- unique((ratio_outliers_Ti_high[,c(9:16)]))

gene_dist <- dist(test1)
gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("Ti_", "Ti ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "Ti"),]$genotype <- "Ti"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"


freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("Ti F0", "Ti F1", "Ti F2", "Ti F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"),
                          labels = c("F0 Ti", "F1 Ti", "F2 TI", "F3 Ti", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p5 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#a6cee3","#54990F"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("Ti High ratios")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+
  geom_text( data    = dat_text,
          mapping = aes(x = -Inf, y = -Inf, label = label),
          hjust   = -0.1,
          vjust   = -0.5,
          size=5 )+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 2))
p5


#####################################

ratio_outliers_S10_high <- filter(ratios_split3_DE, S10_outliers_High == "TRUE")

test1 <- unique((ratio_outliers_S10_high[,c(5:12)]))

gene_dist <- dist(test1)
gene_hclust <- hclust(gene_dist, method = "complete")

plot(gene_hclust, labels = FALSE)
abline(h = 15, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

gene_cluster <- cutree(gene_hclust, k = 3) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  dplyr::rename(cluster = value, gene = name)

test1$gene <- rownames(test1)

clustering <- test1 %>% 
  inner_join(gene_cluster, by = "gene")

clustering <- pivot_longer(clustering, cols = c(1:(ncol(clustering)-2)))
#clustering$group <- str_sub(clustering$name, 1, -2)
clustering$name <- gsub("S10_", "S10 ", clustering$name)
clustering$name <- gsub("TP_TO_", "Prog ", clustering$name)

clustering$genotype <- "empty"
clustering[str_which(clustering$name, "S10"),]$genotype <- "S10"
clustering[str_which(clustering$name, "Prog"),]$genotype <- "Prog"
clustering$genotype <- factor(clustering$genotype, levels = c("S10", "Prog"))

freq <- as.data.frame(table(gene_cluster$cluster))
freq$new_clusters <- paste0(freq$Freq, " transcripts")

clustering$name <- factor(clustering$name, levels = c("S10 F0", "S10 F1", "S10 F2", "S10 F3", "Prog F0", "Prog F1", "Prog F2", "Prog F3"), 
                          labels = c("F0 S10", "F1 S10", "F2 S10", "F3 S10", "F0 Prog", "F1 Prog", "F2 Prog", "F3 Prog"))

dat_text <- data.frame(
  cluster   = c("1", "2", "3"),
  label = c(freq$new_clusters[1],freq$new_clusters[2], freq$new_clusters[3]),
  cyl   = c("1", "2", "3")
)

p6 <- clustering %>% 
  ggplot(aes(name, value)) +
  geom_line(aes(group = gene, color=genotype)) +
  scale_colour_manual(values = c("#a6cee3","#A3CC7A"))+
  geom_line(stat = "summary", fun = "median", colour = c("#042652"), size = 1.5, 
            aes(group = 1)) +
  facet_wrap(~cluster) +
  theme_minimal()+
  xlab("")+
  ylab(bquote("log"['e']*"(Tr"['To']*"/Tr"['Tp']*")"))+
  ggtitle("S10 High Ratios")+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "white"))+
  geom_text(data    = dat_text,
          mapping = aes(x = -Inf, y = -Inf, label = label),
          hjust   = -0.1,
          vjust   = -0.5,
          size=5)+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 3))
p6


p3/p2/p1+ plot_layout(guides = "collect")
p6/p5/p4+ plot_layout(guides = "collect")


png("Figures/Clustering_low_Ratios_F1-F3_New.png", res = 400, width = 3600, height =3600)
p3/p2/p1+ plot_layout(guides = "collect")
dev.off()
# 
png("Figures/Clustering_high_RatiosF1-F3_New.png", res = 400, width = 3600, height =3600)
p6/p5/p4+ plot_layout(guides = "collect")
dev.off()


