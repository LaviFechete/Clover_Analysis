
####Network analysis for white clover cold expression using the WGCNA package###

####Creating gene expression modules###

library(WGCNA)

#Load the differentially expressed genes

DE_Genes <- readxl::read_xlsx("DE_GLM_Up_no_amb_separated_Ratios_Annotations1.xlsx")
length(unique(DE_Genes$gene_ID))

#Data normalized separately for the two subgenomes

module_logcpm <- as.data.frame(log2(d.norm3_1+1))

plot(colSums(module_logcpm))

module_logcpm2 <- module_logcpm[,colnames(module_logcpm) %in% DE_Genes$gene_ID]

#Calculate the soft power
sft <- pickSoftThreshold(module_logcpm2,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed")

#Calculate the R squared
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)


ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(sft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

#Chose a softpower that has an R2 of at least 0.8 to avoid noisy results

net <- blockwiseModules(module_logcpm2, power = 18,
                        corType = "bicor", # use robust correlation
                        networkType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.15,
                        numericLabels = T, pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        saveTOMFileBase = "TOM",
                        verbose = 3)

#Find the number of genes in each module
table(net$colors)

# Convert labels to colors for plotting
mergedColors = net$colors
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

pheatmap::pheatmap(MEs, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(MEs)


net$colors[names(net$colors) %in% "jg35812.t1"]

TP_test <- net$colors[net$colors %in% "6"]

test <- module_logcpm[, colnames(module_logcpm) %in% names(TP_test)]


#correlate with the traits

metadata <- readxl::read_xlsx("Metadata_Network.xlsx")

table(metadata$Sample==rownames(module_logcpm2))

traits <- model.matrix(~0+Condition, data = metadata)
traits <- model.matrix(~0+Genotype, data = metadata)
traits <- model.matrix(~0+Sample1, data = metadata)
traits <- model.matrix(~0+Subgenome, data = metadata)


# Define numbers of genes and samples
nGenes = ncol(module_logcpm2);
nSamples = nrow(module_logcpm2);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(module_logcpm2, moduleLabels)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p");
moduleTraitCor <- moduleTraitCor[order(rownames(moduleTraitCor)),]
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pheatmap(moduleTraitCor, cluster_rows = F, cluster_cols = F, display_numbers = textMatrix, main = "Module-trait relationships",
         color = blueWhiteRed(50), angle_col = 45, fontsize = 20, cellwidth = 100, cellheight = 30, fontsize_number = 12)

moduleTraitCor <- as.data.frame(moduleTraitCor)
moduleTraitCor$order <- as.numeric(gsub("ME", "", rownames(moduleTraitCor)))
moduleTraitCor <- moduleTraitCor[order(as.numeric(moduleTraitCor$order)),]
moduleTraitCor <- as.matrix(moduleTraitCor[, c(1:2)])

#plot the correalation values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)


p<- pheatmap(moduleTraitCor, cluster_rows = F, cluster_cols = F, display_numbers = textMatrix, main = "Module-trait relationships",
             color = blueWhiteRed(50), fontsize_number = 18, fontsize_row = 20, fontsize_col = 20,  cellwidth = 100, cellheight = 50, angle_col = 45, fontsize = 20)
p

pdf(file = "Module_correlation_Subgenome.pdf", width = 20, height = 29)
p
dev.off()

####



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
  
  pdf(paste0(name, ".pdf"), width = 7, height =7)
  treemapPlot(reducedTerms)
  dev.off()
  return(treemapPlot(reducedTerms))
}

for(i in unique(net$colors)){
  DE_Genes_Module <- DE_Genes %>% filter(gene_ID %in% names(net$colors[net$colors %in% i]))
  results <- Gene_enrichment(DE_Genes_Module$gene_ID, paste0("Modules_DE/Module_", i))
  Gene_enrichment_fig(results, paste0("Modules_DE/Module_", i))
}