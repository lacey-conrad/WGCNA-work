library(edgeR)
library(WGCNA)
library(magrittr)
library(reshape)
library(Hmisc)
library(lme4)
library(lmerTest)
library(multcomp)
library(gprofiler2)
library(igraph)
library(plyr)
library(vegan)
library(DESeq2)
library(tidyverse) 
#library(anRichment)
library(gprofiler2)



options(stringsAsFactors = FALSE);
allowWGCNAThreads()

## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
## Dataframe prep - ALL DATA               ##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
thermo_counts <- read.table("thermo_capacity_counts.txt", header = TRUE, sep = "\t")
datExpr <- thermo_counts[, -c(2:6)] 

id <- read.csv("thermo_capacity_rnaseq_samples.csv")
id <- id[order(id$mouse_id),] #order sample ID's alphabetically

rownames(datExpr) <- datExpr$Geneid #replace rownames with GeneIDs
datExpr$Geneid <- NULL #rownames MUST be GeneIDs. Columns must ONLY be individuals
colnames(datExpr) <- id$ID #WARNING. This replaces the individuals IDs, so make sure they are in order
# ------------------------------------------ #
##############################################


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
## GASTROC ONLY ANALYSIS                   ##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
gastrocID <- subset(id, tissue=="gas") #subset ids by tissue
gastrocID <- gastrocID[order(gastrocID$mouse_id),]
gasID <- gastrocID$ID #make a vector containing gastroc mouse ids
gastroc <- datExpr[,colnames(datExpr) %in% gasID] #subset ontogeny data by gastroc mouse ids
colnames(gastroc) <- gastrocID$mouse_id 
#gastroc[1:5,1:10]
#dim(gastroc)
write.csv(gastroc, "gastroc_raw_counts.csv")
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
gastroc$mean <- rowMeans(gastroc) #rowMeans takes means of each row
keep_gastroc <- subset(gastroc, mean >= 10) #filter by means
keep_gastroc$mean <- NULL #clean up dataset
# ------------------------------------------ #
##############################################


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####
## NEW- normalization using DESeq          ####
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ####


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
de_input = as.matrix(keep_gastroc[,-1])
meta_df <- data.frame(Sample = names(keep_gastroc[-1])) %>%
  mutate(Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .))
# DESeqDataSet Object (DESeq2 object) created from DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(round(de_input), meta_df, design = ~Type)
# DESeq2 function
dds <- DESeq(dds)

# make a transformed count matrix, using variance stabilizing transformation (VST)
vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
#rv_wpn <- rowVars(wpn_vsd)
#summary(rv_wpn)
expr_normalized <- wpn_vsd
# ------------------------------------------ #
##############################################


## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
## WGCNA for Gastroc
## Much of the following code was adapted from:
## https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
input_mat = t(expr_normalized)
input_mat[1:5,1:10]
write.csv(as.data.frame(input_mat),file="norm_expr_data.csv")
# ------------------------------------------ #
##############################################

#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
# cluster the samples to check for outliers
sampleTree = hclust(dist(input_mat), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", 
     xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
picked_power = 8
temp_cor <- cor   
# Force it to use WGCNA cor function (fix a namespace conflict issue)
cor <- WGCNA::cor         
netwk <- blockwiseModules(input_mat,
    # == Adjacency Function ==
    power = picked_power, networkType = "signed",
    # == Tree and Block Options ==
    deepSplit = 2, pamRespectsDendro = F, minModuleSize = 30, maxBlockSize = 4000,
    # == Module Adjustments ==
    reassignThreshold = 0, mergeCutHeight = 0.25,
    # == TOM == Archive the run results in TOM file (saves time)
    saveTOMs = T, saveTOMFileBase = "ER",
    # == Output Options
    numericLabels = T, verbose = 3)
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netwk$dendrograms[[1]], mergedColors[netwk$blockGenes[[1]]],
  "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dim(table(mergedColors))
table(mergedColors)
write.csv(as.data.frame(table(mergedColors)),file="merged_colors_new.csv")
# ------------------------------------------ #
##############################################


#mergedColors
#table(netwk$colors)


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
module_df <- data.frame(gene_id = names(netwk$colors), colors = labels2colors(netwk$colors))
write_delim(module_df, file = "gene_modules.txt", delim = "\t")
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
# Add treatment names
MEs0$treatment = row.names(MEs0)
# tidy & plot data
pdf('corr_plot.pdf', h=4, w=7)
mME = MEs0 %>% pivot_longer(-treatment) %>% mutate(name = gsub("ME", "", name),
    name = factor(name, levels = module_order))
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) + geom_tile() +
  theme_bw() + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
    midpoint = 0, limit = c(-1,1)) + theme(axis.text.x = element_text(angle=90)) + 
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

dev.off()
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
# pick out a few modules of interest here
modules_of_interest = c("tan", "turquoise","blue")
# Pull out list of genes in that module
submod = module_df %>% subset(colors %in% modules_of_interest)
row.names(module_df) = module_df$gene_id
# Get normalized expression for those genes
expr_normalized[1:5,1:10]
subexpr = expr_normalized[submod$gene_id,]
submod_df = data.frame(subexpr) %>% mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id) %>% mutate(module = module_df[gene_id,]$colors)
submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment", y = "normalized expression")
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
genes_of_interest <- module_df %>% subset(colors %in% modules_of_interest)
genes_of_interest
head(genes_of_interest)
genes_go <- genes_of_interest$gene_id
genes_go
dim(genes_of_interest)
expr_of_interest <- expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]
write.csv(data.frame(expr_of_interest),file="gene_expr.csv")
write.csv(data.frame(table(genes_of_interest)),file="genes_of_int.csv")
table(expr_of_interest)
# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest), power = picked_power)
# Add gene names to row and columns
row.names(TOM) <- row.names(expr_of_interest)
colnames(TOM) <- row.names(expr_of_interest)
edge_list <- data.frame(TOM) %>% mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1) %>% dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>% subset(!(gene1==gene2)) %>%
  mutate(module1 = module_df[gene1,]$colors, module2 = module_df[gene2,]$colors)
head(edge_list)
# ------------------------------------------ #
##############################################


#############################################
# +++++++++++++++++++++++++++++++++++++++++ #
# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list, file = "edgelist.tsv", delim = "\t")
# ------------------------------------------ #
##############################################



## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
## Gene functional enrichment              ##
## Peromyscus maniculatus bairdii          ##
## @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ##
# try M. musculus ***
# includes pathways

gostres <- gost(query = genes_go, 
                organism = "pmbairdii", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

gos_result <- gostres$result


head(gos_result)
df_ge <- data.frame(gos_result)
#head(df_ge)
#details(gos_result)
df5 <- apply(df_ge,2,as.character)
write.csv(df5, file="func_enrich.csv", row.names = FALSE)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p
