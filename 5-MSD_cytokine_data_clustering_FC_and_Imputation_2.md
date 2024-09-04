MSD_cytokine_clustering
================
Arun Ghosh
2024-09-04

``` r
library(readxl)
library(WGCNA)
library(plotrix)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vsn)
library(imputeLCMD)


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

cytokine_df = data.frame(read_excel("MSD_cytokines_imputed.xlsx", 
                                    sheet = 2, col_names = TRUE, col_types = NULL, 
                                    na = "", skip = 0))
```

``` r
#######################################################
gsg_2 = goodSamplesGenes(cytokine_df[,1:(ncol(cytokine_df)-1)], verbose = 3);
gsg_2$allOK

sampleTree_2 = hclust(dist(cytokine_df[,1:(ncol(cytokine_df)-1)]), method = "average");
# Plot the sample tree

plot(sampleTree_2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
```

![](README_figs/README5--checking%20for%20data%20quality%20and%20outliers-1.png)<!-- -->

``` r
datExpr_2 = cytokine_df
rownames(datExpr_2) <- datExpr_2[,ncol(datExpr_2)]
datExpr_2 <- datExpr_2[,1:(ncol(datExpr_2)-1)]

nCytokines = ncol(datExpr_2)
nSamples = nrow(datExpr_2)

#upload traits data

traitData_2 = data.frame(read_excel("MSD_cytokines_imputed.xlsx", 
                                    sheet = 3, col_names = TRUE, col_types = NULL, 
                                    na = "", skip = 0))
dim(traitData_2)
names(traitData_2)

# Form a data frame analogous to expression data that will hold the demographic traits.
rownames(datExpr_2)
Samples_2 = rownames(datExpr_2)
traitRows_2 = match(Samples_2, traitData_2$ID)
datTraits_2 = traitData_2[traitRows_2, -5]
rownames(datTraits_2) = traitData_2[, 5]
```

``` r
# Re-cluster samples
sampleTree2_2 = hclust(dist(datExpr_2), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_2 = numbers2colors(datTraits_2, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_2, traitColors_2,
                    groupLabels = names(datTraits_2), 
                    main = "Sample dendrogram and trait heatmap")
```

![](README_figs/README5--Re-cluster%20samples-1.png)<!-- -->

``` r
# Re-cluster samples 2 new
sampleTree2_2 = hclust(dist(datExpr_2), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_2 = numbers2colors(datExpr_2, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_2, traitColors_2,
                    groupLabels = names(datExpr_2), 
                    main = "Sample dendrogram and trait heatmap")
```

![](README_figs/README5--Re-cluster%20samples%202%20new-1.png)<!-- -->

``` r
save(datExpr_2, datTraits_2, file = "01-dataInput.RData")
# Load the data saved in the first part
lnames = load(file = "01-dataInput.RData");
##The variable lnames contains the names of loaded variables.
lnames
# Allow multi-threading within WGCNA.
allowWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr_2, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.2,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

``` r
# tried different combinations
# power 5 and minimoduleSize 3 showed acceptable clustering
#=====================================================================================


net = blockwiseModules(datExpr_2, power = 5,
                       TOMType = "unsigned", minModuleSize = 3,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM", 
                       verbose = 3)
#no modules detected? changed power to 4 and 
#mergeCutHeight to lower and higher values, still no modules detected
table(net$colors)
tableList <- net$colors
```

``` r
color.id("#66ccee")
color.id("#228833")
color.id("#CCBB44")

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors, colorSeq = c("steelblue1","forestgreen","darkkhaki")) # ,"maroon", "grey"

# "grey", ,"maroon"
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
```

    ## Warning in pmin(objHeights[dendro$order][floor(positions)],
    ## objHeights[dendro$order][ceiling(positions)]): an argument will be fractionally
    ## recycled

![](README_figs/README5--colors%20assigned%20to%20clusters-1.png)<!-- -->

``` r
moduleLabels = net$colors
moduleColors =labels2colors(net$colors, zeroIsGrey = TRUE, colorSeq = c("steelblue1","forestgreen","darkkhaki")) # ,"maroon","grey"
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

# Define numbers of genes and samples
ncytokines = ncol(datExpr_2);
nSamples = nrow(datExpr_2);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_2, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits_2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


display.brewer.pal(11, "RdBu")
```

![](README_figs/README5-unnamed-chunk-2-1.png)<!-- -->

``` r
coul <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(50))

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_2),
               yLabels = names(MEs),
               #ySymbols = names(MEs),
               ySymbols = c("steelblue1", "darkkhaki", "forestgreen"),
               colorLabels = FALSE,
               colors = coul,
               textMatrix = textMatrix,
               cex.text = .9,
               setStdMargins = FALSE,
               zlim = c(-1,1),
               legendLabel = "Correlation",
               main = paste("Cluster-Cohort Relationships"))
```

``` r
clst1 <- names(datExpr_2)[moduleColors=="steelblue1"] 
length(clst1)
clst1

clst2 <- names(datExpr_2)[moduleColors=="forestgreen"]
length(clst2)
clst2

clst3 <- names(datExpr_2)[moduleColors=="darkkhaki"]
length(clst3)
clst3
```

``` r
imputed_df = data.frame(read_excel("MSD_cytokines_imputed.xlsx", 
                                    sheet = 1, col_names = TRUE, col_types = NULL, 
                                    na = "", skip = 0))

rownames(imputed_df) <- imputed_df[ ,1 ]
imputed_df <- imputed_df[,-1 ]

normalized_cytokine_df2 = as.data.frame(t(imputed_df))
unique_Cytokine <- colnames(normalized_cytokine_df2)
length(unique_Cytokine)

normalized_cytokine_df2$exposure <- substr(rownames(normalized_cytokine_df2), 
                                        1, (nchar(rownames(normalized_cytokine_df2))-11))

exps <- unique(normalized_cytokine_df2$exposure)
length(exps)

cluster1 <- data.frame(matrix(ncol = length(clst1), nrow = length(exps)))
colnames(cluster1) <- clst1
cluster1$exposure <- exps

cluster2 <- data.frame(matrix(ncol = length(clst2), nrow = length(exps)))
colnames(cluster2) <- clst2
cluster2$exposure <- exps

cluster3 <- data.frame(matrix(ncol = length(clst3), nrow = length(exps)))
colnames(cluster3) <- clst3
cluster3$exposure <- exps

imputedHM <- normalized_cytokine_df2 %>%
  group_by(exposure) %>%
  summarise_all(mean)

imputedHM$exposure <- str_replace_all(imputedHM$exposure, c(CB_F= "Cardboard Flaming",
                                                    CB_S= "Cardboard Smoldering",
                                                    PL_F= "Plastic Flaming",
                                                    PL_S= "Plastic Smoldering",
                                                    PW_F= "Plywood Flaming",
                                                    PW_S= "Plywood Smoldering",
                                                    CL_PBS_P= "Control"))

imputedHM <- data.frame(imputedHM, row.names = imputedHM$exposure)
imputedHM <- imputedHM[, 2:ncol(imputedHM)]
imputedHM <- t(imputedHM)
imputedHM <- imputedHM[,c(3,1,2,6,7,4,5)]

clts <- c(clst1, clst2, clst3) 
imputedHM <- imputedHM[clts,]

pheatmap(imputedHM,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
         display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 7,
         angle_col = c("45"), 
         cellwidth = 20, 
         cellheight = 15,
         border_color = "black", 
         treeheight_col = 6, 
         fontsize_row = 10,
         scale = 'row',
         fontsize_col = 10,
         cluster_cols = FALSE,
         cluster_rows = FALSE) 
```

![](README_figs/README5--reading%20in%20imputed%20data%20and%20generating%20heatmap-1.png)<!-- -->

``` r
for(i in 1:length(unique_Cytokine)){
  filtered_df = normalized_cytokine_df2[,(colnames(normalized_cytokine_df2) == unique_Cytokine[i] |
                                            colnames(normalized_cytokine_df2) == "exposure")]
  if (any(filtered_df[unique_Cytokine[i]]<0)){
  scale_num = abs(min(filtered_df[unique_Cytokine[i]]))
  filtered_df[unique_Cytokine[i]] <- filtered_df[unique_Cytokine[i]] + scale_num}
  
  filtered_df1 <- filtered_df %>%
    group_by(exposure) %>%
    summarise_all(mean)
  ctrl = filtered_df1[grep(pattern = "CL_PBS_P", x = filtered_df1$exposure) , ]
  ctrl <- as.numeric(ctrl[,unique_Cytokine[i]])
  filtered_df1[,unique_Cytokine[i]] <- mapply('/', filtered_df1[,unique_Cytokine[i]], ctrl)
  filtered_df1 <- as.data.frame(filtered_df1)
  rownames(filtered_df1) <- filtered_df1[,1]
  #filtered_df1 <- filtered_df1[,-1]
  if(unique_Cytokine[i] %in% colnames(cluster1))
    {cluster1 <- merge(cluster1, filtered_df1, by.x = "exposure", by.y = "exposure")}
  if(unique_Cytokine[i] %in% colnames(cluster2))
    {cluster2 <- merge(cluster2, filtered_df1, by.x = "exposure", by.y = "exposure")}
  if(unique_Cytokine[i] %in% colnames(cluster3))
    {cluster3 <- merge(cluster3, filtered_df1, by.x = "exposure", by.y = "exposure")}
}

cluster1 <- cluster1 %>% select(where(~!all(is.na(.x))))
rownames(cluster1) <- cluster1[,1]
cluster1 <- cluster1[,-1]
cluster1$c1 <- rowMeans(cluster1)
cluster1$exposure <- rownames(cluster1)
cluster1$c1l <- log2(cluster1$c1)

cluster2 <- cluster2 %>% select(where(~!all(is.na(.x))))
rownames(cluster2) <- cluster2[,1]
cluster2 <- cluster2[,-1]
cluster2$c2 <- rowMeans(cluster2)
cluster2$exposure <- rownames(cluster2)
cluster2$c2l <- log2(cluster2$c2)

cluster3 <- cluster3 %>% select(where(~!all(is.na(.x))))
rownames(cluster3) <- cluster3[,1]
cluster3 <- cluster3[,-1]
cluster3$c3 <- rowMeans(cluster3)
cluster3$exposure <- rownames(cluster3)
cluster3$c3l <- log2(cluster3$c3)
```

``` r
cluster <- data.frame(matrix(ncol = 3, nrow = length(exps)))
cluster$exposure <- exps
cluster <- merge(cluster, cluster1[,c("exposure", "c1l")], by.x = "exposure", by.y = "exposure")
cluster <- merge(cluster, cluster2[,c("exposure", "c2l")], by.x = "exposure", by.y = "exposure")
cluster <- merge(cluster, cluster3[,c("exposure", "c3l")], by.x = "exposure", by.y = "exposure")
cluster <- cluster %>% select(where(~!all(is.na(.x))))
rownames(cluster) <- cluster[,1]
cluster <- cluster[,-1]
cluster <- cluster %>% filter(!grepl('PBS', rownames(cluster)))

rownames(cluster) <- str_replace_all(rownames(cluster), c(CB_F= "Cardboard Flaming",
                                                          CB_S= "Cardboard Smoldering",
                                                          PL_F= "Plastic Flaming",
                                                          PL_S= "Plastic Smoldering",
                                                          PW_F= "Plywood Flaming",
                                                          PW_S= "Plywood Smoldering"))

colnames(cluster) <- str_replace_all(colnames(cluster), c(c1l= "cluster 1",
                                                          c2l= "cluster 2",
                                                          c3l= "cluster 3"))

cluster <- cluster[c(1,2,5,6,3,4),]

cluster = cluster %>%
  rownames_to_column(var = "exposure") %>%
  pivot_longer(cols = 2:4, names_to = "cluster", values_to = "score")

q <- ggplot(cluster,aes(x = exposure, y =score, fill = cluster))+ 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c('steelblue1','forestgreen','darkkhaki')) +
  coord_cartesian(ylim = c(-0.2, 0.4))
q + theme(axis.text=element_text(size=12, colour="black"),
          axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_text(angle=45, hjust=1))
```

![](README_figs/README5--Plotting%20cluster%20scores%20as%20bar%20diagram-1.png)<!-- -->
