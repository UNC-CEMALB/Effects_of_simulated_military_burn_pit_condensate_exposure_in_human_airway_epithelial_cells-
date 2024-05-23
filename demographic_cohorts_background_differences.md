demographic_cohorts_background_differences
================
Arun Ghosh
2024-05-21

BACKGROUND DIFFERENCES based on DEMOGRPAHIC COHORTS

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(AnnotationDbi)
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:limma':
    ## 
    ##     plotMA

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

``` r
library(org.Hs.eg.db)
```

    ## 

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
```

    ## Loading required package: ggrepel

``` r
library(eulerr)
library(purrr)
```

    ## 
    ## Attaching package: 'purrr'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     reduce

``` r
library(ggVennDiagram)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ readr     2.1.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%() masks IRanges::%within%()
    ## ✖ dplyr::collapse()     masks IRanges::collapse()
    ## ✖ dplyr::combine()      masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ dplyr::desc()         masks IRanges::desc()
    ## ✖ tidyr::expand()       masks S4Vectors::expand()
    ## ✖ dplyr::filter()       masks stats::filter()
    ## ✖ dplyr::first()        masks S4Vectors::first()
    ## ✖ dplyr::lag()          masks stats::lag()
    ## ✖ ggplot2::Position()   masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()       masks IRanges::reduce()
    ## ✖ dplyr::rename()       masks S4Vectors::rename()
    ## ✖ lubridate::second()   masks S4Vectors::second()
    ## ✖ lubridate::second<-() masks S4Vectors::second<-()
    ## ✖ dplyr::select()       masks AnnotationDbi::select()
    ## ✖ dplyr::slice()        masks IRanges::slice()
    ## ✖ tidyr::unite()        masks ggVennDiagram::unite()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(pheatmap)
library(RColorBrewer) 
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)
```

    ## 
    ## clusterProfiler v4.10.0  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
    ## 
    ## If you use clusterProfiler in published research, please cite:
    ## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:IRanges':
    ## 
    ##     slice
    ## 
    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     rename
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(readxl)
library(psych)
```

    ## 
    ## Attaching package: 'psych'
    ## 
    ## The following objects are masked from 'package:ggplot2':
    ## 
    ##     %+%, alpha
    ## 
    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     distance, reflect

``` r
library(corrplot)
```

    ## corrplot 0.92 loaded

``` r
library(enrichplot)
library(openxlsx)

counts <- read.delim("counts_7G.txt", row.names = 1) # reading in RNAseq data

###############################################################################
##  Analysis without donor specific correction  ###############################

d0 <- DGEList(counts)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 10 # genes expressed in at least 10 samples to be included
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
```

    ## [1] 11374    84

``` r
snames <- colnames(counts) # Sample names

# for identification of non-smoker and smoker donors
nvs <- substr(snames, 1, (nchar(snames)-6))
nvs <- as.factor(nvs)
nvs
```

    ##  [1] N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N N
    ## [39] N N N N S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
    ## [77] S S S S S S S S
    ## Levels: N S

``` r
# for identification of female and male donors
fvm <- substr(snames, 2, (nchar(snames)-5))
fvm <- as.factor(fvm)
fvm
```

    ##  [1] F F F F F F F F F F F F F F F F F F F F F M M M M M M M M M M M M M M M M M
    ## [39] M M M M F F F F F F F F F F F F F F F F F F F F F M M M M M M M M M M M M M
    ## [77] M M M M M M M M
    ## Levels: F M

``` r
# for exposure groups
exposure <- substr(snames, nchar(snames)-3, (nchar(snames)-1))
exposure
```

    ##  [1] "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS"
    ## [13] "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS"
    ## [25] "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS"
    ## [37] "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF"
    ## [49] "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF"
    ## [61] "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF"
    ## [73] "PWS" "PWF" "PLS" "PLF" "CTR" "CBS" "CBF" "PWS" "PWF" "PLS" "PLF" "CTR"

``` r
# for comparison between non-smoker and smoker donors
group2 <- interaction(nvs, exposure)
group2
```

    ##  [1] N.CBS N.CBF N.PWS N.PWF N.PLS N.PLF N.CTR N.CBS N.CBF N.PWS N.PWF N.PLS
    ## [13] N.PLF N.CTR N.CBS N.CBF N.PWS N.PWF N.PLS N.PLF N.CTR N.CBS N.CBF N.PWS
    ## [25] N.PWF N.PLS N.PLF N.CTR N.CBS N.CBF N.PWS N.PWF N.PLS N.PLF N.CTR N.CBS
    ## [37] N.CBF N.PWS N.PWF N.PLS N.PLF N.CTR S.CBS S.CBF S.PWS S.PWF S.PLS S.PLF
    ## [49] S.CTR S.CBS S.CBF S.PWS S.PWF S.PLS S.PLF S.CTR S.CBS S.CBF S.PWS S.PWF
    ## [61] S.PLS S.PLF S.CTR S.CBS S.CBF S.PWS S.PWF S.PLS S.PLF S.CTR S.CBS S.CBF
    ## [73] S.PWS S.PWF S.PLS S.PLF S.CTR S.CBS S.CBF S.PWS S.PWF S.PLS S.PLF S.CTR
    ## 14 Levels: N.CBF S.CBF N.CBS S.CBS N.CTR S.CTR N.PLF S.PLF N.PLS ... S.PWS

``` r
# for comparison between female and male donors
group3 <- interaction(fvm, exposure)
group3
```

    ##  [1] F.CBS F.CBF F.PWS F.PWF F.PLS F.PLF F.CTR F.CBS F.CBF F.PWS F.PWF F.PLS
    ## [13] F.PLF F.CTR F.CBS F.CBF F.PWS F.PWF F.PLS F.PLF F.CTR M.CBS M.CBF M.PWS
    ## [25] M.PWF M.PLS M.PLF M.CTR M.CBS M.CBF M.PWS M.PWF M.PLS M.PLF M.CTR M.CBS
    ## [37] M.CBF M.PWS M.PWF M.PLS M.PLF M.CTR F.CBS F.CBF F.PWS F.PWF F.PLS F.PLF
    ## [49] F.CTR F.CBS F.CBF F.PWS F.PWF F.PLS F.PLF F.CTR F.CBS F.CBF F.PWS F.PWF
    ## [61] F.PLS F.PLF F.CTR M.CBS M.CBF M.PWS M.PWF M.PLS M.PLF M.CTR M.CBS M.CBF
    ## [73] M.PWS M.PWF M.PLS M.PLF M.CTR M.CBS M.CBF M.PWS M.PWF M.PLS M.PLF M.CTR
    ## 14 Levels: F.CBF M.CBF F.CBS M.CBS F.CTR M.CTR F.PLF M.PLF F.PLS ... M.PWS

``` r
#----------------------------------------------------#

mm9 <- model.matrix(~0+group2) # non-smoker and smoker donors
y9 <- voom(d, mm9, plot = T)
```

![](README_figs/README-unnamed-chunk-2-1.png)<!-- -->

FIGURE E5

``` r
# plotting based on smoking status
df1 <- as.data.frame(y9$E)
df1 <- df1 %>% select(ends_with("CTRL"))
nvsLabels <- substr(colnames(df1), 1, (nchar(colnames(df1))-6))
nvsLabels <- str_replace_all(nvsLabels, c(N= "Non-smoker",
                                    S= "Smoker"))
mds <- plotMDS(df1,  pch = 19, cex = 1.2, labels = nvsLabels, gene.selection = "pairwise",
               col=c(rep("darkgreen",6), rep("blue",6)))
```

![](README_figs/README--%20MDS%20plot%20of%20non-smoker%20and%20smoker%20donors-1.png)<!-- -->

``` r
fit9 <- lmFit(y9, mm9) # Fitting linear models in limma

head(coef(fit9))
```

    ##                 group2N.CBF group2S.CBF group2N.CBS group2S.CBS group2N.CTR
    ## ENSG00000225630    4.006262    3.699495    3.673070    3.582915    3.915629
    ## ENSG00000237973    4.155299    4.270518    3.851067    4.158068    4.221899
    ## ENSG00000248527    7.753553    7.722042    7.400547    7.603693    7.858496
    ## ENSG00000228794    3.187638    3.197616    3.203418    3.158135    3.222162
    ## ENSG00000188976    6.351796    6.339733    6.238388    6.317171    6.328424
    ## ENSG00000187961    2.719308    2.952042    2.716415    2.882653    2.742451
    ##                 group2S.CTR group2N.PLF group2S.PLF group2N.PLS group2S.PLS
    ## ENSG00000225630    3.613046    3.953880    3.667804    3.956133    3.580303
    ## ENSG00000237973    4.365459    4.278049    4.338109    4.191064    4.224371
    ## ENSG00000248527    7.784199    7.904980    7.848808    7.819325    7.754334
    ## ENSG00000228794    3.322231    3.154299    3.262144    3.160380    3.168738
    ## ENSG00000188976    6.302884    6.354978    6.291076    6.264854    6.263579
    ## ENSG00000187961    2.793754    2.778630    2.760405    2.720480    3.029797
    ##                 group2N.PWF group2S.PWF group2N.PWS group2S.PWS
    ## ENSG00000225630    3.827800    3.689006    3.935592    3.708112
    ## ENSG00000237973    4.111816    4.235257    4.207138    4.296902
    ## ENSG00000248527    7.681668    7.786258    7.708394    7.762249
    ## ENSG00000228794    3.209780    3.177514    3.090701    3.133401
    ## ENSG00000188976    6.422703    6.453286    6.280424    6.293926
    ## ENSG00000187961    2.737077    3.059765    2.721751    2.930705

``` r
x <- colnames(coef(fit9))
length(x)
```

    ## [1] 14

``` r
x # to see the groups
```

    ##  [1] "group2N.CBF" "group2S.CBF" "group2N.CBS" "group2S.CBS" "group2N.CTR"
    ##  [6] "group2S.CTR" "group2N.PLF" "group2S.PLF" "group2N.PLS" "group2S.PLS"
    ## [11] "group2N.PWF" "group2S.PWF" "group2N.PWS" "group2S.PWS"

``` r
b9 <- list() # for storing analyzed data for significantly altered genes


  contr <- makeContrasts(group2S.CTR-group2N.CTR, levels = colnames(coef(fit9)))
  tmp <- contrasts.fit(fit9, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = row.names(top.table), 
                                 keytype = "ENSEMBL", column = "SYMBOL", 
                                 multiVals="first")) #adding gene names 
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
  top.table <- subset(top.table, top.table$symbol != 'NA')
  top.table<- top.table %>% 
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val < 0.1 & 
                                  abs(top.table$logFC) > 0.5)),]
  top.table <- top.table [(c(7,1,2,3,4,5,6,8))]
  b9 <- top.table
  
###############################################################

counts_CTRL <- counts %>% select(ends_with("CTRL")  ) # selecting control samples

#  Non-smokers vs Smokers
  
  counts_CTRL_NvS <- subset(counts_CTRL, rownames(counts_CTRL) %in% rownames(b9))
  
  counts_CTRL_NvS <- as.data.frame(counts_CTRL_NvS, header = TRUE)
  counts_CTRL_NvS$symbol <- mapIds(org.Hs.eg.db, keys = row.names(counts_CTRL_NvS), 
                                   keytype = "ENSEMBL", column = "SYMBOL", 
                                   multiVals="first") #adding gene names 
```

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
  counts_CTRL_NvS <- subset(counts_CTRL_NvS, counts_CTRL_NvS$symbol != 'NA')
  ncol(counts_CTRL_NvS)
```

    ## [1] 13

``` r
  rownames(counts_CTRL_NvS) <- NULL
  rownames(counts_CTRL_NvS) <- counts_CTRL_NvS$symbol
  counts_CTRL_NvS <- counts_CTRL_NvS[,1:(ncol(counts_CTRL_NvS)-1)]
  
  colnames(counts_CTRL_NvS) <-  sub("[[:digit:]]+","", colnames(counts_CTRL_NvS))
  
  substr(colnames(counts_CTRL_NvS), 2, 3) <- '-'
  colnames(counts_CTRL_NvS) <- gsub('-', '', colnames(counts_CTRL_NvS))
  
  names(counts_CTRL_NvS) <- str_replace_all(names(counts_CTRL_NvS), 
                                            c(NCTRL= "Non-smoker",
                                              SCTRL= "Smoker"))
  temp3 <- as.matrix(counts_CTRL_NvS)
  pheatmap(temp3, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
           display_numbers = FALSE,
           number_color = "black",
           fontsize_number = 5,
           angle_col = c("45"), 
           cellwidth = 10, 
           cellheight = 10, 
           border_color = "black", 
           treeheight_col = 1, 
           fontsize_row = 7, 
           scale = 'row', 
           fontsize_col = 7, 
           cutree_rows = 2, 
           cluster_cols = FALSE, 
           cluster_rows = TRUE)
```

![](README_figs/README--heatmap%20of%20background%20differences%20between%20non-smoker%20and%20smoker%20donors-1.png)<!-- -->

FIGURE 5

``` r
mm10 <- model.matrix(~0+group3) # female and male donors
y10 <- voom(d, mm10, plot = T)
```

![](README_figs/README--MDS%20plot%20of%20female%20and%20male%20donors-1.png)<!-- -->

``` r
# plotting based on sex
df2 <- as.data.frame(y10$E)
df2 <- df2 %>% select(ends_with("CTRL"))
fvmLabels <- substr(colnames(df2), 2, (nchar(colnames(df2))-5))
fvmLabels <- str_replace_all(fvmLabels, c(F= "Female",
                                    M= "Male"))
mds <- plotMDS(df2,  pch = 19, cex = 1.2, labels = fvmLabels, gene.selection = "pairwise",
               col=c(rep("darkgreen",3), rep("blue",3)))
```

![](README_figs/README--MDS%20plot%20of%20female%20and%20male%20donors-2.png)<!-- -->

``` r
fit10 <- lmFit(y10, mm10) # Fitting linear models in limma
  
head(coef(fit10))
```

    ##                 group3F.CBF group3M.CBF group3F.CBS group3M.CBS group3F.CTR
    ## ENSG00000225630    3.832885    3.874558    3.398665    3.845409    3.711615
    ## ENSG00000237973    4.314138    4.111639    3.941017    4.061306    4.404235
    ## ENSG00000248527    7.668478    7.807318    7.344344    7.661091    7.776949
    ## ENSG00000228794    3.204776    3.180181    3.219185    3.145779    3.229154
    ## ENSG00000188976    6.360628    6.330999    6.287337    6.267965    6.298449
    ## ENSG00000187961    2.853720    2.815284    2.810843    2.785328    2.807784
    ##                 group3M.CTR group3F.PLF group3M.PLF group3F.PLS group3M.PLS
    ## ENSG00000225630    3.813962    3.779449    3.846626    3.669809    3.877100
    ## ENSG00000237973    4.182914    4.356591    4.259038    4.277803    4.137622
    ## ENSG00000248527    7.865569    7.811088    7.942766    7.701241    7.872210
    ## ENSG00000228794    3.316203    3.194780    3.219866    3.204732    3.124240
    ## ENSG00000188976    6.332667    6.351065    6.295194    6.291140    6.237291
    ## ENSG00000187961    2.728159    2.774096    2.765375    2.878271    2.861721
    ##                 group3F.PWF group3M.PWF group3F.PWS group3M.PWS
    ## ENSG00000225630    3.813667    3.706802    3.800055    3.848555
    ## ENSG00000237973    4.331295    4.025235    4.457614    4.041321
    ## ENSG00000248527    7.700481    7.767743    7.700833    7.769954
    ## ENSG00000228794    3.234679    3.156555    3.118918    3.103905
    ## ENSG00000188976    6.450750    6.425492    6.287408    6.286971
    ## ENSG00000187961    2.919635    2.882031    2.828911    2.818954

``` r
x <- colnames(coef(fit10))
length(x)
```

    ## [1] 14

``` r
x # to see the groups
```

    ##  [1] "group3F.CBF" "group3M.CBF" "group3F.CBS" "group3M.CBS" "group3F.CTR"
    ##  [6] "group3M.CTR" "group3F.PLF" "group3M.PLF" "group3F.PLS" "group3M.PLS"
    ## [11] "group3F.PWF" "group3M.PWF" "group3F.PWS" "group3M.PWS"

``` r
b10 <- list() # for storing analyzed data fro significantly altered genes
c10 <- list() # for getting the list of all genes

  contr <- makeContrasts(group3F.CTR-group3M.CTR, levels = colnames(coef(fit10)))
  tmp <- contrasts.fit(fit10, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = row.names(top.table), 
                                 keytype = "ENSEMBL", column = "SYMBOL", 
                                 multiVals="first")) #adding gene names 
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
  top.table <- subset(top.table, top.table$symbol != 'NA')
  c10 <- top.table
  
  top.table<- top.table %>% 
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val < 0.1 & 
                                  abs(top.table$logFC) > 0.5)),]
  top.table <- top.table [(c(7,1,2,3,4,5,6,8))]
  b10 <- top.table


  ###############################################################
  #  Female vs Male
  
  counts_CTRL_FvM <- subset(counts_CTRL, rownames(counts_CTRL) %in% rownames(b10))
  
  counts_CTRL_FvM <- as.data.frame(counts_CTRL_FvM, header = TRUE)
  counts_CTRL_FvM$symbol <- mapIds(org.Hs.eg.db, keys = row.names(counts_CTRL_FvM), 
                                   keytype = "ENSEMBL", column = "SYMBOL", 
                                   multiVals="first") #adding gene names 
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
  counts_CTRL_FvM <- subset(counts_CTRL_FvM, counts_CTRL_FvM$symbol != 'NA')
  ncol(counts_CTRL_FvM)
```

    ## [1] 13

``` r
  rownames(counts_CTRL_FvM) <- NULL
  rownames(counts_CTRL_FvM) <- counts_CTRL_FvM$symbol
  counts_CTRL_FvM <- counts_CTRL_FvM[,1:(ncol(counts_CTRL_FvM)-1)]
  
  colnames(counts_CTRL_FvM) <-  sub("[[:digit:]]+","", colnames(counts_CTRL_FvM))
  names(counts_CTRL_FvM) <- substring(names(counts_CTRL_FvM),2,6)
  
  names(counts_CTRL_FvM) <- str_replace_all(names(counts_CTRL_FvM), c(FCTRL= "Female",
                                                                      MCTRL= "Male"))
  temp4 <- as.matrix(counts_CTRL_FvM)
  pheatmap(temp4, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
           display_numbers = FALSE,
           number_color = "black",
           fontsize_number = 5,
           angle_col = c("45"), 
           cellwidth = 10, 
           cellheight = 7, 
           border_color = "black", 
           treeheight_col = 1, 
           fontsize_row = 7, 
           scale = 'row', 
           fontsize_col = 7, 
           cutree_rows =  2, 
           cluster_cols = TRUE, 
           cluster_rows = TRUE)
```

![](README_figs/README--heatmap%20of%20background%20differences%20between%20female%20and%20male%20donors-1.png)<!-- -->

``` r
GdataFvM <- c10 

##########################################################################
#---------------------------------------------#
# comparing between background expression between female and male donors

GdataFvM <- GdataFvM %>% mutate(ProbeID = rownames(GdataFvM))   

GdataFvM <- select(GdataFvM, ProbeID, logFC)
rownames(GdataFvM) <- NULL

#making ranked gene list
genelist_GdataFvM = GdataFvM[,2] #numeric vector
names(genelist_GdataFvM) = as.character(GdataFvM[,1]) #named vector
genelist_GdataFvM = sort(genelist_GdataFvM, decreasing = TRUE) #must sort in descending order

#Performing GSEA analysis
#Gene Ontology (GO) 
gseGO_FvM_ALL <- gseGO(geneList=genelist_GdataFvM, 
                       ont ="ALL", 
                       keyType = "ENSEMBL", 
                       minGSSize = 10, #min size of gene sets for analysis
                       maxGSSize = 500, #max size of gene sets for analysis
                       pvalueCutoff = 0.05, 
                       eps = 0,
                       verbose = TRUE, 
                       OrgDb = organism, 
                       pAdjustMethod = "BH")
```

    ## preparing geneSet collections...

    ## GSEA analysis...

    ## leading edge analysis...

    ## done...

``` r
enrichplot::dotplot(gseGO_FvM_ALL, split=".sign", title = "Female vs Male",
                    showCategory = 10, font.size = 15, label_format = 70) + facet_grid(.~.sign)
```

![](README_figs/README--GSEA%20analysis-1.png)<!-- -->

``` r
enrichplot::dotplot(gseGO_FvM_ALL, split=".sign", title = "Female vs Male",
                    showCategory = 10, font.size = 15, label_format = 70)+ facet_grid(.~.sign)+ 
  scale_x_continuous(breaks = seq(0, 1, by = 1), limits=c(0,1))+
  theme(axis.text.y = element_text(lineheight = 0.80, size = 15),
        title = element_text(size = 15, face="bold"),
        plot.title = element_text(hjust=0.5))
```

![](README_figs/README--GSEA%20analysis-2.png)<!-- -->

FIGURE E6

``` r
#GO over-representation analysis in female donors
geneFvM <- names(genelist_GdataFvM[genelist_GdataFvM[] > 0.5]) # logFC > 0.5
enGOFvM <- enrichGO(gene         = geneFvM,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)
head(enGOFvM)
```

    ##            ONTOLOGY         ID                                     Description
    ## GO:0014829       BP GO:0014829   vascular associated smooth muscle contraction
    ## GO:0001666       BP GO:0001666                             response to hypoxia
    ## GO:0052548       BP GO:0052548            regulation of endopeptidase activity
    ## GO:0036293       BP GO:0036293             response to decreased oxygen levels
    ## GO:0031223       BP GO:0031223                               auditory behavior
    ## GO:1903053       BP GO:1903053 regulation of extracellular matrix organization
    ##            GeneRatio   BgRatio       pvalue  p.adjust     qvalue
    ## GO:0014829     4/137  28/21261 2.995787e-05 0.0325928 0.02907813
    ## GO:0001666    10/137 328/21261 5.414181e-05 0.0325928 0.02907813
    ## GO:0052548    10/137 345/21261 8.239295e-05 0.0325928 0.02907813
    ## GO:0036293    10/137 347/21261 8.641937e-05 0.0325928 0.02907813
    ## GO:0031223     3/137  15/21261 1.125269e-04 0.0325928 0.02907813
    ## GO:1903053     5/137  75/21261 1.240918e-04 0.0325928 0.02907813
    ##                                                                                                                                                                     geneID
    ## GO:0014829                                                                                                 ENSG00000274286/ENSG00000151617/ENSG00000127129/ENSG00000070961
    ## GO:0001666 ENSG00000197635/ENSG00000151617/ENSG00000167772/ENSG00000159167/ENSG00000132170/ENSG00000156113/ENSG00000161544/ENSG00000148926/ENSG00000129521/ENSG00000168140
    ## GO:0052548 ENSG00000132170/ENSG00000206072/ENSG00000196104/ENSG00000169242/ENSG00000106366/ENSG00000022556/ENSG00000158125/ENSG00000152377/ENSG00000129521/ENSG00000215301
    ## GO:0036293 ENSG00000197635/ENSG00000151617/ENSG00000167772/ENSG00000159167/ENSG00000132170/ENSG00000156113/ENSG00000161544/ENSG00000148926/ENSG00000129521/ENSG00000168140
    ## GO:0031223                                                                                                                 ENSG00000079215/ENSG00000174469/ENSG00000184564
    ## GO:1903053                                                                                 ENSG00000197635/ENSG00000225614/ENSG00000125398/ENSG00000170961/ENSG00000168646
    ##            Count
    ## GO:0014829     4
    ## GO:0001666    10
    ## GO:0052548    10
    ## GO:0036293    10
    ## GO:0031223     3
    ## GO:1903053     5

``` r
#barplot
barplot(enGOFvM, showCategory=10, font = 15, title = "Female vs Male")+
  scale_x_continuous(breaks = seq(0, 40, by = 10), limits=c(0,40))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))
```

![](README_figs/README-unnamed-chunk-3-1.png)<!-- -->

``` r
#------------------------------------------------------------------#
#GO over-representation analysis in male donors
geneMvF <- names(genelist_GdataFvM[genelist_GdataFvM[] < -0.5]) # logFC < -0.5
enGOMvF <- enrichGO(gene         = geneMvF,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)
head(enGOMvF)
```

    ##            ONTOLOGY         ID                          Description GeneRatio
    ## GO:0044782       BP GO:0044782                  cilium organization    37/241
    ## GO:0098813       BP GO:0098813       nuclear chromosome segregation    34/241
    ## GO:0000070       BP GO:0000070 mitotic sister chromatid segregation    27/241
    ## GO:0000819       BP GO:0000819         sister chromatid segregation    29/241
    ## GO:0060271       BP GO:0060271                      cilium assembly    35/241
    ## GO:0007059       BP GO:0007059               chromosome segregation    38/241
    ##              BgRatio       pvalue     p.adjust       qvalue
    ## GO:0044782 425/21261 3.767311e-22 9.697058e-19 8.716367e-19
    ## GO:0098813 358/21261 1.295220e-21 1.666948e-18 1.498365e-18
    ## GO:0000070 201/21261 2.722377e-21 1.728638e-18 1.553816e-18
    ## GO:0000819 244/21261 2.823656e-21 1.728638e-18 1.553816e-18
    ## GO:0060271 396/21261 3.357882e-21 1.728638e-18 1.553816e-18
    ## GO:0007059 493/21261 7.106953e-21 3.048883e-18 2.740541e-18
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     geneID
    ## GO:0044782                 ENSG00000170264/ENSG00000048342/ENSG00000092850/ENSG00000162643/ENSG00000138587/ENSG00000165115/ENSG00000106477/ENSG00000198003/ENSG00000087903/ENSG00000256061/ENSG00000158113/ENSG00000188316/ENSG00000105792/ENSG00000155530/ENSG00000197748/ENSG00000123977/ENSG00000152763/ENSG00000153789/ENSG00000163075/ENSG00000113966/ENSG00000183833/ENSG00000186710/ENSG00000165698/ENSG00000204815/ENSG00000127399/ENSG00000120051/ENSG00000103021/ENSG00000156206/ENSG00000173588/ENSG00000163060/ENSG00000146038/ENSG00000164675/ENSG00000087586/ENSG00000142731/ENSG00000169126/ENSG00000103494/ENSG00000111834
    ## GO:0098813                                                                 ENSG00000071539/ENSG00000166851/ENSG00000114346/ENSG00000122952/ENSG00000100918/ENSG00000087586/ENSG00000066279/ENSG00000198901/ENSG00000013810/ENSG00000164109/ENSG00000076382/ENSG00000138160/ENSG00000237649/ENSG00000117724/ENSG00000134057/ENSG00000117650/ENSG00000112742/ENSG00000123219/ENSG00000142945/ENSG00000137804/ENSG00000089685/ENSG00000170312/ENSG00000169679/ENSG00000088325/ENSG00000121152/ENSG00000131747/ENSG00000157456/ENSG00000138778/ENSG00000175063/ENSG00000156970/ENSG00000090889/ENSG00000126787/ENSG00000101057/ENSG00000117399
    ## GO:0000070                                                                                                                                                                                 ENSG00000071539/ENSG00000166851/ENSG00000122952/ENSG00000198901/ENSG00000164109/ENSG00000076382/ENSG00000138160/ENSG00000237649/ENSG00000117724/ENSG00000134057/ENSG00000117650/ENSG00000112742/ENSG00000123219/ENSG00000142945/ENSG00000137804/ENSG00000089685/ENSG00000170312/ENSG00000169679/ENSG00000088325/ENSG00000121152/ENSG00000138778/ENSG00000175063/ENSG00000156970/ENSG00000090889/ENSG00000126787/ENSG00000101057/ENSG00000117399
    ## GO:0000819                                                                                                                                                 ENSG00000071539/ENSG00000166851/ENSG00000122952/ENSG00000198901/ENSG00000013810/ENSG00000164109/ENSG00000076382/ENSG00000138160/ENSG00000237649/ENSG00000117724/ENSG00000134057/ENSG00000117650/ENSG00000112742/ENSG00000123219/ENSG00000142945/ENSG00000137804/ENSG00000089685/ENSG00000170312/ENSG00000169679/ENSG00000088325/ENSG00000121152/ENSG00000131747/ENSG00000138778/ENSG00000175063/ENSG00000156970/ENSG00000090889/ENSG00000126787/ENSG00000101057/ENSG00000117399
    ## GO:0060271                                                 ENSG00000170264/ENSG00000048342/ENSG00000092850/ENSG00000162643/ENSG00000138587/ENSG00000165115/ENSG00000106477/ENSG00000198003/ENSG00000087903/ENSG00000256061/ENSG00000158113/ENSG00000105792/ENSG00000155530/ENSG00000197748/ENSG00000123977/ENSG00000152763/ENSG00000153789/ENSG00000163075/ENSG00000113966/ENSG00000183833/ENSG00000186710/ENSG00000165698/ENSG00000204815/ENSG00000127399/ENSG00000120051/ENSG00000103021/ENSG00000156206/ENSG00000173588/ENSG00000163060/ENSG00000146038/ENSG00000164675/ENSG00000142731/ENSG00000169126/ENSG00000103494/ENSG00000111834
    ## GO:0007059 ENSG00000169220/ENSG00000071539/ENSG00000166851/ENSG00000012048/ENSG00000114346/ENSG00000122952/ENSG00000100918/ENSG00000087586/ENSG00000066279/ENSG00000198901/ENSG00000013810/ENSG00000164109/ENSG00000076382/ENSG00000138160/ENSG00000237649/ENSG00000117724/ENSG00000148773/ENSG00000134057/ENSG00000117650/ENSG00000112742/ENSG00000123219/ENSG00000163535/ENSG00000142945/ENSG00000137804/ENSG00000089685/ENSG00000170312/ENSG00000169679/ENSG00000088325/ENSG00000121152/ENSG00000131747/ENSG00000157456/ENSG00000138778/ENSG00000175063/ENSG00000156970/ENSG00000090889/ENSG00000126787/ENSG00000101057/ENSG00000117399
    ##            Count
    ## GO:0044782    37
    ## GO:0098813    34
    ## GO:0000070    27
    ## GO:0000819    29
    ## GO:0060271    35
    ## GO:0007059    38

``` r
#barplot
barplot(enGOMvF, showCategory=10, font = 15, title = "Male vs Female")+
  scale_x_continuous(breaks = seq(0, 40, by = 10), limits=c(0,40))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))
```

![](README_figs/README-unnamed-chunk-3-2.png)<!-- -->

``` r
#-----------------------------------------------------------------#  
  
write.table(b9, "background_NvS.csv", sep = ",", row.names = F)
write.table(b10, "background_FvM.csv", sep = ",", row.names = F)

#-----------------------------------------------------------------# 
```
