MSD_mulitplex_immunoassay_analysis
================
Arun Ghosh
2024-09-04

``` r
library(readxl)
library(dplyr)
library(vsn)
library(imputeLCMD)
library(tidyverse)
library(rstatix)
library(openxlsx)


cytokine_df = data.frame(read.csv("MSD_cytokines_data.csv"))
rownames(cytokine_df) <- cytokine_df[ ,1 ]
cytokine_df <- cytokine_df[,-1 ]

##############################################################
# Printing percentage of missing values (NAs) in original data 

NAs = (sum(is.na(cytokine_df))/prod(dim(cytokine_df)))*100
NAs = round(NAs, 2)
print(paste( "Percentage of missing values (NAs) in data frame : ", NAs))
```

``` r
#----------------------------------------#

normalized_cytokine_matrix = justvsn(as.matrix(cytokine_df), 
                                     minDataPointsPerStratum = 17)
```

``` r
# plot of row standard deviations versus row means of original data
msd <- meanSdPlot(normalized_cytokine_matrix)
```

    ## Warning: Computation failed in `stat_binhex()`.
    ## Caused by error in `compute_group()`:
    ## ! The package "hexbin" is required for `stat_bin_hex()`.

![](README_figs/README4--to%20verify%20the%20fit-1.png)<!-- -->

``` r
preimputed_df = data.frame(normalized_cytokine_matrix)
imputed_QRILC_object = impute.QRILC(preimputed_df, tune.sigma = 0.1)
imputed_df = data.frame(imputed_QRILC_object[1])
```

``` r
###########################################################################
# Preparing table in long-format for statistical analysis

imputed_df5 <- imputed_df

imputed_df6 = data.frame(imputed_df5) %>%
  rownames_to_column(var = "Cytokine") %>%
  pivot_longer(cols = 2:85, names_to = "Sample_No", values_to = "Norm_Cytokine_Conc")

head(imputed_df6)
```

``` r
##########################################################################
#  Friedman's test for overall comparison among exposure groups
#  for each cytokine
##########################################################################

friedman_df = data.frame()

unique_Cytokine <- unique(imputed_df6$Cytokine)

length(unique_Cytokine)


for(i in 1:length(unique_Cytokine)){
  
  filtered_df = imputed_df6 [(imputed_df6$Cytokine == unique_Cytokine[i]),]
  filtered_df = filtered_df %>% 
    separate(Sample_No, into = c("Condensate", "Burn_Condition", "Condensate_Conc", 
                                 "Donor_NvS", "NvS_No", "Donor_FvM", "FvM_No"))
  filtered_df <- filtered_df %>% 
    unite(Exp_Gr, "Condensate", "Burn_Condition", sep = "_", remove = TRUE, na.rm = FALSE)
  filtered_df$Exp_Gr <- as.factor(filtered_df$Exp_Gr)
  filtered_df <- filtered_df %>% 
    unite(Sample_ID, "Donor_NvS", "NvS_No", "Donor_FvM", "FvM_No", sep = "_", remove = TRUE, na.rm = FALSE)
  filtered_df$Sample_ID <- as.factor(filtered_df$Sample_ID)
  filtered_df = filtered_df [ ,c(5,4,2)]
  FT <- friedman_test(filtered_df, Norm_Cytokine_Conc ~ Exp_Gr | Sample_ID)
  friedman_vector = cbind(unique_Cytokine[i], FT$method, FT$n, FT$statistic, FT$p)
  friedman_df = rbind(friedman_df, friedman_vector)
}

colnames(friedman_df) = c("Cytokine", "method", "sample no (n)", "Statistic", "p-Value")

friedman_df$Statistic <- round(as.numeric(friedman_df$Statistic), 4)
friedman_df$'p-Value' <- round(as.numeric(friedman_df$'p-Value'), 9)

friedman_df <- friedman_df[order(friedman_df$'p-Value'),]

head(friedman_df)
```

``` r
######################################################################
# Dunn's Test for Multiple Comparisons
# for paired comparison between exposure groups for each cytokine
######################################################################

dunn_df = data.frame()

unique_Cytokine <- unique(imputed_df6$Cytokine)
length(unique_Cytokine)


for(i in 1:length(unique_Cytokine)){
  
  filtered_df = imputed_df6 [(imputed_df6$Cytokine == unique_Cytokine[i]),]
  filtered_df = filtered_df %>% 
    separate(Sample_No, into = c("Condensate", "Burn_Condition", "Condensate_Conc", 
                                 "Donor_NvS", "NvS_No", "Donor_FvM", "FvM_No"))
  filtered_df <- filtered_df %>% 
    unite(Exp_Gr, "Condensate", "Burn_Condition", sep = "_", remove = TRUE, na.rm = FALSE)
  filtered_df$Exp_Gr <- as.factor(filtered_df$Exp_Gr)
  filtered_df <- filtered_df %>% 
    unite(Sample_ID, "Donor_NvS", "NvS_No", "Donor_FvM", "FvM_No", sep = "_", remove = TRUE, na.rm = FALSE)
  filtered_df$Sample_ID <- as.factor(filtered_df$Sample_ID)
  filtered_df = filtered_df [ ,c(5,4,2)]
  DT <- dunn_test(filtered_df, Norm_Cytokine_Conc ~ Exp_Gr , p.adjust.method = "holm")
  DT <- as.data.frame(DT)
  filtered_df = filtered_df [ ,c(3,1)]
  scale_num = abs(min(filtered_df$Norm_Cytokine_Conc))
  filtered_df <- filtered_df %>%  group_by(Exp_Gr) %>%  summarise_all(mean)
  DT$group1mean = filtered_df$Norm_Cytokine_Conc[match(DT$group1, filtered_df$Exp_Gr)]
  DT$group2mean = filtered_df$Norm_Cytokine_Conc[match(DT$group2, filtered_df$Exp_Gr)]
  DT <- DT %>% rowwise() %>% mutate(logFC = c(log2((group1mean + scale_num)/(group2mean+ scale_num))))
  DT_vector = cbind(unique_Cytokine[i], DT$group1, DT$group1mean, DT$n1, 
                    DT$group2, DT$group2mean, DT$n2, DT$logFC, 
                    DT$statistic, DT$p, DT$p.adj, DT$p.adj.signif)
  dunn_df = rbind(dunn_df, DT_vector)
}

colnames(dunn_df) = c("Cytokine", "Group1",  "Group1-Mean", "Group1 (n)", 
                      "Group2", "Group2-Mean", "Group2 (n)", "logFC", 
                      "Statistic", "pValue", "adj-p Value", "adj-p Value significance")


dunn_df$Group1 <- str_replace_all(dunn_df$Group1, c(CB_F= "Cardboard Flaming",
                                                    CB_S= "Cardboard Smoldering",
                                                    PL_F= "Plastic Flaming",
                                                    PL_S= "Plastic Smoldering",
                                                    PW_F= "Plywood Flaming",
                                                    PW_S= "Plywood Smoldering",
                                                    CL_PBS= "Control"))

dunn_df$Group2 <- str_replace_all(dunn_df$Group2, c(CB_F= "Cardboard Flaming",
                                                    CB_S= "Cardboard Smoldering",
                                                    PL_F= "Plastic Flaming",
                                                    PL_S= "Plastic Smoldering",
                                                    PW_F= "Plywood Flaming",
                                                    PW_S= "Plywood Smoldering",
                                                    CL_PBS= "Control"))

dunn_df <- dunn_df[order(as.numeric(dunn_df$pValue)),]

head(dunn_df)
```

``` r
###############################################################################

MSD <- list()

MSD$Friedman_test <- friedman_df 
MSD$Dunn_test <-  dunn_df
  

blank_excel <- createWorkbook()

Map(function(df, tab_name){     
  
  addWorksheet(blank_excel, tab_name)
  writeData(blank_excel, tab_name, df)
}, 

MSD, names(MSD)
)

saveWorkbook(blank_excel, file = "Table S8.xlsx", overwrite = TRUE)
```
