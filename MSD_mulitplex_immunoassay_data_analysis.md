MSD_mulitplex_immunoassay_analysis
================
Arun Ghosh
2024-05-22

``` r
suppressWarnings(suppressMessages(library(readxl)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(vsn)))
suppressWarnings(suppressMessages(library(imputeLCMD)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(rstatix)))
suppressWarnings(suppressMessages(library(openxlsx)))


cytokine_df = data.frame(read.csv("MSD_cytokines_data.csv"))
rownames(cytokine_df) <- cytokine_df[ ,1 ]
cytokine_df <- cytokine_df[,-1 ]

##############################################################
# Printing percentage of missing values (NAs) in original data 

NAs = (sum(is.na(cytokine_df))/prod(dim(cytokine_df)))*100
NAs = round(NAs, 2)
print(paste( "Percentage of missing values (NAs) in data frame : ", NAs))
```

    ## [1] "Percentage of missing values (NAs) in data frame :  1.19"

``` r
#----------------------------------------#

normalized_cytokine_matrix = justvsn(as.matrix(cytokine_df), 
                                     minDataPointsPerStratum = 17)
```

``` r
# plot of row standard deviations versus row means of original data
msd <- meanSdPlot(normalized_cytokine_matrix)
```

![](README_figs/README--to%20verify%20the%20fit-1.png)<!-- -->

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

    ## # A tibble: 6 × 3
    ##   Cytokine Sample_No       Norm_Cytokine_Conc
    ##   <chr>    <chr>                        <dbl>
    ## 1 GM.CSF   CB_F_H_NS_1_F_1              2.82 
    ## 2 GM.CSF   CB_F_H_NS_1_F_2              1.22 
    ## 3 GM.CSF   CB_F_H_NS_3_M_1              2.79 
    ## 4 GM.CSF   CB_F_H_NS_4_F_3              0.281
    ## 5 GM.CSF   CB_F_H_NS_5_M_2              2.75 
    ## 6 GM.CSF   CB_F_H_NS_6_M_3              0.614

``` r
##########################################################################
#  Friedman's test for overall comparison among exposure groups
#  for each cytokine
##########################################################################

friedman_df = data.frame()

unique_Cytokine <- unique(imputed_df6$Cytokine)

length(unique_Cytokine)
```

    ## [1] 17

``` r
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

    ##                         Cytokine        method sample no (n) Statistic
    ## Friedman chi-squared14     IP.10 Friedman test            12   38.6429
    ## Friedman chi-squared15     MCP.1 Friedman test            12   27.7500
    ## Friedman chi-squared2       IL.7 Friedman test            12   22.6429
    ## Friedman chi-squared12 Eotaxin.3 Friedman test            12   21.5714
    ## Friedman chi-squared4      IFN.g Friedman test            12   21.4286
    ## Friedman chi-squared8      IL.1b Friedman test            12   18.2857
    ##                            p-Value
    ## Friedman chi-squared14 0.000000841
    ## Friedman chi-squared15 0.000104714
    ## Friedman chi-squared2  0.000925357
    ## Friedman chi-squared12 0.001447508
    ## Friedman chi-squared4  0.001536030
    ## Friedman chi-squared8  0.005556474

``` r
######################################################################
# Dunn's Test for Multiple Comparisons
# for paired comparison between exposure groups for each cytokine
######################################################################

dunn_df = data.frame()

unique_Cytokine <- unique(imputed_df6$Cytokine)
length(unique_Cytokine)
```

    ## [1] 17

``` r
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

    ##      Cytokine               Group1      Group1-Mean Group1 (n)
    ## 308     IP.10              Control  3.8865611663222         12
    ## 262 Eotaxin.3 Cardboard Smoldering 1.41349624017621         12
    ## 311     IP.10      Plastic Flaming 3.11845547697867         12
    ## 296     IP.10    Cardboard Flaming 2.64782035655393         12
    ## 42      IL.1a      Plywood Flaming 1.00083964957413         12
    ## 56       IL.7              Control 1.24877394743327         12
    ##                 Group2        Group2-Mean Group2 (n)             logFC
    ## 308    Plywood Flaming   1.89060721884534         12 0.769853630308083
    ## 262    Plywood Flaming 0.0240968877085607         12 0.603205766520793
    ## 311    Plywood Flaming   1.89060721884534         12 0.519796930510837
    ## 296            Control    3.8865611663222         12 -0.42787309164076
    ## 42  Plywood Smoldering -0.545080374922614         12 0.837592038763722
    ## 56     Plywood Flaming  0.140426054780445         12  1.24625760357641
    ##             Statistic               pValue         adj-p Value
    ## 308 -5.40590506056197 6.44819472815336e-08 1.3541208929122e-06
    ## 262 -3.55651648721182 0.000375804905680061 0.00789190301928128
    ## 311 -3.48120202042381 0.000499168847203548 0.00998337694407096
    ## 296  3.13810278283396  0.00170045237720148  0.0323085951668282
    ## 42  -3.06278831604595  0.00219285119997425  0.0460498751994593
    ## 56  -3.03768349378328  0.00238404209647853  0.0500648840260492
    ##     adj-p Value significance
    ## 308                     ****
    ## 262                       **
    ## 311                       **
    ## 296                        *
    ## 42                         *
    ## 56                        ns

``` r
###############################################################################
#wb = xlsx::createWorkbook()

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
```

    ## $Friedman_test
    ## [1] 0
    ## 
    ## $Dunn_test
    ## [1] 0

``` r
saveWorkbook(blank_excel, file = "Table E8.xlsx", overwrite = TRUE)
```
