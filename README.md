# LIME

**LIME** is a tool that performs generalized linear mixed model based differential expression analysis on single-cell data with multiple samples. Unlike existing DE methods, we suggest to use raw UMI counts without pre-normalization on data and any normalization technique in the model. We provided two models, Poisson-glmm and Binomial-glmm, modeling UMI counts and zero proportion respectively.

In this vignette, we will demonstrate how to install and use this package.

## Installation

**LIME** can be installed from github directly as follows:
```{r}
devtools::install_github("C-HW/LIME")
library(LIME)
```

## Read Data

We provided a subset of a `SingleCellExperiment` object from Kang et al. 2018 in the package. The original dataset comprised 10X droplet-based scRNA-seq PBCM data from eight Lupus patients obtained before and after 6h-treatment with IFN-β. The subset `Bcells_sce` contains only B cells with 1488 cells in control group and 1392 cells in stimulated group. 

The original dataset is available via `ExperimentHub`.

```{r}
# Download data from ExperimentHub
library(ExperimentHub)
eh = ExperimentHub()
query(eh, "Kang")
sce = eh[["EH2259"]]

# remove undetected genes
sce <- sce[rowSums(counts(sce) > 0) > 0, ]

# remove lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)
```
After removing undetected and lowly expressed genes, the dataset consisted of 29065 cells and 7661 genes. 

```{r}
table(sce$cell, sce$stim)
```

For this vignette, we focus on a cell type, B cells. These cells were contributed from eight different individuals.
```{r}
celltype = "B cells"
Bcells_sce = sce[, sce$cell %in% celltype]
table(Bcells_sce$ind)
```

`Bcells_sce` is also accessible with `data(Bcells_sce)`.
```{r}
data(Bcells_sce)
```

## DE analysis

For this dataset, the goal of the DE analysis is to determine the DEGs between two different conditions 
(groups) within B cells. The function `poisson_glmm_DE` and `binomial_glmm_DE` are designed to perform 
DE analysis. Here we present the workflow of **LIME** by `poisson_glmm_DE` first.

The input requires a `SingleCellExperiment` object `sce` with UMI counts retrievable in 
`sce@assays@data$counts`. The parameter `comparison` indicates the group-of-interest, and 
`replicates` represents the donor of each cell. If there are additional fixed effect and experimental batch 
(such as sample ID), `other_fixed` and `exp_batch` can be used to specify the variables accordingly. 

For `Bcells_sce`, the compared groups and donors is stored in `stim` and `ind` respectively. The experimental 
batch is the combination of `stim` and `ind`. Note that if the donor ID has already represented the experimental 
batch, there's no need to add extra experimental batch variable. And for this example we don't adjust for other 
fixed effects.

```{r}
# create sampleID
Bcells_sce$sampleID = paste0(Bcells_sce$stim,"_", Bcells_sce$ind)
# running time is about 8 minutes.
pois_glmm_df = poisson_glmm_DE(sce = Bcells_sce, comparison = "stim", replicates = "ind", exp_batch = "sampleID")
```

## Determine DEGs

To determine DEGs, we implemented `identifyDEGs` to perform our new criteria, which is a modification of the convention one.  
The conventional criteria selects genes satisfying its adjusted p-value passes 0.05 (default) and its absolute value of log2 
fold change passes log2(1.5) (default). The new criteria are based on the convention plus the gene mean and the difference in mean.

If the log2 gene mean in two groups is lower than a certain value (-2.25 by default) and the log2 mean difference is smaller than a 
threshold (-1 by default), the gene would not be considered as a DEG. Both thresholds are adjustable, depending on the dataset's performance 
and characteristics. More details can be found [here](https://c-hw.github.io/DEanalysis/new_criteria.html).

This function's default setting is applying new criteria. You may set `newcriteria = F` for the conventional rule.
```{r}
# new criteria
pois_glmm_df$new_DEGs = identifyDEGs(adj_pval = pois_glmm_df$BH, 
                                     log2FC = pois_glmm_df$log2FC, 
                                     log2mean = pois_glmm_df$log2mean, 
                                     log2meandiff = pois_glmm_df$log2meandiff, 
                                     log2FCcutoff = 1,
                                     log2meancutoff = -2.25,
                                     log2meandiffcutoff = -1,
                                     newcriteria = T)
```

The volcano plot and heatmap below demonstrate the DEGs.
```{r}
library(ggplot2)
pvalcutoff = 0.05
log2FCcutoff = 1    
ggplot(na.omit(pois_glmm_df), aes(x = log2FC, y = -log10(BH), colour = new_DEGs)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_hline(yintercept=-log10(pvalcutoff),linetype=2) +
    geom_vline(xintercept=log2FCcutoff,linetype=2) + 
    geom_vline(xintercept=-log2FCcutoff,linetype=2) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    scale_color_manual(values = c("gray", "blue")) +
    xlim(c(-3,3)) + 
    ylim(c(0,50)) + 
    xlab("Log2 Fold Change") + 
    ylab("-Log10 (adj.pval)")
```

```{r}
library(pheatmap)

sort_log2FC = sort(abs(pois_glmm_df$log2FC[pois_glmm_df$new_DEGs]), index.return = T, decreasing = T)
sorted_DEGs = pois_glmm_df$genes[which(pois_glmm_df$new_DEGs)][sort_log2FC$ix]
mat = counts(Bcells_sce)[sorted_DEGs,]

#row annotation
annotation_log2FC = data.frame(abs_log2FC = sort_log2FC$x)
rownames(annotation_log2FC) = rownames(mat)

#col annotation
annotation_df = data.frame(Donor = paste0("D",Bcells_sce$ind), Group = Bcells_sce$stim)
rownames(annotation_df) = colnames(mat)
annotation_df = annotation_df[with(annotation_df, order(Group, Donor)), ]

pheatmap(mat[,rownames(annotation_df)], 
         main = paste0("Poisson-glmm DEGs\nin UMI counts (", nrow(mat)," DEGs)"),
         annotation_col = annotation_df, 
         annotation_row = annotation_log2FC,
         cluster_rows=F, cluster_cols=F, 
         show_colnames = F, show_rownames = F,
         annotation_names_row = F,
         color=colorRampPalette(c("navy", "white", "red"))(10), breaks = seq(0,3, length.out = 10)
         )
```

## Alternatives

Some studies have shown that the zero proportion of a gene is an indicator for distinguishing different cell types. Hence, we also implemented Binomial-glmm modeling the zero proportion. The result is similar to Poisson-glmm when the DEGs are chosen from two different cell types. See more examples in the case study 1 [here](https://c-hw.github.io/DEanalysis/index.html#Case_study_1).

The DE result of Binomial-glmm can be obtained by the following command.
```{r}
binomial_glmm_df = binomial_glmm_DE(sce = Bcells_sce, comparison = "stim", replicates = "ind", exp_batch = "sampleID")
```

The function `identifyDEGs` also provides another option for conventional criteria.
```{r}
# conventional criteria
pois_glmm_df$conv_DEGs = identifyDEGs(adj_pval = pois_glmm_df$BH, 
                                     log2FC = pois_glmm_df$log2FC, 
                                     log2FCcutoff = 1,
                                     newcriteria = F)
```

However, the heatmap below shows that the convention one may includes plenty of lowly expressed genes. We still recommend to perform the new criteria.
```{r}
# conventional criteria

sort_log2FC = sort(abs(pois_glmm_df$log2FC[pois_glmm_df$conv_DEGs]), index.return = T, decreasing = T)
sorted_DEGs = pois_glmm_df$genes[which(pois_glmm_df$conv_DEGs)][sort_log2FC$ix]
mat = counts(Bcells_sce)[sorted_DEGs,]

#row annotation
annotation_log2FC = data.frame(abs_log2FC = sort_log2FC$x)
rownames(annotation_log2FC) = rownames(mat)

#col annotation
annotation_df = data.frame(Donor = paste0("D",Bcells_sce$ind), Group = Bcells_sce$stim)
rownames(annotation_df) = colnames(mat)
annotation_df = annotation_df[with(annotation_df, order(Group, Donor)), ]

pheatmap(mat[,rownames(annotation_df)], 
         main = paste0("Poisson-glmm DEGs\nin UMI counts (", nrow(mat)," DEGs)"),
         annotation_col = annotation_df, 
         annotation_row = annotation_log2FC,
         cluster_rows=F, cluster_cols=F, 
         show_colnames = F, show_rownames = F,
         annotation_names_row = F,
         color=colorRampPalette(c("navy", "white", "red"))(10), breaks = seq(0,3, length.out = 10)
         )
```

## Author
Chih-Hsuan Wu (U Chicago)    
Bug report, comments or questions please send to chihhsuan@uchicago.edu.
