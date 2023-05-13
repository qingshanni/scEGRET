## EGRET: a simple method for trajectory evaluation 

An R package for evaluating the continuity of inferred single-cell trajectories using multiple metrics, including differential gene expression testing, average expression correlation, and gene expression enthalpy. Associated visualizations for pairwise bin comparisons are included. May be useful for selecting trajectories for downstream validation in complex biological systems, as well as comparisons of trajectories inferred by different algorithms. 

<p align="center">
  <img height="400" src="images/schematic.png" />
</p>


## Installation

To install from GitHub:
```
library(devtools)
install_github("qingshanni/scEGRET", dependencies = T)
```

## An introduction to scEGRET package

As a demonstration, we have included here a full walkthrough in which we start from a gene expression matrix, calculate a pseudotime trajectory, and then assess the trajectory using scEGRET. In this exercise, we will show how our assessment can interact with single-cell data stored as a Seurat object, and trajectories calculated using monocle2. However, alternative trajectories calculated using other packages(e.g. slingshot), and other forms of single-cell data objects (e.g. SingleCellExperiment)  can also be used, as long as the relevant data slots can be accessed and loaded.

Load libraries
```
library(Seurat)
library(monocle)
library(scEGRET)
```

We will analyze a simulated dataset, stored inside the package, containing 1,000 cells and 10,000 genes. This dataset has been simulated using splatter as three distinct groups falling along a single, simple linear trajectory.

Create a Seurat object 
```
sim.seur <- CreateSeuratObject(counts = sim.expr,
                               meta.data = sim.meta)
```

We can first visualize what this dataset looks like using UMAP by processing it through Seurat with default parameters. 

```
sim.seur <- NormalizeData(sim.seur)
sim.seur <- FindVariableFeatures(sim.seur, nfeatures = 3000)
sim.seur <- ScaleData(sim.seur)
sim.seur <- RunPCA(sim.seur)
sim.seur <- RunUMAP(sim.seur, reduction = 'pca', dims = 1:30)
```

We can see that the three groups fall into three continuous, but non-overlapping, sections of a single arc-like shape in UMAP.

```
DimPlot(sim.seur, group.by = 'Group')
```
<p align="center">
  <img height="400" src="images/01.png" />
</p>

To continue our analysis, we then infer a pseudotime trajectory for this dataset using monocle2. For convenience, we have included two wrapper functions for processing.

```
sim.mono <- seurat2Monocle(sim.seur)
sim.mono <- RunMonocle2(sim.mono, 
                varGenes = sim.seur@assays$RNA@var.features)
plot(sim.mono$Pseudotime, sim.mono$TrueStep)
```
<p align="center">
  <img height="400" src="images/02.png" />
</p>

We can also visualize how this looks in monocle2 itself.

```
plot_cell_trajectory(sim.mono, color_by = 'Pseudotime')
```

<p align="center">
  <img height="400" src="images/03.png" />
</p>

```
sim.mono$Group <- as.factor(sim.mono$Group)
plot_cell_trajectory(sim.mono, color_by = 'Group')
```
<p align="center">
  <img height="400" src="images/04.png" />
</p>

We can now move onto scEGRET to assess the continuity of this trajectory,We first collect the relevant data from their respective slots.

```
input.pseudotime <- sim.mono$Pseudotime
input.df <- as.data.frame(sim.seur@assays$RNA@counts)
input.meta <- sim.seur@meta.data
```

...and now we can generate a scEGRET with a flap of its wings
```
sim.egret <- egret.flap(df = input.df,
                        metadata = input.meta,
                        pseudotimes = input.pseudotime)
```
Calculations can be run using a single command. we recommend using the rice rule for histogram bin number to determine the number of bins to set of course, one should explore the impact of modifying this parameter on their own dataset prior to making a final determination

```
sim.egret <- egret.fly(sim.egret, steps = 20)
```

Congratulations! you're done.

To visualize the output, we can use a number of different visualizations. For instance, to see how step enthalpy changes

```
Soar_StepDist(sim.egret)
```

<p align="center">
  <img height="400" src="images/05.png" />
</p>

It may be helpful to include a cutoff line for reference Because our graphics are depicted using ggplot2, we can do this through 
```
Soar_StepDist(sim.egret)+geom_hline(yintercept = 0.5)
```


<p align="center">
  <img height="400" src="images/06.png" />
</p>

While we have found 0.5 is a useful reference, more adjustments can be made based that take into consideration the minimum step change found across multiple bins and so on. To help with this, we can visualize all possible enthalpy differences between bins as a pyramid

```
Soar_Mountain_StepChange(sim.egret)
```

<p align="center">
  <img height="400" src="images/07.png" />
</p>


sessionInfo

```
sessionInfo()
```

```
R version 4.0.5 (2021-03-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.3.so

locale:
 [1] LC_CTYPE=zh_CN.UTF-8       LC_NUMERIC=C               LC_TIME=zh_CN.UTF-8       
 [4] LC_COLLATE=zh_CN.UTF-8     LC_MONETARY=zh_CN.UTF-8    LC_MESSAGES=zh_CN.UTF-8   
 [7] LC_PAPER=zh_CN.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=zh_CN.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets 
 [9] methods   base     

other attached packages:
 [1] scEGRET_0.1.0       monocle_2.18.0      DDRTree_0.1.5       irlba_2.3.3        
 [5] VGAM_1.1-5          ggplot2_3.3.5       Biobase_2.50.0      BiocGenerics_0.36.1
 [9] Matrix_1.3-4        SeuratObject_4.0.4  Seurat_4.0.5       

loaded via a namespace (and not attached):
  [1] Rtsne_0.15            colorspace_2.0-2      deldir_1.0-6         
  [4] ellipsis_0.3.2        ggridges_0.5.3        proxy_0.4-26         
  [7] rstudioapi_0.13       spatstat.data_3.0-1   farver_2.1.0         
 [10] leiden_0.3.9          listenv_0.8.0         ggrepel_0.9.1        
 [13] RSpectra_0.16-0       fansi_0.5.0           codetools_0.2-18     
 [16] docopt_0.7.1          knitr_1.36            polyclip_1.10-0      
 [19] ggstream_0.1.0        jsonlite_1.7.2        ica_1.0-2            
 [22] cluster_2.1.1         png_0.1-7             pheatmap_1.0.12      
 [25] uwot_0.1.10           shiny_1.7.1           sctransform_0.3.2    
 [28] spatstat.sparse_3.0-1 compiler_4.0.5        httr_1.4.2           
 [31] assertthat_0.2.1      fastmap_1.1.0         lazyeval_0.2.2       
 [34] limma_3.46.0          later_1.3.0           htmltools_0.5.2      
 [37] tools_4.0.5           igraph_1.2.9          gtable_0.3.0         
 [40] glue_1.5.0            RANN_2.6.1            reshape2_1.4.4       
 [43] dplyr_1.0.7           Rcpp_1.0.7            scattermore_0.7      
 [46] slam_0.1-49           vctrs_0.3.8           nlme_3.1-152         
 [49] lmtest_0.9-39         xfun_0.39             stringr_1.4.0        
 [52] globals_0.14.0        mime_0.12             miniUI_0.1.1.1       
 [55] lifecycle_1.0.1       goftest_1.2-3         future_1.23.0        
 [58] MASS_7.3-53.1         zoo_1.8-9             scales_1.1.1         
 [61] spatstat.core_2.3-2   promises_1.2.0.1      spatstat.utils_3.0-3 
 [64] RColorBrewer_1.1-2    yaml_2.2.1            reticulate_1.22      
 [67] pbapply_1.5-0         gridExtra_2.3         rpart_4.1-15         
 [70] fastICA_1.2-3         stringi_1.7.5         densityClust_0.3     
 [73] rlang_0.4.12          pkgconfig_2.0.3       matrixStats_0.61.0   
 [76] qlcMatrix_0.9.7       evaluate_0.14         lattice_0.20-41      
 [79] ROCR_1.0-11           purrr_0.3.4           tensor_1.5           
 [82] labeling_0.4.2        patchwork_1.1.1       htmlwidgets_1.5.4    
 [85] cowplot_1.1.1         tidyselect_1.1.1      parallelly_1.29.0    
 [88] RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_2.0.1       
 [91] R6_2.5.1              generics_0.1.1        combinat_0.0-8       
 [94] DBI_1.1.1             pillar_1.6.4          withr_2.4.2          
 [97] mgcv_1.8-34           fitdistrplus_1.1-6    survival_3.2-10      
[100] abind_1.4-5           tibble_3.1.6          future.apply_1.8.1   
[103] crayon_1.4.2          KernSmooth_2.23-18    utf8_1.2.2           
[106] spatstat.geom_3.2-1   plotly_4.10.0         rmarkdown_2.11       
[109] viridis_0.6.2         grid_4.0.5            data.table_1.14.2    
[112] FNN_1.1.3             matrixTests_0.2.2     sparsesvd_0.2        
[115] HSMMSingleCell_1.10.0 digest_0.6.28         xtable_1.8-4         
[118] tidyr_1.1.4           httpuv_1.6.3          munsell_0.5.0        
[121] viridisLite_0.4.0 
```

