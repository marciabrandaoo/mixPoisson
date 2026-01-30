# mixPoisson
Long-term Survival Models Based on a New Mixed Poisson Distribution for Competing Causes
by Márcia Brandão, Diego I. Gallardo and Marcelo Bourguignon 

This repository contains all data and R scripts required to fully reproduce the simulation study and real data application presented in the article:

For questions, comments or remarks about the code, please contact Márcia Brandão (mbrandao@ufam.edu.br).

The codes have been written using the R software (platform: Windows 11, 64-bit, version 23H2) and organized to ensure transparency, clarity, and reproducibility of all numerical results, tables, and figures reported in the paper.
To reproduce the corresponding results of application and simulation studies in the manuscript, please follow the steps as below. 



## Reproducibility – Simulation Study (Section 4)
 
The simulation study described in Section 4 of the article can be reproduced by running the scripts below in the specified order.

Step 1. Set your working directory to the folder where the file auxiliary_functions.R is located.

Step 2. Source the auxiliary functions file:

```r
source("auxiliary_functions.R")

```


Reproducing Tables 1 and 2 – Section 4.1 (Performance of the ML estimators in finite samples)

Step 3 (option A – full simulation). Source the file:
```r
source("simulation/Simulator - Performance of the ML estimators in finite samples.R")

```

Step 3 (option B – precomputed results). Alternatively, load the RData file provided by the authors:
```r
load("simulation/MC-Study-MPIG_MNBcr.RData")

```

Then run:
```r
sim_MC_results1
sim_MC_results2

```

Reproducing Table 3 – Section 4.2 (Misspecification model)

Step 4 (option A – full simulation). Source the file:
```r
source("simulation/Simulator - Misspecification model.R")

```

Step 4 (option B – precomputed results). Alternatively, load the RData file provided by the authors:
```r
load("simulation/MC-Study-Misspecification-MIG-MGA-POI.RData")

```

Then run:
```r
resumo_all

```

## Random Seeds and Reproducibility
All simulation scripts set fixed random seeds at the beginning of execution to ensure full reproducibility of the results. The seed values are explicitly defined within each script.




## Reproducibility – Application (Section 5)
 
The application results presented in Section 5 of the article can be reproduced by running the scripts below in the specified order.

Step 1. Set your working directory to the folder where the file auxiliary functions.R and the dataset data_breast_cancer.csv are located.

Step 2. Source the auxiliary functions file:
```r
source("auxiliary_functions.R")

```

Step 3. Load the dataset:
```r
data_BC <- read.csv("data/data_breast_cancer.csv", sep = ",")

```

Reproducing Table 4 – Section 5 (Application to medical dataset)

Step 4 (option A – full analysis).
Source the file:
```r
source("application/application_breast_cancer_table_selection_criteria.R")

```

Step 4 (option B – precomputed results).

Alternatively, load the RData file provided by the authors:
```r
load("application/Application-Table-Selection_criteria.RData")

```
Then run:
```r
criterios_df

```

Reproducing Tables 5, 6 and 7, and Figures 2, 3, 4 and 5 – Section 5

Step 5 (option A – full analysis).
Source the file:
```r
source("application/Application_breast_cancer-MPIGcr and PTcr.R")

```

Step 6.
Source the file:
```r
source("application/Application_breast_cancer_tables_figures.R")

```
or

Step 5 (option B – precomputed results).
Alternatively, load the RData file provided by the authors:
```r
load("application/Application_breast_cancer_MPIG_POISSON.RData")

```

Step 6.
Source the file:
```r
source("application/Application_breast_cancer_tables_figures.R")

```

All Tables (5, 6 and 7) and Figures (2, 3, 4 and 5) will be saved in the working directory where the files auxiliary_functions.R and data_breast_cancer.csv are located.




## Computational Requirements
R version: R- 4.5.1 (or higher) 
Operating system: Windows 
The analyses were conducted using base R and freely available CRAN packages.

 
## Required R Packages
The following R packages are required to run the codes:
R version 4.5.1 (2025-06-13 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Matrix products: default
  LAPACK version 3.12.1

locale:
[1] LC_COLLATE=Portuguese_Brazil.utf8  LC_CTYPE=Portuguese_Brazil.utf8    LC_MONETARY=Portuguese_Brazil.utf8 LC_NUMERIC=C                      
[5] LC_TIME=Portuguese_Brazil.utf8    

time zone: Etc/GMT+4
tzcode source: internal

attached base packages:
 [1] splines   stats4    compiler  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] VGAM_1.1-13      Bessel_0.6-1     quantreg_6.1     SparseM_1.84-2   nortest_1.0-4    timereg_2.0.7    ggsurvfit_1.2.0  survminer_0.5.1 
 [9] ggpubr_0.6.2     survival_3.8-3   ggfortify_0.4.19 foreign_0.8-90   lubridate_1.9.4  forcats_1.0.0    stringr_1.5.2    dplyr_1.1.4     
[17] purrr_1.1.0      readr_2.1.5      tidyr_1.3.1      tibble_3.3.0     ggplot2_4.0.0    tidyverse_2.0.0  devtools_2.4.6   usethis_3.2.1   
[25] readxl_1.4.5     PScr_1.1         asaur_0.50       pracma_2.4.4    

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1    Rmpfr_1.1-1         farver_2.1.2        S7_0.2.0            fastmap_1.2.0       digest_0.6.37       timechange_0.3.0   
 [8] lifecycle_1.0.4     ellipsis_0.3.2      magrittr_2.0.4      rlang_1.1.6         tools_4.5.1         data.table_1.17.8   knitr_1.50         
[15] ggsignif_0.6.4      labeling_0.4.3      pkgbuild_1.4.8      RColorBrewer_1.1-3  pkgload_1.4.1       abind_1.4-8         withr_3.0.2        
[22] numDeriv_2016.8-1.1 grid_4.5.1          xtable_1.8-4        future_1.67.0       MASS_7.3-65         globals_0.18.0      scales_1.4.0       
[29] cli_3.6.5           crayon_1.5.3        generics_0.1.4      remotes_2.5.0       rstudioapi_0.17.1   future.apply_1.20.0 km.ci_0.5-6        
[36] tzdb_0.5.0          sessioninfo_1.2.3   cachem_1.1.0        parallel_4.5.1      cellranger_1.1.0    survMisc_0.5.6      vctrs_0.6.5        
[43] Matrix_1.7-3        carData_3.0-5       car_3.1-3           hms_1.1.3           rstatix_0.7.3       Formula_1.2-5       listenv_0.9.1      
[50] parallelly_1.45.1   glue_1.8.0          codetools_0.2-20    stringi_1.8.7       gtable_0.3.6        gmp_0.7-5           pillar_1.11.1      
[57] lava_1.8.1          R6_2.6.1            KMsurv_0.1-6        evaluate_1.0.5      lattice_0.22-7      backports_1.5.0     memoise_2.0.1      
[64] broom_1.0.10        MatrixModels_0.5-4  gridExtra_2.3       xfun_0.54           fs_1.6.6            zoo_1.8-14          pkgconfig_2.0.3 



## Versioning

The results presented in the manuscript correspond to **Release v1.0** of this repository.
