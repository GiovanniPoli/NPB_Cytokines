# High-dimensional Bayesian semi-parametric models for small samples: a principled approach to the analysis of cytokine expression data
Welcome to the GitHub repository dedicated to the scientific draft titled "High-dimensional Bayesian semi-parametric models for small samples: a principled approach to the analysis of cytokine expression data" by Poli, Argiento, Amedei and Stingo. 
In this repository, you will find the files and data used in the article, as well as the source code used for statistical analysis and generating graphs.

![alt text](https://github.com/GiovanniPoli/NPB_Cytokines/blob/progetto/Files/FIG2_page-0001.jpg?raw=true)

# Guide to the code
- Packages used and version of R are detailed in Section 1.
- Section 2 contains a guide to fully reproduce the results of simulations and analysis on real data. This may require <b>several days of running</b> on a desktop PC. Alternatively, results (i.e., the MCMC chains and simulated data) can be downloaded following procedure described in the last Subsection of Section 2.
- Steps described in Section 3 reproduce all graphs and tables included in the article.
- Checking the reproducibility of the downloaded results is possible following the steps described in Section 4.
- The remaining files in this folder are described in Section 5.

# Section 1: 

The required packages are listed below. 
Most of them are available on `CRAN` and can be easly installed using `install.packages("name-of-package")`.
Full guide to the installation of packages that are not on `CRAN` can be found in the next Subsection.

### Pakages required 


```{r c00, eval=TRUE, message=FALSE, warning=FALSE, echo=TRUE}
library(mclust)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(lme4)
library(rstudioapi)
library(ROCR)
library(devtools)
library(mcclust.ext)
library(BiocManager)
library(limma)
library(statmod)
library(kableExtra)
library(formattable)
library(knitr)
library(tidyverse)
library(reshape2)
library(dplyr)
library(plyr)
library(readr)
library(car)
library(MASS)
library(ggplot2)
library(tidyr)
library(viridis)
library(scales)
library(plotly)
library(gridExtra)
library(dendextend)
library(ggdendro)
library(RColorBrewer)
library(cowplot)
library(latex2exp)
library(sn)
library(hrbrthemes)
library(paletteer)
library(coda)
```

### Not on `CRAN` packages 

3 packages are not available on `CRAN`; i.e., `limma` and `statmode` and  `mcclust.ext`.
The first two can be installed using the `BiocManager`  package.  
The last one must be installed from Github and thus requires `devtools`.


```{r c0111, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
library(BiocManager)
library(devtools)

BiocManager::install("limma")
BiocManager::install("statmode")

devtools::install_github("sarawade/mcclust.ext")
```



# Section 2:  

This section describes the steps to fully reproduce the results of the study, using the original code.
This can be time-consuming.
To reproduce the plots, the results can also be downloaded following the steps described in the last Subsection of this Section.

The steps described in this Section require the installation of packages listed in Section 1.
Computation time depends on the cores available for parallel computation. 
By default the script will use the number of detected logical cores minus 2.

### Reproduce MCMC chains simulations 
Make sure the  `Scripts` and `Data` folders and the other `.R` files  are in the same folder. 
After that, manually (using source button) or via console (using your own correct path) run `1A_source_file_simulations.R`.
This file creates the folders needed to systematically save the chains and will set the working directory as the folder where files are stored. E.g:  
```{r c1_1, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/1A_source_file_simulations.R", echo=TRUE)
```
After that, you can run the simulation scripts. 
Order is not relevant.
```{r c1_2, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("S1_mod_hierarchical.R")
source("S1_mod_hierarchical_ExFDR.R")
source("S1_mod_v3.R")   
source("S1_mod_v3_ExFDR.R")
source("S2_S4_mod_hierarchical_FDR.R")        
source("S2_S4_mod_hierarchical_Unif.R")  
source("S2_S4_mod_v3_FDR.R")
source("S2_S4_mod_v3_Unif.R") 
source("S3_mod_hierachical.R")
source("S3_mod_v3.R")
```
### Reproduce MCMC chains for real data and sensitivity 
Running both the data and simulations scripts is not required.
As in previous section, make sure the `Scripts` and `Data` folders and the other `.R` files are stored in the same folder. 
Then, simply run the script `2A_source_file_realdata`.
This file creates the folders needed to systematically save the chains and will set the working directory as the folder where the files are stored. 

E.g.:
```{r c1_3, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/2A_source_file_realdata.R", echo=TRUE)
```
After that run the script `Real_data_semi_parametric_and_sensitivity.R`.
```{r c1_4, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("Real_data_semi_parametric_and_sensitivity.R")
```
### Download results      
Intermediate results can be downloaded at this [link](https://drive.google.com/drive/folders/1LQyEs6UILxygsH_7WtMAKp40xIdzgbsh?usp=sharing)

It is possible to directly download the shared folder which contains the shared files and the results of those with a long computational time.

Alternatively, you can download only a subset of the folders that store the results for target checks.
In that case, make sure to follow the correct folder structure created using files `1A_source_file_simulations.R` and `2A_source_file_realdata.R` or following the structure of the shared folder.

To reproduce Figure 1 and Table 2 in the manuscript, folders for scenarios results (i.e. `Scenario_I`, `Scenario_II`, `Scenario_III` and `Scenario_IV`) must be downloaded and placed in the same folder as the other files. 

To reproduce Figure 2, Table 3 and Supporting Informations C in the manuscript, the folder that stores the chains for real data (i.e. `Chains_real_data`) must be downloaded and placed in the same folder as the other files. 
We remark that the final setting must be similar to the shared folder, which has been created following the procedures described with this `README.htlm`.

# Section 3:
### Simulations results
Make sure that:

* Folders `Scenario_I`, `Scenario_II`, `Scenario_III` and `Scenario_IV` contain all 4800 chains.
* The structure of those folders recalls the structure of the [shared folder](https://drive.google.com/drive/folders/1LQyEs6UILxygsH_7WtMAKp40xIdzgbsh?usp=sharing).
* All packages are properly installed.

Now simply run manually (via source button) or via console (using your own corect path) the `1B_Post_processing_simulations.R` script. E.g.:

```{r c1_5, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/1B_Post_processing_simulations.R", echo=TRUE)
```

Object `TAB_2` will show the values in Table 2. File `FIG1.PDF` will be saved from object `plot.export` and corresponds to Figure 1 in the manuscript.

### Real data results
Make sure that

* `Chains_real_data` folder contains all  80 chains. 
* The structure of that folder recalls the structure from the  [shared folder](https://drive.google.com/drive/folders/1LQyEs6UILxygsH_7WtMAKp40xIdzgbsh?usp=sharing). 
* All packages are properly installed.

Now simply run `2B_Post_processing_realdata.R` script. E.g.:

```{r c1_6, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/2B_Post_processing_realdata.R", echo=TRUE)
```

Object `TAB_3` will store the values in Table 3. 
File `FIG2.PDF` will be saved from object `plot.export2` and corresponds to Figure 2 in the manuscript. 
Objects `C.1`,`C.2` and `C.3` store the values of the corresponding Tables in the Supporting Information Section C.
Similarly the script will save files `FIG_C1.PDF`, `FIG_C2.PDF`, and `FIG_C3.PDF` that correspond to Figures in Supporting Information Section C.

# Section 4:

The scripts in this subsection are used to check reproducibility of intermediate results; i.e. chains downloaded from the shared folder.
These scripts use the same functions used to create shared folder, and obtain the first sample stored from the  chain using the same simulated data and the same seed.
Note that the process is not immediate since the chains save samples after a burn-in period (125 iterations for simulations, 525 for the real data).   
Reproducing the first sample for all simulations with script `3A_simulations_test_ALL.R` takes about 5 hours.
In contrast, a check on a target folder with `3B_simulations_test_TARGET.R` takes under 10 minutes and checking the real data chains `3C_real_data_test.R` takes about 30 minutes.
Using  `sample = 2500` option would produce identical results. 

### Check results for simulations      

Check that the downloaded folders  `Scenario_I`, `Scenario_II`, `Scenario_III` and `Scenario_IV` are in the same folder of script `3A_simulations_test_ALL.R`.
Then simply run manually (via source button) or via console (using your own corect path) the script.
Object `results` will be a list of 4 `TRUE`.

```{r c1_13, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/3A_simulations_test_ALL.R", echo=TRUE)
```

### Check results for target simulation 

Check that the target downloaded folder is in the same folder of script `3B_simulations_test_TARGET.R` and at least one file is in the corresponding folder. 
Select element of vectors `scenarios`, `var_levels` and `models` (script lines 22-24).
E.g.:
```{r c1_test, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
scenario = scenarios[1]
var      = var_levels[1]
model    = models[1]
```
Without manipulations the script will check the first sample for chains of the Scenario I and  low-variances cases, for the Semi-parametric model with uniform priors.
Always use the structure of the shared folder i.e. `Scenario_I/Mod_v3/low/dowloaded_file_xx.rds`.
Then simply run the script `3B_simulations_test_TARGET.R`.
Object `results` will be a list of `TRUE` elements of length equal to the number of files stored (e.g. 100 if a full folder is downloaded).
```{r c1_14, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/3B_simulations_test_TARGET.R", echo=TRUE)
```

### Check results for real data       

Check that the downloaded folder is in the same folder of script `3C_real_data_test.R`.
Then simply run script.
Object `results` shuold be a list of 80 `TRUE`.
```{r c1_15, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/CytoDPM/folder_to_push/3C_real_data_test.R", echo=TRUE)
```


# Section 5: 
Here we we briefly describe other files in the folder:

- File `functions_v3.0.R` in folder `Scripts` contains main functions used for simulations and analysis with a brief description.
- Files `My_GG_Scripts.R` in folder `Scripts` contains functions for plots.
- `Markdonw_Sensitivity.RMD` produce `Markdonw_Sensitivity.pdf`. The pdf has plots for convergences checks of all real data chains.
