# Supplementary-Material-for-``Spatial Interference Detection in Treatment Effect Model''
Supplementary codes and data for ``Spatial Interference Detection in Treatment Effect Model''

## Data and software availability
- Previously published data CMFD (China meteorological forcing dataset) and CHAP (ChinaHighAirPollutants) were used. One can refer to https://doi.org/10.11888/AtmosPhys.tpe.249493.file and https://doi.org/10.5281/zenodo.3752465 for accessing and enabling data download or detailed description of the two datasets. These two datasets are public but they should be cited with the guides provided by the reference sites.

## Code for ``Spatial Interference Detection in Treatment Effect Model''
### Overview

1. Directory *Data Generator* contains the main functions for generating datasets used in simulation and application. Specifically, 
- The R script ***Generate_Hebei_Data*** contains the code for generating dataset *Simulations/CMFD_Hebei.RData*.
- The R script ***Generate_Huabei_Data*** contains the code for generating dataset *Applications/CMFD_PM10_Huabei.RData*.

2. Directory *Simulations* contains the main functions for conducting different detection methods and different simulation settings. Specifically,
- The R script ***BootTest*** contains the code for the estimation and global testing proposed in the main paper.
- The R script ***DebiasedTest*** contains the code for conducting global testing by the debiased lasso test.
- The R script ***BiRS_TEM*** contains the detection code for BiRS and step-down.
- The R script ***Knockoff_TEM*** contains the detection code for the knockoff method.
- The R script ***SimuSize_Independent*** contains the code for conducting simulations under ``Independent Design'' to calculate the empirical size.
- The R script ***SimuSize_Correlated*** contains the code for conducting simulations under ``Correlated Design'' to calculate the empirical size.
- The R script ***SimulationTEM_Independent*** contains the code for conducting simulations under ``Independent Design'' to calculate the empirical TPR and FDR.
- The R script ***SimulationTEM_Correlated*** contains the code for conducting simulations under ``Correlated Design'' to calculate the empirical TPR and FDR.
- The R script ***Detection_plot*** contains the code for visualizing the detection results and ploting Figure 1 and 2 in the main paper.
- The R scripts ***SimulationEST_Independent*** and ***SimulationEST_Correlated*** contains the code for computing the post-detection ATEs and ploting Figure 3 in the main paper.
- The R scripts ***SimulationTEM_MF_Independent*** and ***SimulationTEM_MF_Correlated*** contains the code for obtaining the detection results of our method under Mean-Field assumption.
- The R scripts ***SimulationEST_MF_Independent*** and ***SimulationEST_MF_Correlated*** contains the code for computing the post-detection ATEs under Mean-Field assumption and ploting Figure 4 in the main paper.

3. Directory *Applications* contains the main functions for conducting real data analysis on PM10 data in Huabei, China. Specifically, 
- The R scripts ***BootTest*** and ***BiRS_TEM*** are the same as the two codes with the same names in directory *Simulations*.
- The R script ***Application_TEM*** contains the code for conducting real data analysis on the PM10 emission in Huabei, China.
- The R scrtpt ***plot_application*** contains the code for visualizing the detection results and ploting Figure 5 in the main paper. We note here that the backgroud map in the Figure 5 is obtained from ggmap package, but it could be unavaiable sometimes according to the google API service.

### Workflows
#### Generate Data
1. Run ***Data Generator/Generate_Hebei_Data*** to get the R data *CMFD_Hebei.RData* and put it in the directory *Simulations*.
2. Run ***Data Generator/Generate_Huabei_Data*** to get the R data *CMFD_PM10_Huabei.RData* and put it in the directory *Appliations*.

#### Simulation for comparing detection and estimation performance of different methods
1. Please create folders *Simulations/Results_Independent* and *Simulations/Results_Correlated* for saving the simulation results.
2. Run ***Simulations/SimulationSize_Independent.R*** and ***Simulations/SimulationSize_Correlated.R*** for getting the size results.
3. Run ***Simulations/SimulationTEM_Independent.R*** and ***Simulations/SimulationTEM_Correlated.R*** for getting the detection results.
4. Run ***Simulations/SimulationEST_Independent.R*** and ***Simulations/SimulationEST_Correlated.R*** for getting the estimation results under hetergeneous interference.
5. Run ***Simulations/SimulationEST_MF_Independent.R*** and ***Simulations/SimulationEST_MF_Correlated.R*** for getting the estimation results under mean field interference.
6. Run ***Simulations/Detection_Plot*** to get the figures.
   
#### Data application
1. Run ***Applications/Application_TEM.R*** to get the detection results.
2. Run ***Applications/plot_application*** to get the figures in application.
