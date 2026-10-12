This repository contains the code and data used for the MS: 

Where curse of the Pharaoh meets landscape of disgust: host behavior and pathogen release mode jointly shape disease dynamics and spatial distribution

Amélie C. Dolfi _a_, Marie L. J. Gilbertson _a_, Lauren A. White _b_, Wendy C. Turner _c_

_a_ Wisconsin Cooperative Wildlife Research Unit, Department of Forest and Wildlife Ecology, University of Wisconsin-Madison, Madison, WI 53706, USA
_b_ National Socio-Environmental Synthesis Center, University of Maryland, Annapolis, MD, USA
_c_ US Geological Survey, Wisconsin Cooperative Wildlife Research Unit, Department of Forest and Wildlife Ecology, University of Wisconsin-Madison, Madison, WI 53706, USA


In this repository you will find:
In **Code**: sixe R codes
  - Landscape_Parameters.R: The code used to create the four landscape used in the simualtion and to obtain the parameters
  - Simulation_model.R : Source code to run the model
  - Data_Summary.R : The code used to summarised the simulations
  - Satscan_run.R: The code used to run satscan
  - Satscan_summary.R: The code used to summarise the satscan results
  - Analysis_Figures.R: The code used to analysis the data and obtain the results and figures in the MS and the SI

In **Data_input** : 7 files
  -Land1.csv; Land2.csv; Land3.csv; Land4.csv : The files with landscape matrix
  - Parameters_CS.csv; Parameters_OK.csv; Parameters_Sk.csv; The parameters for each release model (CS: Continuous-shedder; OK: Obligate-killer; SK: Shedder-killer)

In **Data_output** : 14 files
  - Numbers_means_CS.csv; Numbers_means_OK.csv; Numbers_means_SK.csv : Summaries at the simulation level
  - Numbers_final_CS.csv; Numbers_final_OK.csv; Numbers_final_SK.csv: Summaries at the parameter set level
  - Unique_clusters_CS.RData; Unique_clusters_OK.RData; Unique_clusters_SK.RData: Files containing all the unique clusters found by satscan
  - All_CS_time_dynamics.RData; All_OK_time_dynamics.RData; All_SK_time_dynamics.RData: files containing the disease dynamic across time
  - HabRisk_6_Sim19_Para2.csv and Satscan_6_Para_2_OK.RData: files containing results information used to produce figure 1 of the MS

In **Model_output_249**: 330 files
- HabRisk_+_Sim+_Para249_+.csv: 300 files containing the model outputs for the parameter set 249
- Inds_+_Para249_+.csv : 15 files containing the output for each individual in the simulation
- Sum_+_Para249_+.csv : 15 files containing the summary of each simulation

In **Satscan_output_249**:
  - Satscan_+_Para249_+.RData: Satscan results summary
