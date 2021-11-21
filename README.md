# multimeta
The files and code contained here are associated with Pantel et al. 202x, Metapopulation dynamics of multiple species in a heterogeneous landscape.

# Citations
Pantel, J.H., T. Lamy, M. Dubart, J.-P. Pointier, P. Jarne, and P. David. 202x. Metapopulation dynamics of multiple species in a heterogeneous landscape. Ecological Monographs. In press.

# Table of Contents
*data_for_submission* contains the files included in the manuscript supplemened DataS1 and Data S2. DataS1 is the folder with the data + model without covariates. Each directory contains the raw data files (dry_final.txt, presence.txt, state.txt) needed to run the R script JHP_JAGS.R (which calls the JAGS model model_full.txt. DataS2 is the folder with the data + model including covariates. Each directory contains the raw data files (dry_final.txt, presence.txt, state.txt, colsource.txt, cov.txt, habit.txt, pluvio.txt) needed to run the R script JHP_JAGS_covar.R (which calls the JAGS model model_2_cov.txt).

- |data_for_submission
- | |-DataS1
- | |-DataS2
- | |-MetadataS1.txt
- | |-MetadataS2.txt

*results* contains the files used in DataS1 (for no_covar) and DataS2 (for covar) to produce JAGS results files, and all of the JAGS results files (fjsample.RData). All of the jsample.RData results files contained here are the results obtained and analyzed in Pantel et al. 202x. They are called in the main posthoc analysis script script_V01.R.

-|results
-| |-covar
-| |-no_covar

*code* contains the main R script used for analysis of all of the results of the Bayesian metapopulation model script_v01.R. It also contains the R environment produced at the conclusion of running the entire R script, 21_11_2021.RData. It also contains the Mathematica notebooks used for the metapopulation simulations including different combinations of covariates in the folder *simulation*.

-|-code
-| |-script_v01.R
-| |-21_11_2021.RData
-| |-simulation
-| | |-simulmetapop_V04.nb
-| | |-simulmetapop_V04.01.nb
-| | |-simulmetapop_V04.02.nb
-| | |-simulmetapop_V04.03.nb
-| | |-simulmetapop_V04.04.nb
-| | |-simulmetapop_V05.nb


*data* contains various input files needed throughout the running of script_v01.R. It also contains a folder, *simulation*, which contains the input files needed to run the Mathematica simulations, and contains some empty folders needed to hold more input files needed for the simulations (that are produced in script_v01.R).
-|-data
-| |-Visite_offline.txt
-| |-final_site_list_coordinates.csv
-| |-dry_reduced_modified.txt
-| |-island_site_id.txt
-| |-final_island_list.txt
-| |-simulation
-| | |-species_obs.csv
-| | |-no_covar
-| | | |-parameters
-| | |-parameters

*raw_output* contains many empty folders to hold the results of the figures and tabled produced by script_v01.R. It also contains the actual results of the Mathematica simulation that are included in Pantel et al. 202x (in the folder *simulation* and contains some of the data needed to create the QGIS maps (in the folder *maps*).

-|-raw_output
-| |-Figure_1
-| |-Figure_3
-| |-Figure_4
-| |-Figure_5
-| |-S2_Figure_S3
-| |-S2_Figure_S4
-| |-Table_1
-| |-Table_S5
-| |-Table_S6
-| |-maps
-| | |-GLM_adm_shp
-| | |-Ec_maps
-| | |-pstar_sim_maps
-| |-simulation
-| | |-simulmetapop_V04
-| | |-simulmetapop_V04.01
-| | |-simulmetapop_V04.02
-| | |-simulmetapop_V04.03
-| | |-simulmetapop_V04.04
-| | |-simulmetapop_V05

