## Code for paper "Predicted Norovirus Resurgence in 2021-2022 Due to the Relaxation of Nonpharmaceutical Interventions Associated with COVID-19 Restrictions in England: A Mathematical Modelling Study"

This is currently under review and available on [MedrXiv](https://www.medrxiv.org/content/10.1101/2021.07.09.21260277v1)

The code provided in this repo was used to generate the results in the above paper. The steps taken were in the following order:

1. Construct a mathematical model for norovirus, based on the model developed by Ben Lopman (see refs in the paper)
- the ODE's are found in "noro_functions_mcmc_clean.R" using the function "ja.multistage.model.ii.seas", the other function includes specific assumptions about births.
- This code was developed starting with a super-useful introduction by Aaron King and Helen Wearing (https://ms.mcmaster.ca/~bolker/eeid/2011_eco/waifw.pdf) 
2. Find best fitting values of the parameter "q" when the model is compared to data from witihin the Harris study. Repeat this for several different models that cover some of the uncertainties we still have about norovirus transmission.
- this was done by making minor changed in the code witihin noro_functions_mcmc_clean.R and 
- polymod data was used during this fitting part, using the data within the package "socialmixr"
3. Use these values for each of the simulations and run through 3 years of simulations to generate the prediction incidence 2019-2023
- The code for this is in "noro_age_structure_sims08_seas_clean.R"
- simulations have been saved using rdat files
4. Plot the simulations using the file "plot_models_noro_24May21_clean.R"
- not that the PHE is not available for use but the models can be run without this data