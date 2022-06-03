# killifish-meso
## Size-mediated density-dependence in killifish ecotypes

This project looks at the effects of size-mediated density-dependence on two distinct ecotype of Hart’s killifish, Rivulus hartii. The killifish-only (KO) ecotype is adapted to living in the absence of other fish species. The low-predation ecotype (LP) is adapted to living with one other fish species, the Trinidadian guppy. Guppies and killifish engage in intra-guild predation: while juvenile killifish and guppies compete for similar food resources, adult killifish eat guppies, and adult guppies eat killifish eggs and neonates.

Here, we consider the effects of intraspecific competition on KO and LP killifish, i.e. in the absence of guppies. We ran experiments where we manipulated the population density and population size-structure of killifish in artificial	 streams (mesocosms). We collected data on size, growth, fecundity, and offspring size for each fish in the experiment. Each mesocosm contained killifish from a single population (i.e. KO or LP) at either low density (6 fish standard size distribution), high density (12 fish, standard size distribution), large size structure (9 fish, distribution with large mean), or small size structure (9 fish, distribution with small mean). We replicated the study across two (putatively) independent origins of the KO/LP ecotypes: the Guanapo river and the Naranjo river.

In this repo, we use this experimental data to estimate vital rates (i.e. growth, reproduction, offspring size) as a function of the body size of focal individuals and those of its (varying number of) competitors. We assess support for size-structured density dependence vs classical density dependence. We then use these vital rate functions to build integral projection models (IPMs), which describe the change in the number and size-distribution of individuals over time. Finally, we conduct various analyses of these IPMs to determine how density-dependent fitness varies between killifish ecotypes. 


## Scripts:

**PRIMARY SCRIPTS**:

‘killifish_vital_rates.R’

- Run this to call the vital rate scripts, fit the models, and run diagnostics

“killifish_ipms.R”

- Run this to build and analyse IPMs, and (in future) to perform perturbation analyses.

**DATA UPLOAD**:

‘upload_killifish_data.R’ 

- A single function that selects and uploads the relevant data from the mesocosm experiments. Creates the “gdata” data.table for use in vital rates analyses. Can select specific drainages i.e. “naranjo”, “guanapo”, or “both”.

**VITAL RATE MODELS**:

These scripts define the models and priors for the vital rates. All include random effects of channel on growth.

‘nimble_killifish_model_null.R’

- Code for the null model, i.e. with phi and rho set to zero - equivalent to classical density dependent model

‘nimble_killifish_model_full.R’

- Code for the full model, i.e. KO and LP-specific values of phi and rho are estimated

‘nimble_killifish_model_common surface.R’

- Code for the common surface model, i.e. phi and rho are estimated but assumed to be the same for KO and LP

**IPM FUNCTIONS**:

‘ipm_functions.R’

 - Code to create the IPM kernels and iterate IPMs 

## Data objects:

‘all_meso_data.csv’

- the raw mescosm experimental data file

## Model objects:

MCMC samples: These are the stored MCMC chains for the drainage-specific model runs.

For Naranjo:

- ‘samples_null_naranjo.Rda’

- ‘samples_full_naranjo.Rda’

- ‘samples_common_naranjo.Rda’

For Guanapo:

- ‘samples_null_guanapo.Rda’

- ‘samples_full_guanapo.Rda’

- ‘samples_common_guanapo.Rda’

