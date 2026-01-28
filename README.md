
# Size structure as kernel mixture models

## Purpose

Inference from size structured data via [kernel mixture models](https://bitbucket.org/autocatalysis/model_size) and [Conditional Autoregressive SpatioTemporal Models](https://github.com/jae0/carstm).

## Rationale

A **size-frequency distribution** of some **system** of interest is usually obtained from the observation of size measurements at various locations and times. To represent the overall size distribution these observations must be summarized in some manner. The usual approach is to simply aggregate the counts at some arbitrary intervals of size. Implicit in this naive aggregation is that theses samples are **representative** for an explicit unit of **time** and **space**. For this to be true, each observation must have an equal amount of influence or weight in representing the system (that is, independently and identically distributed, or *iid*). When samples are completely random, this approach is fine. 

However, in dealing with the abundance and spatial and temporal distribution of organisms that are:

  - selective in their habitat due to metabolic constraints and behavioural traits (aggregation/clustering/schooling) and generally demonstrating ontogenetic shift in such preferences
  - living in areas of dynamically changing environmental and ecosystem conditions (i.e., **viable habitat**) which is incongruent/mismatched with a usually static observation/sampling domain
  - observed/sampled with bias due to certain environments not being observable due to imperfect sampling:
     
    - nets not able to access rocky, heterogenous, or deep environments
    - not able to access areas of high currents or rapid changes in bathymetry (cliffs and rocky protuberances)
    - sampling of one vertical stratum (size of sampling gear) and missing of mobile organisms that can escape by swimming over or around the nets or others that can burrow into sediments or shelf beside rocky outcrops below the nets.
    - sampling mesh too large to capture small organisms

The usual recourse is to assume some operationally simplistic random stratified design (usually based on depth) will account for most/all of these factors. Armed with this shield of presumed unbiased sampling compute naive size frequency distributions and subsequent analyses. 

It is of course, nearly impossible to correctly design or account for these dynamic preferences. When, as is often the case, the controlling factors are associated with environmental variability that itself changes unpredictably by its nature across both time and space, it is nearly impossible to correctly define these strata and obtain an unbiased design.

What this means is that observations may not be of equal importance in describing the system. For example, if an observations occurs with more frequently an area that is not representative of the system then it’s influence would be more elevated than it should be. Similarly, if an observation occurs in an area that is highly representative of the system, then it’s influence should be higher than it may be considered under the assumption of equal influence. In other words, simply adding the observations together becomes inappropriate. 

A model based approach is one possible recourse to adjust and discover the relative importance of each observation in representing the overall system and ultimately reduce some of theses biases as much as possible (Post-stratified weights).

Here we explore the use of kernel mixture modelling to decompose the size structure at the level of each observation event. The modes and distributional parameters of the kernel mixture model are then modelled in space and time using a spatio-temporal CAR model. These are then reconstituted in the sampling domain to provide a model-based size structure, that accounts for these biases as much as possible through the use of environmental covariates as well as spatial and spatiotemporal random effects to absorb the unmeasured biases.

## Main steps

- From simple kernel density representations of size frequency at small area unit scale, determine the magnitudes of  modes. This is done as there may be regional and time-dependent changes in modal sizes (year-classes)

- Using the most frequently encountered modal groups, identify the main growth groups (instars) as a function of maturity and sex and then define a growth model assuming exponential growth by fitting a log-linear model.

- From the above approximate modal growth model, define a latent kernel mixture model with initial guesses centered on these modal growth groups. This process decomposes the observed size structure into approximate composition/representation of **growth groups** (instar, maturity, sex): **alpha** (relative composition of each growth group), **sigma_mean** (standard deviation of sizes in each growth group), and **imode** (the latent modes of each growth group). These are determined from Bayesian kernel mixture models. 

- Temporal and spatial modelling for each growth group of the latent modes, alpha, and sigma_mean using spatial and spatiotemporal CAR models with depth, bottom temperatures as covariates.

- The posterior distributions of each of these parameters are then used to reconstitute/reconstruct a size distribution that can be assessed in any sub-domain, for description and further analysis (such as growth and dynamics). 

and make use of Julia/Turing and R/INLA for modelling.


## More info

[See size structure documentation](docs/size_structured_model_snowcrab.md)

