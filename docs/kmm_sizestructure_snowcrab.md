---
title: "Size structure as kernel mixture models: snow crab"
author: "Jae S. Choi"
date: last-modified
date-format: "YYYY-MM-D"
toc: true
toc-depth: 4
number-sections: true
highlight-style: pygments
editor:
  render-on-save: false
  markdown:
    references:
      location: document
execute:
  echo: true
format:
  html: 
    code-fold: true
    code-overflow: wrap
    html-math-method: katex
    self-contained: true
    embed-resources: true
params:
  year_assessment: 2024
  model_variation: Maximum_likelhood_basic
---
<!-- Preamble


This is a Markdown document ... To create HTML or PDF, etc, run: 

As it runs R and Julia, Quarto is probably the better tool to render. Alternatively, you can run everything directly and then render just the figures by turning off the R and Julia code chunks. To create an HTML via Quarto:

>> cd ~/projects/model_size/scripts

>> make quarto FN=kmm_sizestructure_snowcrab YR=2024 SOURCE=~/projects/model_size/scripts WK=~/projects/model_size/outputs DOCEXTENSION=html PARAMS="-P year_assessment:2024 -P model_variation:Maximum_likelhood_basic"
 
  
Alter year and directories to reflect setup or copy Makefile and alter defaults to your needs.
   
-->

# Size structure as kernel mixture models: background

## Purpose

Inference from size structured data via [kernel mixture models](https://bitbucket.org/autocatalysis/model_size) and [Conditional Autoregressive SpatioTemporal Models](https://github.com/jae0/carstm).


See also, model_covariance::docs::covariance.md


## Rationale

A **size-frequency distribution** of some **system** of interest is usually obtained from the observation of size measurements at various locations and times. To represent the overall size distribution these observations must be summarized in some manner.

If the samples are derived from a truly random sample from the space-time-size domain, then the usual approach of simply aggregating counts is sufficient. This is because these samples are **representative** and each observation has an equal amount of influence or weight in representing the system (that is, independently and identically distributed, or *iid*).

However, in dealing with the abundance and spatial and temporal distribution of organisms that are:

- selective in their habitat due to metabolic constraints and behavioural traits (aggregation/clustering/schooling)
- demonstrate ontogenetic (size and age-based) shifts in such habitat preferences
- living in areas of dynamically changing environmental and ecosystem conditions (i.e., **habitat**) which is incongruent/mismatched with a usually static observation/sampling domain
- observed/sampled with bias due to certain environments not being observable due to imperfect sampling:

  - nets not able to access rocky, heterogenous, or deep environments
  - not able to access areas of high currents or rapid changes in bathymetry (cliffs and rocky protuberances)
  - sampling of one vertical stratum (size of sampling gear) and missing mobile organisms that can escape by swimming over or around the nets or others that can burrow into sediments or shelf beside rocky outcrops below the nets.
  - sampling speed to capture rapidly moving organisms (escapement)
  - sampling mesh too large to capture small organisms (escapement)

The usual recourse is to assume some operationally simplistic random stratified design (usually on depth) and **hope** the stratification will account for most/all of these above factors. Armed with this shield of **presumed** unbiased sampling within a stratum, one then computes naive size frequency distributions, sometimes weighted by the surface area of each stratum and subsequent analyses based upon them.

This latter random-hopium stratified approach is the defacto standard in most settings. This is the case because, it is generally, nearly impossible, to correctly design or account for these dynamic, context and size/age-dependent preferences. When, as is often the case, the controlling factors are associated with environmental variability that itself changes unpredictably by its nature across both time and space, it is indeed a challenge to correctly define these strata and obtain an unbiased design.

What this means is that observations may not be of equal importance in describing the system. For example, if an observations occurs with more frequently an area that is not representative of the system or stratum, then it’s influence would be elevated more than it should be. Similarly, if an observation occurs in an area that is highly representative of the system, then it’s influence would be higher than it may be. In other words, simply adding the observations together becomes inappropriate.

A model-based approach is one possible recourse (post-hoc) to adjust and discover the relative importance of each observation in representing the overall system and ultimately reduce some of theses biases as much as possible.

Here we explore the use of **Kernel Mixture Models (KMMs)** to decompose the size structure at the level of each observation event. The modes and distributional parameters of the kernel mixture model are then modelled in space and time using a statistical spatio-temporal model. These are then reconstituted in the sampling domain to provide a model-based (weighted) size structure, that accounts for these biases as much as possible through the use of environmental covariates as well as spatial and spatiotemporal random effects to absorb the unmeasured biases.

## Notes:

Estimate size structure and related parameters (growth, maturity, etc.).

From simple kernel density representations of size frequency at small area unit scale, determine the magnitudes of  modes. This is done as there may be regional and time-dependent changes in modal sizes (year-classes)

Inference of modal sizes for each sex, instar, maturity group in order to develop a growth model.

A number of methods to explore:

- Direct via counting (**Simple sums** ... traditional)
- Direct via **areal density** (arithmetic and geometric)
- **Modelled solutions**: mostly too large of a problem to use
- Kernel density areal density (arithmetic and geometric)
- Size at maturity
- Modal size groups

---

# Size structure as kernel mixture models: computing

## Base data

Most/all functions for the initial data extraction are found in bio.snowcrab and aegis packages.

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup-R-environment
 
require(bio.snowcrab)

year.assessment = 2024 
p = snowcrab_parameters(
    project_class="carstm", 
    yrs=1999:year.assessment,   
    areal_units_type="tesselation",
    carstm_model_label=  paste( "default", "fb", sep="_" )  # default for fb (fishable biomass)
)

require(ggplot2)
require(data.table)

loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")

# overwrite/ over-ride defaults
p$project_name ="model_size"
p$project_directory = file.path( homedir, "projects", "model_size" ) 
p$data_root = p$project_directory
p$datadir = p$data_root   # all unprocessed inputs (and simple manipulations) ..   #  usually the datadir is a subdirectory: "data" of data_root as in snowcrab.db, .. might cause problems
p$project.outputdir = file.path( p$data_root, "outputs" ) #required for interpolations and mapping
p$modeldir = file.path( p$data_root, "outputs", "modelled" )  # all model outputs

source( file.path(p$project_directory, "R", "sizestructure_db.R") )

rawdatadir = file.path(p$project_directory, "outputs", "size_structure", "modes_kernel_mixture_models_set")

ss_outdir=file.path(p$project.outputdir, "size_structure")


if ( !file.exists(p$project_name) ) dir.create( p$project_name, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$project.outputdir) ) dir.create( p$project.outputdir, showWarnings=FALSE, recursive=TRUE )

# dimensionality and labels set up for carstm
p$dimensionality="space-time"

# stage labels
stages = list(
  mi = c("m|i|05", "m|i|06", "m|i|07", "m|i|08", "m|i|09", "m|i|10",  "m|i|11",  "m|i|12" ),
  mm = c( "m|m|10",  "m|m|11",  "m|m|12" ),
  fi = c("f|i|05", "f|i|06", "f|i|07", "f|i|08","f|i|09", "f|i|10" ), 
  fm = c( "m|m|9",  "m|m|10",  "m|m|11" )
)

# sex codes
# male = 0
# female = 1
# sex.unknown = 2

# # maturity codes
# immature = 0
# mature = 1
# mat.unknown = 2

# mapping background
additional_features = features_to_add( 
    p=p, 
    isobaths=c( 100, 200, 300, 400, 500  ), 
    xlim=c(-80,-40), 
    ylim=c(38, 60) 
)

# Get data and format based upon parameters:

survey_size_freq_dir = file.path( p$annual.results, "figures", "size.freq", "survey")

years = as.character(1996: year.assessment)
regions=c("cfanorth", "cfasouth", "cfa4x")
 

# recreate_polygons = TRUE 
recreate_polygons = FALSE 
pg = sizestructure_db(p=p, "arealunits", redo=recreate_polygons)
dim(pg)
  
nb = attributes(pg)$nb$nbs
au =pg$AUID

# recreate_polygon = TRUE  # only if you want to redo all analyses and recreate all categories
recreate_polygon = FALSE
SS = sizestructure_db( p=p, DS="carstm_inputs", datasource="snowcrab", rawdatadir=rawdatadir, redo=recreate_polygon  )  
 
p$space_name = pg$AUID 
p$space_id = 1:nrow(pg)

p$time_name = as.character(p$yrs)
p$time_id =  1:p$ny

p$cyclic_name = as.character(p$cyclic_levels)
p$cyclic_id = 1:p$nw

p$stages = names(SS$sk)
p$stages = p$stages[ -grep("04", p$stages)] # data density of instar 4 seems to be too sparse  and results in unstable solutions 

p$varsnames = c( "imodes", "sigmasq_mean",  "alpha_mean",  "Nkmm" ) 

p$nposteriors = 2000

xrange = c(8, 170)  # size range (CW)

dx = 2 #  width of carapace with discretization to produce "cwd"

 

```

## Create base data

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: base-data
 
# merge set and individual level data for further processing
M = size_distributions(p=p, toget="base_data", pg=pg, xrange=xrange, dx=dx, outdir=ss_outdir, redo=TRUE)

# save useful data to file for use in Julia
out = list( Y=M, nb=nb, au=au )  # saving areal units in case we do carstm in julia .. for now not used
fn = file.path( p$project_directory, "outputs", "size_structure", "base_data_julia.rdata" )  
save( out, file=fn )
 

```

## Simple sums

Now that we have data ready, we can do some simple tabulations using data.tables (for speed) and produce the usual size-frequency distributions (already computed in main assessment process):

```{r}
#| eval: false
#| output: false
#| echo: false
#| label: base-counts

# tabulate... non-zero counts ... must use add_zeros=TRUE to add them, on the fly
M = size_distributions(p=p, toget="tabulated_data", xrange=xrange, dx=dx, outdir=ss_outdir, redo=TRUE)
M$cwd = as.numeric( as.character( M$cwd) ) 

hist(M[sex==0, cwd])

```

## Numerical Density per unit area (simple direct or crude, geometric and arithmetic)

```{r}
#| eval: false
#| output: false
#| echo: false
#| label: base-counts

# NOTE: already created in snow crab assessment process ...

years_ss = as.character( c(-11:0) + year.assessment )

# Method 1 ..  equivalent to full factorial model without intercept.
# directly compute areal densities (from above tabulations) 
M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years_ss,  outdir=ss_outdir,redo=TRUE)
#NOTE: "crude" is same as "simple_direct" but without pg and unrolled

# take subset in years
M = size_distributions(p=p, toget="simple_direct", xrange=xrange, dx=dx, Y=years_ss , outdir=ss_outdir )

outdir=file.path( survey_size_freq_dir, "direct")

plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
    plot_sex="female", 
    yvar="den",  # den=arithmetic mean density, denl = geometric mean density  
    outdir=outdir 
) 

plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
    plot_sex="female", 
    yvar="denl",  # den=arithmetic mean density, denl = geometric mean density  
    outdir=outdir 
)

plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
    plot_sex="male", 
    yvar="den",   # den=arithmetic mean density, denl = geometric mean density  
    outdir=outdir 
)

plot_histogram_carapace_width( M=M, years=years_ss, regions=regions, 
    plot_sex="male", 
    yvar="denl",   # den=arithmetic mean density, denl = geometric mean density  
    outdir=outdir 
)


if (0) {
    ss = M[ region=="cfanorth" & year== 2017, which=TRUE]

    ggplot( M[ ss ,], aes(cwd, den, fill=mat, colour=sex) ) +
        #geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
        geom_bar(stat = "identity") +
        labs(x="cw", y="density", size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 
}


```

### Hot spots of Numerical densities

Areal densities are simply computed as well, but they need to make sure zero-valued results are included. And biases due to extreme events are corrected modelled or at least, reduced or robustified.

```{r}
#| eval: false
#| output: false
#| echo: false
#| label: spatial-counts

# directly from databases
M = snowcrab.db( DS="set.complete", p=p ) # note depth is log transformed here
setDT(M)

# high density locations  
i = which(M$totno.all > 2*10^5)
H = M[i, .(uid, plon, plat, towquality, dist, distance, surfacearea, vessel, yr, z, julian, no.male.all, no.female.all, cw.mean, totno.all, totno.male.imm, totno.male.mat, totno.female.imm, totno.female.mat, totno.female.primiparous, totno.female.multiparous, totno.female.berried)]

H$log10density = log10(H$totno.all)

library(ggplot2)

cst = coastline_db( p=p, project_to=st_crs(pg) ) 
  
isodepths = c(100, 200, 300)
isob = isobath_db( DS="isobath", depths=isodepths, project_to=st_crs(pg))
isob$level = as.factor( isob$level)
  
plt = ggplot() +
    geom_sf( data=cst, show.legend=FALSE ) +
    geom_sf( data=isob, aes( alpha=0.1, fill=level), lwd=0.1, show.legend=FALSE) +
    geom_point(data=H, aes(x=plon, y=plat, colour=log10density), size=5) +
    coord_sf(xlim = c(270, 940 ), ylim = c(4780, 5200 )) +
    theme(legend.position="inside", legend.position.inside=c(0.08, 0.8), 
        axis.title.x=element_blank(),
        axis.title.y=element_blank())

fn = file.path( p$project.outputdir, "maps",  "map_highdensity_locations.png" )
png(filename=fn, width=1000,height=600, res=144)
    (plt)
dev.off()

```

Direct arithmetic and geometric means, as above are simple. But to account for environmental covariates, a model-based approach is more flexible. Unfortunately the models below **do not work** due to large problem size and corresponding RAM/CPU bottlenecks.

```{r}
#| eval: false
#| output: false
#| echo: false
#| label: model-based-counts

# ---------------------
# method 2: simple linear (gaussian) model via biglm .. too slow to use

    O = size_distributions(p=p, toget="linear_model" , outdir=ss_outdir)
 
    ss = O[ region=="cfanorth" & year== 2017, which=TRUE]
 
    ggplot( O[ ss, ], aes(cwd, den, fill=mat, colour=mat) ) +
        # geom_ribbon(aes(ymin=density_lb, max=density_ub), alpha=0.2, colour=NA) +
        # geom_line() +
        geom_bar(stat = "identity") +
        labs(x="cw", y="density", size = rel(1.5)) +
        # scale_y_continuous( limits=c(0, 300) )  
        theme_light( base_size = 22 ) 
 

# ---------------------
# method 3: poisson model  via biglm .. problem is too large to compute
  
    # too slow to complete
    O = size_distributions(p=p, toget="poisson_glm" ,  outdir=ss_outdir)
 
    outdir=file.path( survey_size_freq_dir, "poisson_glm")
  
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )
 

# ---------------------
# method 4: poisson via inla .. problem is too large to compute
   
    # adjust based upon RAM requirements and ncores
    require(INLA)
    inla.setOption(num.threads= floor( parallel::detectCores() / 2) )
  
    O = size_distributions(p=p, toget="poisson_inla",  outdir=ss_outdir )
 
    outdir=file.path( survey_size_freq_dir, "poisson_inla")
  
    regions = "cfanorth"
    plot_histogram_carapace_width( M=O$P, years=years, regions=regions, 
        plot_sex="male", 
        yvar="N",  # den=arithmetic mean density, denl = geometric mean density  
        outdir=outdir 
    )

# --- 
# method 5: model using CARSTM ?  lots of tedious computations ...


# --- 
# method 6: model using JuliaGLM ? might be faster...



```

So the above modelling attempts do not work (operationally). The trick is to find an approach that will.

Giving up for now to create a size-space-time model ...  but try Julia-GLM

## Size structure in continuous form: Kernel-density based methods

Continuing with size analysis, we can use kernel density estimates of specific components: sex, maturity, year, time of year (season in quarters), region and set (sid).
First we construct kernel density estimates using a bandwidth of 0.025 units on a logarithmic scale. This corresponds to about 4 dx (mm), where dx is the increment width of discretization.

These results are normalized by swept area to provide a density per unit area (1 km$^{-2}$).

**NOTE: Consider moving this to Julia**

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: size-normalization

# Normalization via weighted kernel density data from "base-data"
 
# key defaults that define kernal densities:
np = 512  # # discretizations in fft

xr = round( log(xrange), digits=2 ) 
ldx = diff(xr)/(np-1)  # 0.005988
xvals = seq( xr[1], xr[2], by=ldx )

# bw is on log scale ... approx (log) SD for each interval  
#  data_resolution is ~1 to 2 mm (observation error)
#  but some catgories are wider than others :
#   male imm  = 10 to 100 mm
#   male mat  = 60 to 168 mm
#   female imm = 10 to 45 mm
#   female mat = 35 to 65 mm

# bw = 0.1 # ~ 20 ldx ~ overly smooth
# bw = 0.05  # used for modal analysis
# bw = 0.025  # optimal for sparse data  <<<<<< DEFAULT for histograms >>>>>>
# bw = 0.01 # 2 ldx is too noisy for arithmetic but good for geometric means

bw =list( 
  "0"=list("0"=0.05, "1"=0.05), #male( imm, mat)
  "1"=list("0"=0.04, "1"=0.04 ) #female( imm, mat)
)

# years of interest:
years = as.character( c(1996:year.assessment) )

ti_window=c(-4,4)  # include data from +/1 4 weeks 
sigdigits = 3

lowpassfilter=0.001
lowpassfilter2=0.001

# strata = "smryzt"  # with temp and depth strata too  ...
strata = "yasm"    # moving average in space and time .. sa-weighted kernel density by sid , sex, mat (with au and quarter)

redo =TRUE
# redo = FALSE

# compute kernel density estimates
M = size_distributions(p=p, toget="kernel_density_weighted", 
  bw=bw, np=np, ldx=ldx, xrange=xrange, 
  Y=years, strata=strata, pg=pg, sigdigits=sigdigits, ti_window=ti_window,  
  outdir=ss_outdir, redo=redo ) 
 

nmin = 3
bw2 =list( 
  "0"=list("0"=0.03, "1"=0.03 ), #male( imm, mat)
  "1"=list("0"=0.03, "1"=0.03 ) #female( imm, mat)
)

# identify modes from kernel density estimates when there is enough data
M = size_distributions(p=p, toget="kernel_density_modes", 
  strata=strata, bw=bw2, np=np, 
  Y=years, pg=pg, sigdigits=sigdigits, n_min=nmin,
  lowpassfilter=lowpassfilter, lowpassfilter2=lowpassfilter2,  outdir=ss_outdir,
  redo=redo )

    
# find most frequent peaks and valleys .. if recreating this, stepping through the function is probably best 
# as there are many decisions that need to be made   

bw3 =list( 
  "0"=list("0"=0.03, "1"=0.03 ), #male( imm, mat)
  "1"=list("0"=0.03, "1"=0.03 ) #female( imm, mat)
)

mds = size_distributions(p=p, toget="modal_groups", strata=strata, bw=bw3, np=np, ldx=ldx, 
  M=M, 
  sigdigits=sigdigits, lowpassfilter2=lowpassfilter2, outdir=ss_outdir, redo=redo )


```

Using these KD estimates with associated weights, we can flexibly aggregate across any strata ...

Here we use a geometric mean after adding a small positive valued offset (smallest non-zero-value) that mimics the scale (magnitude) of observation errors. Here is amounts to 100 individuals /km^2, below which density is not defined.

Using the kernel density approach, we can compute on the normalized densities to identify size modes.

The algorithm is simply to go to every polygon and identify samples from it and potentially surrounding areas and times and potentially surrounding time windows and then compute modes and troughs from first and second order differentials of the smoothed kernel densities. Currently, the default is to use a small local spatial group with a time window (strata = "yasm") ...

However, we can use other factors to stratify (strata = "smryzt"), such as depth (zlevels are left bounds in meters) and temperature (tlevels are left bounds of temperature in Celcius). Here simple, low/med/high are implemented, but not used as the above small-area-based methods work well.

```code
Useful figures (generated above):

[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_imm_allsolutions.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_mat_allsolutions.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_mat_all_solutions.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_female_mat_all_solutions.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_imm.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_mat.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_female_imm.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_female_mat.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_female_growth_trajectory_empirical.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_growth_trajectory_empirical.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_female_growth_trajectory_empirical_tweaked.png"
[1] "/home/jae/bio.data/bio.snowcrab/assessments/2024/figures/size.freq/survey/modes_male_growth_trajectory_empirical_tweaked.png"
```

These are the results:

Female growth of modes (t vs t-1)

```output

summary( lm(formula = logcw ~ logcw0 * mat, data = mds[sex == "f", ], 
    na.action = "na.omit") )

Residuals:
     Min       1Q   Median       3Q      Max 
-0.06112 -0.01220  0.00077  0.01119  0.05588 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   0.2558     0.0819    3.12   0.0206
logcw0        1.0167     0.0248   40.93  1.4e-08
matm         -1.7146     0.3592   -4.77   0.0031
logcw0:matm   0.4026     0.0929    4.33   0.0049

Residual standard error: 0.0367 on 6 degrees of freedom
  (4 observations deleted due to missingness)
Multiple R-squared:  0.997,	Adjusted R-squared:  0.996 
F-statistic:  779 on 3 and 6 DF,  p-value: 3.67e-08
 
```

Male growth of modes  (t vs t-1)

```output

summary( lm(formula = logcw ~ logcw0 * mat, data = mds[sex == "m", ], 
    na.action = "na.omit") )

Residuals:
     Min       1Q   Median       3Q      Max 
-0.07578 -0.01237  0.00164  0.01968  0.05031 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.16958    0.07278    2.33    0.053
logcw0       1.03569    0.02079   49.81  3.4e-10
matm        -0.17502    0.37497   -0.47    0.655
logcw0:matm  0.00467    0.09023    0.05    0.960

Residual standard error: 0.04 on 7 degrees of freedom
  (5 observations deleted due to missingness)
Multiple R-squared:  0.998,	Adjusted R-squared:  0.997 
F-statistic: 1.05e+03 on 3 and 7 DF,  p-value: 1.2e-09

```

Female mature simple model: cw vs instar

```output
 
summary( lm(formula = logcw ~ instar, data = mds[sex == "f" & mat == "m", 
    ], na.action = "na.omit") )

Residuals:
      1       2       3 
-0.0253  0.0507 -0.0253 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  -0.0157     0.4402   -0.04    0.977
instar        0.4090     0.0439    9.32    0.068

Residual standard error: 0.0621 on 1 degrees of freedom
Multiple R-squared:  0.989,	Adjusted R-squared:  0.977 
F-statistic: 86.9 on 1 and 1 DF,  p-value: 0.068
 
```

Female immature simple model: cw vs instar

```output

summary( lm(formula = logcw ~ instar, data = mds[sex == "f" & mat == "i", 
    ], na.action = "na.omit") )

Residuals:
     Min       1Q   Median       3Q      Max 
-0.05154 -0.00868 -0.00699  0.00403  0.06546 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.12368    0.04548    24.7  2.9e-07
instar       0.31059    0.00593    52.4  3.2e-09

Residual standard error: 0.0355 on 6 degrees of freedom
  (3 observations deleted due to missingness)
Multiple R-squared:  0.998,	Adjusted R-squared:  0.997 
F-statistic: 2.75e+03 on 1 and 6 DF,  p-value: 3.24e-09

```

Male mature simple model: cw vs instar

```output

summary( lm(formula = logcw ~ instar, data = mds[sex == "m" & mat == "m", 
    ], na.action = "na.omit") )

Residuals:
       1        2        3 
 0.00283 -0.00567  0.00283 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.65417    0.05413    12.1   0.0526
instar       0.33550    0.00491    68.4   0.0093

Residual standard error: 0.00694 on 1 degrees of freedom
  (1 observation deleted due to missingness)
Multiple R-squared:     1,	Adjusted R-squared:     1 
F-statistic: 4.67e+03 on 1 and 1 DF,  p-value: 0.00931

```

Male immature simple model: cw vs instar

```output

summary( lm(formula = logcw ~ instar, data = mds[sex == "m" & mat == "i", 
    ], na.action = "na.omit") )

Residuals:
     Min       1Q   Median       3Q      Max 
-0.04238 -0.02548 -0.00584  0.02079  0.08236 

Coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)   1.1912     0.0446    26.7  2.6e-08
instar        0.2994     0.0053    56.5  1.4e-10

Residual standard error: 0.0411 on 7 degrees of freedom
  (3 observations deleted due to missingness)
Multiple R-squared:  0.998,	Adjusted R-squared:  0.997 
F-statistic: 3.19e+03 on 1 and 7 DF,  p-value: 1.43e-10
 
```



## Identify modal groups

Using the most frequently encountered modal groups in the above, identify the main size groups (instars) as a function of maturity and sex and then define a growth model assuming exponential growth by fitting a log-linear model.

### Estimation of abundance at modes: knifed edged cuts

Density and variability estimation via Modal Kernel Mixture Models (KMM) is done in Julia: See projects/model_size/kmm_snowcrab.md for more info. NOTE: this approach is still too slow to use operationally at each set level -- but is viable annually. But that would prevent further covariate modelling.

Here instead, to estimate areal density, we use a knife-edged cut at midpoints between modal groups. This is imperfect as large groups can bleed into adjacent smaller groups. But it is consistent and simple.  We will use this for modelling: either via

- carstm (TODO -- that is to update data source as modelling approach is already complete )
- stmv (TODO -- that is to update data source as modelling approach is already complete )
- abm (see: projects/model_agent/julia/ )

### Estimation of abundance at modes: Classify data using KMMs

From the above approximate modal growth model, define a latent kernel mixture model with initial guesses centered on these modal growth groups. This process decomposes the observed size structure into approximate composition/representation of **growth groups** (instar, maturity, sex): **alpha** (relative composition of each growth group), **sigma_mean** (standard deviation of sizes in each growth group), and **imode** (the latent modes of each growth group). These are determined from Bayesian kernel mixture models.

For more details, see: https://journal.r-project.org/articles/RJ-2023-043/

First save a data dump for all size data to be read into Julia: initial estimate of main modes via kernel density estimation at a fine space-time resolution from R (**size_distributions(p=p, toget="base_data", ...)**), and "mds" naive modal estimates.

### Kernel mixture modelling via Julia

First bring in the base data.


#### Data

```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: size-import-data-julia

using DrWatson

quickactivate( joinpath(homedir(), "projects", "model_size") )


include( scriptsdir( "kmm_startup.jl" )) # load libs and local functions

include( srcdir( "kmm_functions.jl" ))

# install_required_packages(pkgs)  # in case you need to install packages

   
nmin=3
sexes = ["m", "f"]
mats = ["i", "m"]
yrs = 1996:2024
regions = ["cfanorth", "cfasouth", "cfa4x"  ]

size_data_directory = projectdir( "outputs", "size_structure" )

include( joinpath( homedir(), "projects", "model_covariance", "src", "car_functions.jl" ))   
 
# load R-saved data:
Y, nb, aus = size_structured_data()  # "base_data_julia" in R
 
xrange = (8, 170)   # mm
Y = Y[ (Y.logcw .>= log(xrange[1])) .&  (Y.logcw .<= log(xrange[2])) ,:]

```

#### Naive approach

Compute modes from Kernel Density Estimated smooths (weighted by swept area) for each areal-unit (small area polygon with nearest neighbours) and time-unit (+/- 4 weeks) level mode estimates are then aggregated across all areas and times

```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: modes-naive

# construct modes for each au and time windows 
# only need to compute if adding new year
# kdem = kernel_density_estimated_modes(Y; toget="compute", nb=nb, aus=aus, yrs=2024, nmin=nmin ) 

# summary of above modes from another kde of the modes 
mds = kernel_density_estimated_modes(Y; toget="summary", yrs=yrs ) 


```

#### Node with growth constraints

```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: modes-constrained

growth with modes from kde



```

### Classify with KMM

#### Basic model: a **determinsitic** kernel mixture model (KMM) with n_imodes gaussian components

For each area of interest (CFA), classify instar relative numbers.

Estimate weights of size/stage structure (sigma, SD) at the scale of large fishing area (CFAs) and year.

This model ignores environmental covariates and does not permit instar modes to change, ... so could be biased.


```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-basic-model

# run mcmc sampling model size distributions across strata (predicted modes are deterministically fixed)
# aggregated by region and year

outdir = projectdir( "outputs", "size_structure", "kernel_mixture_models", "deterministic_modes_region" )

mds = instar_inference() 
 
kmm_chain( Y, mds; modeltype="deterministic", yrs=yrs, sexes=sexes, mats=mats, regions=regions, savedir=outdir, nmin=10 ) 

# example: show decomposition

  # first load model chain of size distributions across specified strata 
  region = "cfanorth"
  yr = 2004
  sex = "m"
  mat = "i"

  chain = kmm_chain_load(  region=region, yr=yr, sex=sex, mat=mat, savedir=outdir )  
  showall( summarystats(chain) )

  imodes = mds[ (mds.sex.==sex) .& (mds.mat.==mat), :predicted ]
  
  # extract posterior samples of size distributions 
  os = kmm_samples(chain, 5000, imodes=imodes, modeltype="deterministic" )  
  
  xr = round.( (minimum(Y.logcw), maximum(Y.logcw)), digits=2 ) # (2.14, 5.14)
  xbins = range(xr[1], stop=xr[2], step=0.05 )   # for plots

  pl = histogram( skipmissing(vec(os)), bins=xbins, color=:white )
  for i in 3:8
    histogram!( skipmissing( os[:, i]), bins=xbins)
  end
  pl
  
  # create a summary and save to file
  # is also saving to RData: using RCall and add flag save_RDS=true
  ks = kmm_summary( mds=mds, yrs=yrs, sexes=sexes, mats=mats, regions=regions, toget="compute", modeltype="deterministic", savedir=outdir )
   
  ks = kmm_summary(; toget="saved_results", yrs=yrs, savedir=outdir ) 

  pl = StatsPlots.scatter(ks.imodes, sqrt.(ks.sigmasq_mean))

  # majority of SD of each mode:
  # < 0.16; all <0.25 (log cw)
   
```

#### **Latent** KMM with approximate modes and observation error

For each sampling location, classify instar relative numbers.

The results can be aggregated to higher levels including full domain via CAR.

Here, the determind modes are starting points (rather than a fixed constraint as in the previous analysis) and the actual latent modes are estimated from the data.

This takes about 6 hours for 1 year of data ...

```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-latent-model
 
outdir = projectdir( "outputs", "size_structure", "kernel_mixture_models", "latent_modes_set" )

nmin = 5  # lower as this is computed at the level each trawl sampling ("set")

# run mcmc sampling model size distributions across strata 
yrs = 2014:2024

xrange = (8, 170)   # mm
mds = instar_inference() 
mds = mds[ (mds.predicted .>= log(xrange[1])) .& (mds.predicted .<= log(xrange[2])) ,:]
Y = Y[ (Y.logcw .>= log(xrange[1])) .&  (Y.logcw .<= log(xrange[2])) ,:]

kmm_chain( Y, mds; modeltype="latent", yrs=yrs, sexes=sexes, mats=mats, nmin=nmin, savedir=outdir )   
 

# extract as an example, a single solution: kmm_S04072000~1_m_i_2000.jl2
  chain = kmm_chain_load( sid="S01102023~6", yr=2023, sex="m", mat="m", savedir=outdir  ) 
  chain = kmm_chain_load(fnout = joinpath(outdir, "2023", "kmm_S01102023~6_m_m_2023.jl2" ))
  os =  kmm_samples(chain, 5000,  modeltype="latent") # obtain samples for it

  imodes = mds[ (mds.sex.=="m") .& (mds.mat.=="m"), :predicted ]
  n_imodes = length(imodes)

  xr = round.( (minimum(Y.logcw), maximum(Y.logcw)), digits=2 ) # (2.14, 5.14)
  xbins = range(xr[1], stop=xr[2], step=0.05 )   # for plots

  # example decomposition of modes (probabilities of each instar)
  pl=plot()
  for i in 1:n_imodes
    pl = histogram!( pl, skipmissing( os[:, i]), bins=xbins, fill=true )
  end
  pl

  # actual distribution
  i = findall( (Y.sid.=="S01102023~6") .&  (Y.sex.=="m")  .& (Y.mat.=="m") )
  Y[i,:]
  plot()
  histogram!( log.(Y[i,:cw]) )

  # create a summary and save to file
  # to save Rdata: using RCall  and flaf save_RDS=true
  # NOTE: <$> activates R REPL <backspace> to return to Julia
  for yr in yrs
    print(yr); print("\n")
    ks = kmm_summary(mds=mds, yrs=yr, sexes=sexes, mats=mats, toget="redo", modeltype="latent", savedir=outdir, save_RDS=true )
  end

  # export to R and jl2 complete
  
  ks = kmm_summary( yrs=2024, toget="saved_results", modeltype="latent", savedir=outdir ) # or if already completed

  pl = StatsPlots.scatter(ks[1].instar, exp.(ks[1].imodes))


  
```

If it ever becomes fast enough to handle the GP modelling of covariates ...

It might be already for a simple spatial model but using R/INLA as infrastructure is already there so using that (now completed).

```{julia}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-guassian-process

# spatial parameters
nAU = size( nb, 1 )    # no of au
nAU_float = convert(Float64, nAU)
auid = parse.(Int, aus) 
node1, node2, scaling_factor = nodes( nb ) # pre-compute required vars from adjacency_matrix outside of modelling step
nnodes = length(node1)

outdir = joinpath( project_directory, "size_structure", "modes_kernel_mixture_models_set" )
 
 
function icar_form(theta, phi, s_sigma, s_rho)
  dphi_s = phi[node1] - phi[node2]
  dot_phi = dot( dphi_s, dphi_s )
  sum_phi = sum(phi) 
  re = s_sigma .* ( sqrt.(1 .- s_rho) .* theta .+ sqrt.(s_rho ./ scaling_factor) .* phi )  
  return dot_phi, sum_phi, re
end


Turing.@model function kmm_spatial(x, imodes, n_imodes, N, sd_imodes )
  
  # icar (spatial effects)
  s_theta ~ filldist( Normal(0.0, 3.0), nAU, nz)  # unstructured (heterogeneous effect)
  s_phi ~ filldist( Normal(0.0, 3.0), nAU, nz) # spatial effects: stan goes from -Inf to Inf .. 
  s_sigma ~ filldist( LogNormal(0.0, 1.0), nz) ; 
  s_rho ~ filldist(Beta(0.5, 0.5), nz);
  
  for z in 1:nz
    # spatial effects (without inverting covariance)  
    dot_phi, sum_phi_s, convolved_re_s = icar_form( s_theta[:,z], s_phi[:,z], s_sigma[z], s_rho[z] )
    Turing.@addlogprob! -0.5 * dot_phi
    sum_phi_s ~ Normal(0, 0.001 * nAU_float);      # soft sum-to-zero constraint on s_phi)
   
    lambda = exp.( X * beta +  convolved_re[auid] + log_offset )
    @. y ~ Poisson( lambda );
  end

  
  # Kernel Mixture model with modes as **latent** pre-specified components
  # alpha = concentration parameter of results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
  sigmasq ~ filldist( truncated( InverseGamma(5.0, 0.1), 0.0, sd_imodes), n_imodes)  # variance prior 
  kernels = map( i -> Normal( imodes[i], sqrt(sigmasq[i]) ), 1:n_imodes ) 
  alpha ~ Dirichlet(n_imodes, 1.0)
  mixdist = MixtureModel(kernels, alpha)
  x ~ filldist(mixdist, N)
 
end


# run mcmc sampling model size distributions across years
   
# log_offset (if any)
data, N, modes_mds, n_imodes, sd_imodes = size_structured_data_subset( Y, mds; yrs=yr, sexes=sex, mats=mat, ids=id)
    N < nmin && continue 
    imodes = modes_mds[:, :predicted_mean]
    M = kmm(data, imodes, n_imodes, N, sd_imodes )
    # chain = sample(M, sampler, nsamples; max_depth=8, init_ϵ=0.0125)  
    chain = sample(M, sampler, MCMCThreads(), nsamples, n_chains; max_depth=8, init_ϵ=0.0125)  
  
kmm_spatial_models( Y, nb, aus, yrs=yrs, sexes=sexes, mats=mats ) 

function kmm_spatial_models( Y, yrs=yrs, sexes=sexes, mats=mats ) 
 ... bym()
end
```

## Spatiotemporal modelling of KMM parameters

Temporal and spatial modelling for each growth group of the latent modes, alpha, and sigma_mean using spatial and spatiotemporal CAR models with depth, bottom temperatures as covariates.

Spatial modelling by year for associated stats using CAR and sa as offsets to give densities.

A spatial CAR for each time slice (year) to spatially smooth local
population size structure.

Time modeling (e.g. AR, RW, etc ) requires lags and as growth is not
deterministic, adding that structure here is probably not a good choice

better left for growth modeling in a separate process stage later.

To consider:

-- assume no movement across au's to compute au-specific sampling selectivity
-- compute growth curves based upon location of latent modes in subsequent time periods
-- cartsm above
-- carstm modes -- environmental challenges -> lower mode values
-- carstm sigmas -- environmental variability -> sigma
-- carstm alphas
    -- alphas are probabilities of latent modes in each sampling location .. so integral of "population" proportions (should sum to 1.0) ..
    -- track peaks ... ??

Prepare inputs

```{r}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-carstm

## NOTE:: Nkmm is specific to each (year, sex, mat, set)  .. but not stage (instar) 

for (stage in p$stages ) {
  # stage = "f|m|11"
  M0 = SS$sk[[stage]]

  for (variabletomodel in p$varsnames) {

    p$carstm_model_label = paste("kmm", stage, variabletomodel, sep="_")
    p$variabletomodel = variabletomodel

    M = rbind(M0[ is.finite(get(variabletomodel)), ], SS$M, fill=TRUE, use.names=TRUE)
    # variabletomodel = "imodes"

    # keep model as simple as possible .. reflecting environmental covariates only
    p$formula = as.formula( paste( p$variabletomodel,
        ' ~ 1',
        ' + f( time, model="iid",  hyper=H$iid  ) ',
        ' + f( cyclic, model="iid", hyper=H$iid  )',
        ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2 ) ',
        ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2 ) ',
      #  ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      #  ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      #  ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
        # ' + f( inla.group( pca3, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
        ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2 ) ',
        ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="iid", hyper=H$iid)) '
    ) )

    p$family =  "gaussian"  
    weights = M$Nkmm
    weights[ is.na(weights)] =1
    weights[ weights==0] =1

    if (variabletomodel=="Nkmm") {
      p$formula = update.formula( p$formula, . ~ . + offset( data_offset) ) # CARSTM does log-transformation internally  
      p$family =  "poisson"
      weights = rep(1, nrow(M))
    }

    res = NULL
    res = carstm_model( p=p, sppoly=pg, data=M,
      nposteriors = p$nposteriors, 
      toget = c("summary", "random_spatial", "predictions"), 
      posterior_simulations_to_retain=c( "summary", "random_spatial", "predictions"), 
      num.threads="4:2",  # very memory intensive ... serial process
      compress=TRUE, 
      # debug="extract",
      control.inla = list(  cmin=0, h=0.05, diagonal=1e-5, fast=FALSE, improved.simplified.laplace=TRUE, restart=1 ), 
      weights=weights,
      verbose=TRUE   
    )

  
    # res = carstm_model(  p=pN, DS="carstm_summary" )  # parameters in p and summary
    res_vars = c( names( res$hypers), names(res$fixed) )
    for (i in 1:length(res_vars) ) {
      o = carstm_prior_posterior_compare( res, vn=res_vars[i] )  
      dev.new(); print(o)
    }   

  }
}

# extract and reformat posterior simulations
psims = list(m=list(), f=list() )
for (stage in p$stages ) {
  
  if (grepl("^f", stage)) sex = "f"
  if (grepl("^m", stage)) sex = "m"
  
  for (variabletomodel in p$varsnames) {
    if (!exists(variabletomodel, psims[[sex]] )) psims[[sex]][[variabletomodel]] = list()
    p$carstm_model_label = paste("kmm", stage, variabletomodel, sep="_")
    p$variabletomodel = variabletomodel
    res = NULL
    res = carstm_model( p=p, DS="carstm_samples", sppoly=pg )
    rsims = NULL
    if (!is.null(res)) rsims = res$predictions
    psims[[sex]][[variabletomodel]][[stage]] = rsims
  }
}

# save by sex and variablename (split as they are large file sizes) for further analysis (12_* : growth and instar dynamics)
for (sex in c("f", "m")) {
for (variabletomodel in p$varsnames) {
  fn = file.path( p$modeldir, paste("psims", "_", sex, "_", variabletomodel, ".RDS", sep="") )
  out = psims[[sex]][[variabletomodel]]
  read_write_fast( data=out, file=fn )
}}
psims = out = NULL ; gc()


# Now some maps and figures

for (stage in p$stages) {

    M0 = SS$sk[[stage]]
  
    for (variabletomodel in p$varsnames) {

      p$carstm_model_label = paste("kmm", stage, variabletomodel, sep="_")
      p$variabletomodel = variabletomodel

      res = NULL
      res = carstm_model( p=p,  DS="carstm_summary" )  # parameters in p and direct summary
  
      if (is.null(res)) next()
            
      oeffdir = file.path(p$data_root, p$carstm_model_label, "figures")
      fn_root_prefix = "Predicted_numerical_abundance"
      carstm_plot_marginaleffects( p, oeffdir, fn_root_prefix ) 
  

      # maps of some of the results
      outputdir = file.path(p$data_root, p$carstm_model_label, "maps" )
      carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="random_spatial", probs=c(0.025, 0.975) ) 
      carstm_plot_map( p, outputdir, fn_root_prefix , additional_features, toplot="predictions", probs=c(0.1, 0.9)) 


      # posterior predictive check
      MM = rbind(M0[ is.finite(get(variabletomodel)), ], SS$M, fill=TRUE, use.names=TRUE)
      iobs = which(MM$tag == "observations")
      vn = variabletomodel

      fit = NULL; gc()
      fit = carstm_model( p=p, DS="modelled_fit") #,  sppoly = pg )

      pld = data.table(
        observed = MM[iobs , ..vn] , 
        fitted = fit$summary.fitted.values[["mean"]] [iobs]
      )
      names(pld) = c("observed", "fitted")
      anno1 = paste( "Pearson correlation: ", round( cor( pld$fitted, pld$observed, use="pairwise.complete.obs" ), 3))
      # cor( fitted, observed, use="pairwise.complete.obs", "spearman" )

      out = ggplot(pld, aes(x =  observed, y = fitted )) +
        geom_abline(slope=1, intercept=0, color="darkgray", lwd=1.4 ) +
        geom_point(color="slategray") +
        labs(caption=anno1, color="slateblue") +
        theme( plot.caption = element_text(hjust = 0, size=12 ) )# move caption to the left 
  
      outputdir = file.path( p$modeldir, p$carstm_model_label )
      fn = file.path(outputdir, paste("posterior_predictive_check_", vn, ".png", sep="") )
      ggsave(filename=fn, plot=out, device="png", width=12, height = 8)
      print(out)  

      fit  = MM = NULL; gc()
  }
}

```

## Reconstruct size structure

The posterior distributions of each of these parameters are then used to reconstitute/reconstruct a size distribution that can be assessed in any sub-domain, for description and further analysis (such as growth and dynamics).

Reconstruct size structure and associated parameters from KMM samples

From the KMM decompostion and spatiotemporal interpolations, we can now reconstruct a few biologically informative features:

- growth increment across cohort as a function of location
- growth rate parameters as a function of location
- size at maturity as a function of location
- min value of selectivity coefficients grwven constant natural mortality (as a function of location)
- etc.

### Growth increments

A simple kernel density approach across was used to get a first approximation of growth patterns by identifying frequently occurring modes from size frequency data at areal unit and annual timescales (above).

Now we can try to make this more precise by tracking instars. The kernel mixture modelling decomposes the size density and estimates latent modes of size withing each approximate range of sizes associated with an instar. The assumption is that genetically determined growth trajectory is strong and approximated by  the population level kernel density estimates, and that  variations/deviations from it these are due to variations in the local environmental conditions (resource availability, stress, etc).

These latent modal size estimates, W = log(carapace width; mm).

Load results from previous step and plot some histograms (density):

```{r}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-growth-stanzas


# histograms of modes
for (sex in c("f", "m")) {
  # sex = "m"

    variabletomodel = "imodes"

    fn = file.path( p$modeldir, paste("psims", "_", sex, "_", variabletomodel, ".RDS", sep="") )
    psims = aegis::read_write_fast( fn )
  
    # flatten for plotting
    sss = rbindlist( lapply( psims, function(x) as.data.table(x)), idcol=TRUE)

    x11()
    o = ggplot( sss, aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

    x11()
    o = ggplot( sss[time=="2020" ,], aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density_2020", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

    x11()
    o = ggplot( sss[time=="2021",], aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density_2021", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

    x11()
    o = ggplot( sss[time=="2020" & space=="500",], aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density_2020_500", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

    x11()
    o = ggplot( sss[time=="2020" & space=="1",], aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density_2020_1", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

    x11()
    o = ggplot( sss[time=="2021" & space=="1",], aes( value )) +
      geom_density(aes(fill =.id ), color=NA, alpha = 0.5,  position = "identity") +
      theme_light( base_size = 22) + 
      theme( legend.title=element_blank()) +
      labs(x="log(cw; mm)", y="Density", size = rel(1.5)) 
    fn = file.path( p$modeldir, "figures", paste("density_2021_1", "_", sex, "_", variabletomodel, ".png", sep="") )
    ggsave( fn, o, width=8, height=4, units="in", dpi=320, create.dir=TRUE )

}


# growth increments
compute_growth_increments = FALSE
if (compute_growth_increments) {
  sex = "m"; variabletomodel = "imodes"
  psims = aegis::read_write_fast( file.path( p$modeldir, paste("psims", "_", sex, "_", variabletomodel, ".RDS", sep="") ) )

  y0 = 1:(p$ny-1)
  y1 = 2:p$ny

  grw = list()  # label is the end state

  grw[["m|i|06"]] = psims[["m|i|06"]][, y1, ] - psims[["m|i|05"]][, y0,] 
  grw[["m|i|07"]] = psims[["m|i|07"]][, y1, ] - psims[["m|i|06"]][, y0,] 
  grw[["m|i|08"]] = psims[["m|i|08"]][, y1, ] - psims[["m|i|07"]][, y0,] 
  grw[["m|i|09"]] = psims[["m|i|09"]][, y1, ] - psims[["m|i|08"]][, y0,] 
  grw[["m|i|10"]] = psims[["m|i|10"]][, y1, ] - psims[["m|i|09"]][, y0,] 
  grw[["m|i|11"]] = psims[["m|i|11"]][, y1, ] - psims[["m|i|10"]][, y0,] 
  grw[["m|i|12"]] = psims[["m|i|12"]][, y1, ] - psims[["m|i|11"]][, y0,]

  grw[["m|m|10"]] = psims[["m|m|10"]][, y1, ] - psims[["m|i|09"]][, y0,] 
  grw[["m|m|11"]] = psims[["m|m|11"]][, y1, ] - psims[["m|i|10"]][, y0,] 
  grw[["m|m|12"]] = psims[["m|m|12"]][, y1, ] - psims[["m|i|11"]][, y0,] 


  # growth increments female
  sex = "f"; variabletomodel = "imodes"
  psims = aegis::read_write_fast( file.path( p$modeldir, paste("psims", "_", sex, "_", variabletomodel, ".RDS", sep="") ) )

  grw[["f|i|06"]] = psims[["f|i|06"]][, y1, ] - psims[["f|i|05"]][, y0,] 
  grw[["f|i|07"]] = psims[["f|i|07"]][, y1, ] - psims[["f|i|06"]][, y0,] 
  grw[["f|i|08"]] = psims[["f|i|08"]][, y1, ] - psims[["f|i|07"]][, y0,] 
  grw[["f|i|09"]] = psims[["f|i|09"]][, y1, ] - psims[["f|i|08"]][, y0,] 
  grw[["f|i|10"]] = psims[["f|i|10"]][, y1, ] - psims[["f|i|09"]][, y0,] 

  grw[["f|m|09"]] = psims[["f|m|09"]][, y1, ] - psims[["f|i|08"]][, y0,] 
  grw[["f|m|10"]] = psims[["f|m|10"]][, y1, ] - psims[["f|i|09"]][, y0,] 
  grw[["f|m|11"]] = psims[["f|m|11"]][, y1, ] - psims[["f|i|10"]][, y0,] 

  read_write_fast( data=grw, file=file.path( p$modeldir, "growth_increments.RDS")  )
}

grw = aegis::read_write_fast( file.path( p$modeldir, "growth_increments.RDS") )


# map the mean growth increments
grwmeans = list()
for (nm in names(grw)) {
  grwmeans[[nm]] = apply(grw[[ nm ]], MARGIN=c(1,2), mean )
}


for (nm in names(grw)) {
  for (yr in p$yrs[-1]) {
    pg$toplot = grwmeans[[nm]][, as.character(yr)]
    plt = carstm_map( sppoly=pg, vn="toplot",
      title= yr, 
      outfilename=file.path( file.path( p$modeldir, "growth_increments" ), paste(nm, "_", yr, ".png", sep="") ),
      # scale=1.5,
      colors=rev(RColorBrewer::brewer.pal(5, "RdYlBu")),
      additional_features=additional_features
    ) 
    plt
  }
}

# examine relationship of growth increments with location and temperature, depth
# first merge environmental data from SS$M, the data from the prediction surface
# note time for increments is for the terminal time point 

dd = degreedays( pg, years, t0)


M = SS$M[, c("space", "year", "z", "t", "log.substrate.grainsize", "pca1" ), with=FALSE] 
setnames(M, "year", "time" )

M$space = as.character(M$space)
M$time = M$time +  1
M$time = as.factor(as.character(M$time))
 
grw = rbindlist( lapply( grwmeans, function(x) as.data.table( as.data.frame.table(x) )), idcol=TRUE)
setnames(grw, ".id", "id" )

grw = M[grw, on=c("space", "time")]
setnames(grw, "Freq", "log_growth_increment" )
grw$space = as.numeric(grw$space)
grw$space_time = grw$space
grw$time_space = grw$time

grw$log.z = log(grw$z)

plot( log_growth_increment ~ z, grw[id=="m|i|06",]  )

o = inla( log_growth_increment ~ id + f(time, model="iid") +
  f(inla.group( log.z, method="quantile", n=7 ), model="rw2", scale.model=TRUE) + 
  f(inla.group( t, method="quantile", n=7 ), model="rw2", scale.model=TRUE) +
  f(inla.group( log.substrate.grainsize, method="quantile", n=7 ), model="rw2", scale.model=TRUE) +
  f(inla.group( pca1, method="quantile", n=7 ), model="rw2", scale.model=TRUE) +
  f( space, model="bym2", graph=slot(pg, "nb"), scale.model=TRUE ) +
  f( space_time, model="bym2", graph=slot(pg, "nb"), scale.model=TRUE, group=time_space,  control.group=list(model="iid"))  , 
  data=as.data.table(grw ) )

cor( o$summary.fitted.values[,"mean"], grw$log_growth_increment )
  # 0.793 spatial model

```

### Reconstruct size structure at arbitrary areal units

```{r}
#| eval: false 
#| output: false
#| echo: false
#| label: kmm-reconstruct

## Next simulate with julia called from R:
if ( !any( grepl("JuliaCall", o[,"Package"] ) ) ) install.packages("JuliaCall")

# load JuliaCall interface
library(JuliaCall)

julia = try( julia_setup( install=FALSE, installJulia=FALSE ) )

if ( inherits(julia, "try-error") ) {
  install_julia()
  julia = try( julia_setup( install=FALSE, installJulia=FALSE ) ) # make sure it works
  if ( inherits(julia, "try-error") )  stop( "Julia install failed, install manually?") 
}

# currentwd = getwd() 

julia_command( "using Distributions" )  # for kmm_sample (defined in kmm_functions.jl)

include( joinpath( homedir(), "projects", "model", "kmm_functions.jl" ))   

julia_assign( "N", N )  
julia_assign( "n_imodes", n_imodes )  

# setwd(p$output.dir)  # julia_source cannot traverse directories .. temporarily switch directory
out = list()

for (stage in p$stages) {
  out[[stage]] = list()

  M0 = SS$sk[[stage]]

  for (variabletomodel in p$varsnames) {

    p$carstm_model_label = paste("kmm", stage, variabletomodel, sep="_")
    p$variabletomodel = variabletomodel

    res = NULL
    res = carstm_model( p=p, DS="carstm_predictions", sppoly=pg ) # to load currently saved results
    if (is.null(res)) next()

    ## transfer params to julia environment
    sims = list()

    for (au in aus) {
      res_pred = res$predictions[aui,,]       
      julia_assign( "res_pred", res_pred )  # copy data into julia session
      julia_command("out = kmm_sample( res, N=1, n_imodes=1 )")
      sims[[au]] = julia_eval("out")  # copy results back to R
    }

    out[[stage]][[variabletomodel]] = sims
  
  }
}

read_write_fast( data=out, file="outputfile.RDS" )


```
