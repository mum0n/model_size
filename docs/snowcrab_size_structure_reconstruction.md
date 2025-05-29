
# Reconstruct size structure

Basic premise:

The probability of observing a snow crab $P_\text{obs}$ is assumed to be composed of independent factors:

$$P_\text{obs} = P_\text{size} \cdot P_\text{sex} \cdot P_\text{maturity} \cdot P_\text{location} \cdot P_\text{time} \cdot P_\text{depth} \cdot P_\text{temperature}$$

These probabilities can be decomposed from a Binomial Bernoulli model of presence-absence of observations by sampling event (survey trawl set).

To reduce the number of possible combinations of the above factors, we can use size bins that span a range and then use a rw2 model to smooth across the sizes.

## Analysis to determine sample weight of each individual observation

The main idea is that size-structure cannot be examined easily as each individual measurement is derived from a sample of unequal sampling intensity (swept are), size-bias of sampling gear and environmental bias conditions. This does not stop people from still aggregating and so ignoring these differences in sampling weights. Accounting for size-based bias is attempted when this bias is well understood, though not always.

Here we attempt to factor in all such biases into one model.


## Load environment

```{r}
#| eval: true
#| output: false
#| echo: false
#| label: setup-R-environment

# homedir="C:/home/jae/"  # on MSwindows

project_name = "model_size"
project_directory = file.path( homedir, "projects", project_name )

source( file.path( project_directory, "R", "loadfunctions.r") )

loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")

# require(bio.snowcrab)
require(ggplot2)
require(data.table)

year_start = 1999
year_assessment = 2024

yrs = year_start:year_assessment

sexid="male"
sexid="female"

p = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,
    areal_units_type="tesselation",
    carstm_model_label= sexid,
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ) ),
      biologicals_using_snowcrab_filter_class=sexid
    )
)

sppoly = areal_units( p=p )


if (0) {
  # direct access:
  aufn = project.datadirectory( "bio.snowcrab", "areal_units", "snowcrab~tesselation~1~snowcrab~8~1~none~snowcrab_managementareas.rdata")
  load(aufn)
}

p$ny = length(p$yrs)
p$nw = 12
p$nt = p$ny * p$nw 
p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components 

p = temporal_parameters(p=p, dimensionality="space-time-cyclic", timezone="America/Halifax")

# override prediction time values (decimal-year), # output timeslices for predictions in decimla years, yes all of them here
tout = expand.grid( yr=p$yrs, dyear=1:p$nw, KEEP.OUT.ATTRS=FALSE )
p$prediction_ts = sort( tout$yr + tout$dyear/p$nw - p$tres/2 )# mid-points


p$space_name = sppoly$AUID
p$space_id = 1:nrow(sppoly)  # must match M$space

p$time_name = as.character(p$yrs)
p$time_id =  1:p$ny

p$cyclic_name = as.character(p$cyclic_levels)
p$cyclic_id = 1:p$nw


# over-ride defaults to get sabove to continue and get size data
p = parameters_add(p, list(
  project_name = project_name,
  project_directory = project_directory,
  data_root = project_directory,
  datadir = project_directory,   # all unprocessed inputs (and simple manipulations) ..
  project.outputdir = file.path( project_directory, "outputs" ), #interpolations and mapping
  modeldir = file.path( project_directory, "outputs", "modelled" )  # all model outputs
) )

source( file.path(project_directory, "R", "sizestructure_db.R") )
source( file.path(project_directory, "R", "model_size_data_carstm.R" ))

xrange = c(5, 170)  # size range (CW)
dx = 5 #  width of carapace with discretization to produce "cwd"

redo_datafiles = FALSE
# redo_datafiles = TRUE


M = model_size_data_carstm(p=p, sppoly=sppoly, xrange=xrange, dx=dx, sexid=sexid,  redo=redo_datafiles )


# M = M[ yr>2010, ]  # for testing

plot(jitter(pa) ~ cw, M, pch=".") # zeros extend beyond to give "prior" info to upper size ranges  (male, female)
plot(jitter(pa) ~ year, M, pch="." )

    # males -- hvariance compression: 2002
    # females -- variance compression: 2000:2003, 2012, 2013, 2018 (low abundance periods)

plot(jitter(pa) ~ dyear, M, pch="." ) # no season-bias (male, female)

plot(dyear ~ year, M, pch="." )  # time-bias up to 1999:2004 (male, female)
plot(t ~ dyear, M, pch="." )  # seasonal temperature bias (male, female)
plot(z ~ dyear, M, pch="." )  # seasonal depth bias (male, female) --shallows in winter Dec-Jan

# observed presence is spanned by observed absence (ie. safely goes beyond distribution ) (male, female)
plot(jitter(pa) ~ t, M, pch=".")
plot(jitter(pa) ~ z, M, pch=".")  # shallow areas sampled in winter (weather)


# , group=mat_group, control.group=list(model="iid" )

# ' + f(inla.group( cw, method="quantile", n=15 ), model="rw2", scale.model=TRUE, group=mat_group, control.group=list(model="iid" ) ) ',
     

p$formula = as.formula( paste(
  ' pa ~ 1 ',
      ' + offset( data_offset ) + mat ',
      ' + f(inla.group( cw, method="quantile", n=11 ), model="rw2", scale.model=TRUE) ',
      ' + f( time, model="ar1",  hyper=H$ar1 ) ',
      # ' + f( cyclic, model="ar1", hyper=H$ar1 )',
      # ' + f( cyclic, model="seasonal", scale.model=TRUE, season.length=12, hyper=H$iid  )',
      ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2)',
      #' + f( inla.group( substrate.grainsize, method="quantile", n=9), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      # ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
      # ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2",  scale.model=TRUE, hyper=H$rw2) ',
      # ' + f( inla.group( pca3, method="quantile", n=9 ), model="rw2",  scale.model=TRUE, hyper=H$rw2) ',
      ' + f( space, model="bym2", graph=slot(sppoly, "nb"),  scale.model=TRUE, hyper=H$bym2 ) ',
      ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"),  scale.model=TRUE, group=time_space, hyper=H$bym2,  control.group=list(model="ar1", hyper=H$ar1_group)) '
) )



# model pa using all data
carstm_model( p=p, data=M, sppoly=sppoly,
  # theta = c( 0.8917, 2.0052, 4.5021, -0.0000, 1.5400, -2.4689, 1.1762, 2.6536, 2.9546, -2.1406, 3.5352, -0.7465, 3.2443, 2.4420 ),
  nposteriors=1,
  redo_fit=TRUE,
  toget = c("summary", "random_spatial", "predictions"),
  posterior_simulations_to_retain = c("summary"),
  family = "binomial",  # "binomial",  # "nbinomial", "betabinomial", "zeroinflatedbinomial0" , "zeroinflatednbinomial0"
  # redo_fit=FALSE,
  # debug = "summary",
  debug = "predictions",
  # control.family=list(control.link=list(model="logit")),  # default for binomial .. no need to specify
  control.inla = list( strategy="adaptive", int.strategy="eb", h=0.05 ),
  # num.threads="4:3",
  verbose=TRUE
)


#  p$carstm_directory = file.path(p$modeldir, p$carstm_model_label)

  # modelinfo = carstm_model(p = p, DS = "carstm_modelinfo")
  fit = carstm_model(p = p, DS = "modelled_fit")
  vn = "pa"

  iobs = which(M$tag == "observations")
  ipreds = which(M$tag == "predictions")

  O = M[iobs, ]
  O$fitted_mean = fit$summary.fitted.values[["mean"]][iobs]
  O$fitted_sd = fit$summary.fitted.values[["sd"]][iobs]

  cor(O$pa, O$fitted_mean)


  P = M[ipreds, ]
  P$prediction_mean = fit$summary.fitted.values[["mean"]][ipreds]
  P$prediction_sd = fit$summary.fitted.values[["mean"]][ipreds]


 #  fit = M = NULL; gc()

P = P[, .(AUID, year, cyclic, cw, mat, prediction_mean, prediction_sd)]

O = O[ pa==1,]


breaks = seq(xrange[1], xrange[2], by=dx)
mids = breaks[-length(breaks)] + dx/2

O$cw = discretize_data( O$cw, brks=breaks, labels=mids, resolution=dx )

O = P[ O, on=.(AUID, year, cyclic, cw, mat)]

O$relative_rate = O$fitted_mean / O$prediction_mean



model_summary = carstm_model(p = p, DS = "carstm_summary")


S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores )

# posterior predictive check
carstm_posterior_predictive_check(p=p, M=M[ , ]  )

# EXAMINE POSTERIORS AND PRIORS
res = carstm_model(  p=p, DS="carstm_summary" )  # parameters in p and summary

outputdir = file.path(p$modeldir, p$carstm_model_label)

res_vars = c( names( res$hypers), names(res$fixed) )
for (i in 1:length(res_vars) ) {
  o = carstm_prior_posterior_compare( res, vn=res_vars[i], outputdir=outputdir  )
  dev.new(); print(o)
}


plot( jitter(M$pa), fit$summary.fitted.values$mean, pch="." )
cor( jitter(M$pa), fit$summary.fitted.values$mean ) # 0.5441


vn = 'inla.group(cw, method = "quantile", n = 15)'
o = fit$summary.random[[vn]]
plot( o$ID, (o$mean) )

vn = 'inla.group(z, method = "quantile", n = 9)'
o = fit$summary.random[[vn]]
plot( o$ID, (o$mean) )

vn = 'inla.group(t, method = "quantile", n = 9)'
o = fit$summary.random[[vn]]
plot( o$ID, (o$mean) )

o = fit$summary.fixed
plot( (o$mean) )

 
# for mapping below, some bathymetry and polygons
additional_features = snowcrab_mapping_features(p)

  ylab = "Probability"
  fn_root_prefix = "Predicted_presence_absence"
  fn_root = "habitat"
  # title= paste( snowcrab_filter_class, "Probability")

outputdir = file.path( p$modeldir, p$carstm_model_label, "figures" )
if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )


    #   # to compute habitat prob
    #   sims = carstm_posterior_simulations( pH=p, pa_threshold=0.05, qmax=0.95 )
    #   SM = aggregate_simulations(
    #     sims=sims,
    #     sppoly=sppoly,
    #     fn=carstm_filenames( p, returnvalue="filename", fn="aggregated_timeseries" ),
    #     yrs=p$yrs,
    #     method="mean",
    #     redo=TRUE
    #   )
    #   outputdir = file.path( carstm_filenames( p, returnvalue="output_directory"), "aggregated_habitat_timeseries" )
    #   ylabel = "Habitat probability"
    #   fn_ts = "habitat_M0.png"
    #   vn = paste("habitat", "predicted", sep=".")
    #   outputdir2 = file.path( carstm_filenames( p, returnvalue="output_directory"), "predicted_habitat" )




# to load currently saved results
res = carstm_model( p=p,  DS="carstm_summary" )  # parameters in p and direct summary
res$direct
res = NULL; gc()


# plots with 95% PI
carstm_plot_marginaleffects( p, outputdir, fn_root )


# maps of some of the results
outputdirmap = file.path(p$modeldir, p$carstm_model_label, "maps" )

carstm_plot_map( p=p, outputdir=outputdirmap, fn_root_prefix=fn_root_prefix , additional_features=additional_features,
  toplot="random_spatial", probs=c(0.025, 0.975) )

carstm_plot_map( p=p, outputdir=outputdirmap, fn_root_prefix=fn_root_prefix , additional_features=additional_features,
  toplot="predictions", probs=c(0.1, 0.9))

}



```




## Analysis of size structure .. tests

```{r}

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


# required

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


# try an individual-based model first:





o = size_distributions(p=p, toget="tabulated_data_by_stage", add_zeros=TRUE )

o = o[-.(region, year)]

o = size_distributions(p=p, toget="tabulated_data_by_stage" )


o = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE )







```


