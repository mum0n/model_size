
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

year_start = 1999
year_assessment = 2024
yrs = year_start:year_assessment

# choose one:
sexid="male"
sexid="female"

carstm_model_label = sexid

selection = list(
  type = "presence_absence",
  biologicals=list( spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ) ),
  biologicals_using_snowcrab_filter_class=sexid
)

source( file.path( project_directory, "scripts", "startup.r") )


sppoly = areal_units( p=p )

p$ny = length(p$yrs)
p$nt = p$ny # must specify, else assumed = 1 (1= no time)  ## nt=ny annual time steps, nt = ny*nw is seassonal
p$nw = 12 
p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
p$dyears = discretize_data( span=c(0, 1, p$nw), toreturn="lower" )  # left breaks .. get 12 intervals of decimal years... fractional year breaks

p$dyear_centre = discretize_data( span=c(0, 1, p$nw), toreturn="midpoints" ) 

# used for creating timeslices and predictions  .. needs to match the values in aegis_parameters()
p$prediction_dyear = lubridate::decimal_date( lubridate::ymd("0000/Sep/01"))
# output timeslices for predictions in decimla years, yes all of them here
p$prediction_ts = p$yrs + p$prediction_dyear

p = temporal_parameters(p=p, dimensionality="space-time-cyclic", timezone="America/Halifax")

p$space_name = sppoly$AUID 
p$space_id = 1:length(p$space_name)  # must match M$space

p$time_name = as.character(p$yrs)
p$time_id =  1:p$ny

cyclic_levels = discretize_data( brks=p$dyears, toreturn="midpoints" )   # default midpoints; same as:
p$cyclic_name = as.character(p$cyclic_levels)
p$cyclic_id = 1:p$nw


# note ranges in CW will be log transformed later
p$span = switch(sexid,
    male   = c( 5, 155, 40),
    female = c( 5, 95,  40)
)

p$formula = as.formula( paste(
' pa ~ 1 ',
    ' + offset( data_offset ) + mat ', 
    ' + f( inla.group( cwd, method="quantile", n=13 ), model="rw2", scale.model=TRUE) ', 
    ' + f( inla.group( cwd2, method="quantile", n=13 ), model="rw2", scale.model=TRUE, group=mat_group, hyper=H$rw2, control.group=list(model="iid", hyper=H$iid)) ', 
    ' + f( time, model="ar1",  hyper=H$ar1 ) ',
    ' + f( cyclic, model="ar1", hyper=H$ar1 )',
    # ' + f( cyclic, model="seasonal", scale.model=TRUE, season.length=12, hyper=H$iid  )',
    ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    # ' + f( inla.group( pca3, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2 ) ',
    ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
    # ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="iid", hyper=H$iid)) '
) )

 

# individual level data from snow crab surveys and prediction surface
# use a log scale for size info .. full range .. it will be reduced to data range in the next step
# male: 30+ bins seem sufficient to describe the shape
# female: 30+ bins seem sufficient to describe the shape

redo_datafiles = FALSE
# redo_datafiles = TRUE

M = model_size_data_carstm(p=p, sexid=sexid, sppoly=sppoly, redo=redo_datafiles, carstm_set_redo=redo_datafiles )  


hist((as.numeric(as.character(M$cwd[M$pa==1]))), "fd")  


plot(jitter(pa) ~ cwd, M, pch=".") # zeros extend beyond to give "prior" info to upper size ranges  (male, female)
plot(jitter(pa) ~ year, M, pch="." )  
    # males -- variance compression: 2002  
    # females -- variance compression: 2000:2003, 2012, 2013, 2018 (low abundance periods)

plot(jitter(pa) ~ dyear, M, pch="." ) # no season-bias (male, female)

plot(dyear ~ year, M, pch="." )  # time-bias up to 1999:2004 (male, female)
plot(t ~ dyear, M, pch="." )  # seasonal temperature bias (male, female)
plot(z ~ dyear, M, pch="." )  # seasonal depth bias (male, female) --shallows in winter Dec-Jan

# observed presence is spanned by observed absence (ie. safely goes beyond distribution ) (male, female)
plot(jitter(pa) ~ t, M, pch=".")
plot(jitter(pa) ~ z, M, pch=".")  # shallow areas sampled in winter (weather)
  

# action = "redo"
action = "load"

fit = model_size_presence_absence( p=p, sexid=sexid, sppoly=sppoly, action=action ) 
summary(fit)

fit = NULL; gc()


# action = "redo"
action = "load"

O = individual_sampling_weights( p=p, sexid=sexid, sppoly=sppoly, action=action ) 



#  p$carstm_directory = file.path(p$modeldir, p$carstm_model_label)

#  modelinfo = carstm_model(p = p, DS = "carstm_modelinfo") 
#  fit = carstm_model(p = p, DS = "modelled_fit")

cor(O$pa, O$fitted_mean)  

    if (0){ 

        male 11 and no group mat 
        cor 0.7329

        Deviance Information Criterion (DIC) ...............: 605046.98
        Deviance Information Criterion (DIC, saturated) ....: 603874.13
        Effective number of parameters .....................: 9584.17

        Watanabe-Akaike information criterion (WAIC) ...: 606534.62
        Effective number of parameters .................: 10811.25

        Marginal log-Likelihood:  -302117.80


        male 13 and group mat:

        cor=0.7804

        Deviance Information Criterion (DIC) ...............: 532406.39
        Deviance Information Criterion (DIC, saturated) ....: 531233.54
        Effective number of parameters .....................: 10986.52

        Watanabe-Akaike information criterion (WAIC) ...: 567799.45
        Effective number of parameters .................: 28804.10

        Marginal log-Likelihood:  -264656.11


        female 13 and group mat

        cor = 0.8136

        Deviance Information Criterion (DIC) ...............: 350406.27
        Deviance Information Criterion (DIC, saturated) ....: 349374.65
        Effective number of parameters .....................: 7448.05

        Watanabe-Akaike information criterion (WAIC) ...: 350224.72
        Effective number of parameters .................: 7056.97

        Marginal log-Likelihood:  -175462.54
         

    }


  hist(O$relative_rate, "fd")
  summary(O$relative_rate)
  plot( O$fitted_mean, O$prediction_mean, pch="." )
  
  plot( exp(O$fitted_mean), exp(O$prediction_mean), pch=".")
  
  cor( (O$fitted_mean), (O$prediction_mean), use="pairwise.complete.obs") 
  cor( exp(O$fitted_mean), exp(O$prediction_mean), use="pairwise.complete.obs") # 0.732 M, -0.00157 F

 
  O$w = O$relative_rate * O$cfanorth_surfacearea
  yrp = "2019"
  x11(); ggplot(O[yr==yrp & w>0,], aes( log(cw), weight=w ) ) + geom_histogram( bins = span[3] )
  x11(); ggplot(O[yr==yrp & w>0,], aes( log(cw) ) ) + geom_histogram( bins = span[3] )


  O$w = O$relative_rate * O$au_sa_km2 
  x11(); ggplot(O[w>0,], aes( log(cw), weight=w ) ) + geom_histogram( bins = span[3] )
  x11(); ggplot(O[w>0,], aes( log(cw) ) ) + geom_histogram( bins = span[3] )



# ----

# snippets and ideas


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
 


## to add spatial effects:
 yrs = 1999:2024

  spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 )
  
  snowcrab_filter_class = "fb"     # fishable biomass (including soft-shelled )  "m.mat" "f.mat" "imm"
   
  carstm_model_label= paste( "male_size_structure" )
  carstm_model_label= paste( "female_size_structure" )
 
 
 
    io = which(M$tag=="observations")
    ip = which(M$tag=="predictions")
    iq = unique( c( which( M$totno > 0), ip ) )
    iw = unique( c( which( M$totno > 5), ip ) )  # need a good sample to estimate mean size
 
    carstm_model( p=p, data=M[ , ], sppoly=sppoly, 
      # theta=c( 2.7291,1.8146,2.9382,0.0132,3.8666,-0.0211,4.2673,5.5037,6.1421,0.2391,4.2522,0.7666,-0.0100,0.8763 ),
      nposteriors=5000,
      toget = c("summary", "random_spatial", "predictions"),
      posterior_simulations_to_retain = c("predictions"),
      family = "poisson",
      verbose=TRUE,
      # redo_fit=FALSE, 
      # debug = "summary",
      # debug = "predictions",
      num.threads="4:3"  
    ) 

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

 

  # for mapping below, some bathymetry and polygons
  additional_features = snowcrab_mapping_features(p)  
  
        ylab = "Probability"
        fn_root_prefix = "Predicted_presence_absence"
        fn_root = "habitat"
        # title= paste( snowcrab_filter_class, "Probability")  
     
    outputdir = file.path( p$modeldir, p$carstm_model_label, "figures" )
      if ( !file.exists(outputdir)) dir.create( outputdir, recursive=TRUE, showWarnings=FALSE )

# plots with 95% PI
carstm_plot_marginaleffects( p, outputdir, fn_root )

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
 
      

# maps of some of the results
outputdirmap = file.path(p$modeldir, p$carstm_model_label, "maps" )

carstm_plot_map( p=p, outputdir=outputdirmap, fn_root_prefix=fn_root_prefix , additional_features=additional_features,
  toplot="random_spatial", probs=c(0.025, 0.975) )

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



# try an individual-based model first:





o = size_distributions(p=p, toget="tabulated_data_by_stage", add_zeros=TRUE )

o = o[-.(region, year)]

o = size_distributions(p=p, toget="tabulated_data_by_stage" )


o = size_distributions(p=p, toget="tabulated_data", add_zeros=TRUE )







```


