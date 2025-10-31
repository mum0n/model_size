 
# startup R environment for model_size 


# main functions
fns = list.files( path=project_directory, pattern="*.R$", full.names=T, recursive=T, ignore.case=T, include.dirs=F )
fns = fns[-which( grepl( "startup.r", fns )) ]
for (fn in fns) source( fn)

# source( file.path( project_directory, "R", "loadfunctions.r") )
loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")


# require(bio.snowcrab)
require(ggplot2)
require(data.table)


p = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,
    areal_units_type="tesselation",
    carstm_model_label = "model_size",
    dimensionality = "space-time-cyclic",
    selection = list(
      type = "presence_absence",
      biologicals=list( spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ) )
    )
)


# over-ride defaults to get above to continue and get size data
p = parameters_add(p, list(
  project_name = project_name,
  project_directory = project_directory,
  data_root = project_directory,
  datadir = project_directory,   # all unprocessed inputs (and simple manipulations) ..    
  project.outputdir = file.path( project_directory, "outputs" ), # interpolations and mapping
  modeldir = file.path( project_directory, "outputs", "modelled" )  # all model outputs
) )

 
if ( !file.exists(p$project_directory) ) dir.create( p$project_directory, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$project.outputdir) ) dir.create( p$project.outputdir, showWarnings=FALSE, recursive=TRUE )



# override parameter defaults: (predict on every time slice)

p$ny = length(p$yrs)
p$nw = 12 
p$nt = p$ny * p$nw  # must specify, else assumed = 1 (1= no time)  ## nt=ny annual time steps, nt = ny*nw is seassonal

p$tres = 1/ p$nw # time resolution .. predictions are made with models that use seasonal components
p$dyears = discretize_data( span=c(0, 1, p$nw), toreturn="lower" )  # left breaks .. get 12 intervals of decimal years... fractional year breaks

# used for creating timeslices and predictions  .. needs to match the values in aegis_parameters()
# output timeslices for predictions in decimla years, yes all of them here
dyears_midpoints = discretize_data( span=c(0, 1, p$nw), toreturn="midpoints" ) 

# output timeslices for predictions in decimal years, yes all of them here
tout = expand.grid( yr=p$yrs, dyear=dyears_midpoints , KEEP.OUT.ATTRS=FALSE )
p$prediction_ts = sort( tout$yr + tout$dyear  ) # mid-points

# a time-slice to predict upon in order to remove seasonal effects
p$prediction_dyear_index = 9  # (month) september

p = temporal_parameters(p=p, dimensionality="space-time-cyclic", timezone="America/Halifax")
 
# define size range in CW .. log transformed internally
p$size_range = function( bioclass ) {
  out = switch(bioclass,
    all    = c( 5,  155),
    male   = c( 5,  155),
    female = c( 5,   95),
    m.imm  = c( 5,  135),
    m.mat  = c( 50, 155),
    f.imm  = c( 5,   80),
    f.mat  = c( 35,  95)
  ) 
  return(out)
}

p$size_bandwidth = 30  # no of size bins (for modelling size-bias only) .. on a log scale internally


# model formula: 

p$formula = as.formula( paste(
' pa ~ 1 ',
    ' + offset( data_offset ) ', 
    ' + f( inla.group( cwd, method="quantile", n=11 ), model="rw2", scale.model=TRUE) ', 
    ' + f( inla.group( cwd2, method="quantile", n=11 ), model="rw2", scale.model=TRUE, group=space2, control.group=list(model="besag", hyper=H$besag, graph=slot(sppoly, "nb") ) ) ',    
    ' + f( inla.group( cwd3, method="quantile", n=11 ), model="rw2", scale.model=TRUE, group=time3, control.group=list(model="ar1", hyper=H$ar1_group) ) ',   
    ' + f( time, model="ar1",  hyper=H$ar1 ) ',
    ' + f( cyclic, model="ar1", hyper=H$ar1 )',
    ' + f( inla.group( t, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( z, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( substrate.grainsize, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( pca1, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( inla.group( pca2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, hyper=H$rw2) ',
    ' + f( space, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, hyper=H$bym2 ) ',
    ' + f( space_time, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, group=time_space, hyper=H$bym2, control.group=list(model="ar1", hyper=H$ar1_group)) '
) )

# alternative:
#    ' + f( space2, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, 
#    group=inla.group( cwd2, method="quantile", n=9 ), hyper=H$besag, control.group=list(model="rw2", hyper=H$rw2)) ', 

p$theta = function(bioclass) {

  # start modes of paramaters closer to known solutions
  
  out = switch( bioclass,
    f.imm = c(
      9.3291, 10.1394, -0.7173, 1.4523, 0.3921, 0.8456, 1.8795, -0.0039, 0.4538, -0.6189, 
      -0.8713, 2.3288, 2.1746, -0.6796, 5.1048, 43.8626, 21.1540, 33.1340
    ),
    f.mat = c(
      6.4364, 14.0605, 11.3605, 7.8415, 4.7670, 4.2591, 4.6586, 2.0921, 4.3206, 1.8505, 
      2.2031, 6.1731, 4.9362, 9.0836, 12.1969, 15.9643, 6.4769, 25.2386
    ),
    m.imm = c(
      9.9191, -2.4676, -1.8416, 1.7711, 1.1651, 0.9926, 1.7656, 0.0257, 0.5771, -0.1422, 
      -1.4755, 3.9797, 1.8674, -2.2559, 4.1655, -1.3974, -0.7044, 0.8213 
    ),
    m.mat = c(
      9.8791, 0.1994, 1.5117, 2.1416, 0.5552, 2.5428, 1.7774, 0.015, 1.9001, -1.9128, 0.2067, 
      3.7526, 4.9292, -2.4348, 4.4517, -1.1324, -0.3161, 1.0651
    ),
    NULL
  )
  return(out) 
}

 