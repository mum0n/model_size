 
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

  # modes of paramaters closer to final solution
  
  out = switch( bioclass,
    f.imm = c(
      9.3291, 10.1394, -0.7173, 1.4523, 0.3921, 0.8456, 1.8795, -0.0039, 0.4538, -0.6189, -0.8713, 2.3288, 2.1746, -0.6796, 5.1048, 43.8626, 21.1540, 33.1340
    ),
    f.mat = c(
      6.0993, 14.3600, 9.3252, 8.8198, 1.9373, 3.0580, 0.6737, 0.9012, 0.7223, -3.4820, -1.9508, 2.4880, 1.3039, 7.4967, 12.6800, 14.5784, 6.0851, 23.0806
    ),
    m.imm = c(
          9.9200, -3.0405, -2.3969, 1.9333, 0.6128, 1.0123, 1.5371, 0.1702, 1.2130, -3.1902, -1.2038, 3.6389, 3.3960, -2.3632, 3.9804, -1.6675, -0.7746, 0.8269),
    m.mat = c(
      9.903, 0.354, 2.297, 1.557, 0.825, 2.174, 1.416, 0.00429, 2.074, -1.795, 2.818, 4.046, 4.5681, -2.496, 4.354, -1.251, -0.343, 1.025
    ),
    NULL
  )
  return(out) 
}


  
