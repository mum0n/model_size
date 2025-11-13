 
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
    ' + f( inla.group( cwd, method="quantile", n=9 ), model="rw2", scale.model=TRUE) ', 
    ' + f( inla.group( cwd2, method="quantile", n=9 ), model="rw2", scale.model=TRUE, group=space2, control.group=list(model="besag", hyper=H$besag, graph=slot(sppoly, "nb") ) ) ',    
    ' + f( inla.group( cwd3, method="quantile", n=9 ), model="rw2", scale.model=TRUE, group=time3, control.group=list(model="ar1", hyper=H$ar1_group) ) ',   
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

# start modes of paramaters closer to known solutions
p$theta = list(
    f.imm = c(
7.1393, 10.3944, -0.9175, 1.2138, -1.9280, 0.2750, 1.6001, -1.4661, 0.1789, -1.1114, 1.4618, 2.4874, 2.8718, 0.5663, 4.9777, 43.4759, 11.1819, 31.3289
    ),
    f.mat = c(
8.8502, -0.7965, 1.3387, 1.8563, -0.6780, 1.5799, 0.0386, 0.1632, 0.6436, -2.9773, -2.3924, 2.2028, 1.0249, -3.8357, 12.1025, -2.8908, -0.6615, 1.2978
# 7.2095, -0.8070, 1.3985, 1.7830, -0.5914, 1.3498, 0.0768, 0.6744, -1.6007, -2.6736, -2.7108, 2.8759, 1.2127, -3.8437, 12.2346, -2.8731, -0.7006, 1.2860
    ),
    m.imm = c(
9.9185, -2.5078, -2.341, 2.1849, 1.5128, 1.3547, 1.502, 0.0052, -0.4524, -0.2573, -1.3309, 4.2263, 2.0461, -2.2527, 4.3384, -1.3941, -0.6868, 0.8241
    ),
    m.mat = c(
9.9007, 0.1963, 1.4610, 2.1988, 0.6028, 2.5789, 1.9486, 0.0005, 0.2798, 1.2776, -0.1185, 3.2845, 4.2272, -2.4450, 4.4546, -1.1401, -0.3474, 1.0919
    ) 
)

 
p$xrange = c(8, 170)  # size range (CW)

figures_dir = file.path( project_directory, "outputs", "figures" )
dx = 2 #  width of carapace with discretization to produce "cwd"
ss_outdir = file.path(p$project.outputdir, "size_structure")
sizedatadir = file.path( homedir, "bio.data", "bio.snowcrab", "output", "size_structure" )


