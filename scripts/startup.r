 
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
9.8791, -2.1682, -1.3834, 1.3522, 0.7025, 0.7836, 2.2636,
4e-04, -0.6663, -0.5588, 0.1819, 3.201, 2.1819, -2.2789, 3.7194,
-1.6712, -0.1614, 0.7537
    ),
    f.mat = c(
9.5326, -0.7771, 1.3083, 1.8652, -0.6052, 1.5363, 0.0204,
0.1819, 0.2174, -3.1332, -2.33, 2.0505, 1.0892, -3.8373, 12.1444,
-2.891, -0.6585, 1.3005
    ),
    m.imm = c(
9.9177, -2.5095, -2.3357, 2.1833, 1.4904, 1.3417, 1.5166,
0.0069, -0.4652, -0.2505, -1.338, 4.2284, 2.0433, -2.2634, 4.336,
-1.3908, -0.7029, 0.8225
    ),
    m.mat = c(
9.9007, 0.2020, 1.4624, 2.1964, 0.5940, 2.5761, 1.9452,
0.0102, 0.2808, 1.2686, -0.1217, 3.2817, 4.2241, -2.4502, 4.4544,
-1.1382, -0.3519, 1.0917
    ) 
)
   

p$xrange = c(8, 170)  # size range (CW)

figures_dir = file.path( project_directory, "outputs", "figures" )
dx = 2 #  width of carapace with discretization to produce "cwd"
ss_outdir = file.path(p$project.outputdir, "size_structure")
sizedatadir = file.path( homedir, "bio.data", "bio.snowcrab", "output", "size_structure" )


