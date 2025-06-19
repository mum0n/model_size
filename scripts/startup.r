 
# startup R environment for model_size 

source( file.path( project_directory, "R", "loadfunctions.r") )

loadfunctions( "aegis")
loadfunctions( "bio.snowcrab")

# main functions
fns = c(
  "sizestructure_db.R",
  "model_size_data_carstm.R",
  "model_size_presence_absence.R",
  "summarize_observations.R",
  "individual_sampling_weights.R"
)

for (fn in fns) source( file.path(project_directory, "R", fn) )


# require(bio.snowcrab)
require(ggplot2)
require(data.table)


p = snowcrab_parameters(
    project_class="carstm",
    yrs=yrs,
    areal_units_type="tesselation",
    carstm_model_label = carstm_model_label,
    dimensionality = "space-time-cyclic",
    selection = selection
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


 
if ( !file.exists(p$project_name) ) dir.create( p$project_name, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$datadir) ) dir.create( p$datadir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$modeldir) ) dir.create( p$modeldir, showWarnings=FALSE, recursive=TRUE )
if ( !file.exists(p$project.outputdir) ) dir.create( p$project.outputdir, showWarnings=FALSE, recursive=TRUE )
