# automatic load of things related to all projects go here
 
current_directory =  @__DIR__() 
print( "Current directory is: ", current_directory, "\n\n" )

if !@isdefined base_directory 
    @static if Sys.iswindows()
        base_directory = joinpath( "C:\\", "home", "jae", "projects" )
    elseif Sys.islinux()
        base_directory = joinpath(homedir(), "projects" )
    elseif Sys.isapple()
    else
    end
end
 
if !@isdefined base_directory 
    project_directory = joinpath( base_directory, "model_size")
end


import Pkg  # or using Pkg
Pkg.activate(project_directory)  # so now you activate the package
Base.active_project()  
push!( LOAD_PATH, project_directory )  # add the directory to the load path, so it can be found


pkgs = [
    "Revise", 
    "MKL", # "Plots", # "OpenBLAS32", 
    "Logging", "Random", "Setfield",  "ForwardDiff",  
    "CodecZstd", "RData",  "RCall",
    "DataFrames", "JLD2", "CSV",
    "PlotThemes", "Colors", "ColorSchemes", "StatsPlots",  "CairoMakie", "MakieExtra",
    "StatsBase", "Statistics", "Distributions", "KernelDensity", "GLM",
    "MultivariateStats", "StaticArrays", "LazyArrays", "FillArrays",
    "Memoize",
    "Graphs", "Distances", "Peaks", "KernelDensity",
    "Interpolations", "LinearAlgebra", "DynamicHMC", "Turing" 
]
      

# load directly can cause conflicts due to same function names 
pkgtoskipload = [  "RCall",   "CairoMakie", "PlotlyJS",  "PlotlyBase",  "PlotlyKaleido", "LazyArrays" ]
 
print( "Loading libraries:\n\n" ) 
 
for pk in pkgs; 
    if !(Base.find_package(pk) === nothing)
        if !(pk in pkgtoskipload)
            @eval using $(Symbol(pk)); 
        end
    end
end


include( srcdir( "shared_functions.jl") )

print( "\nTo (re)-install required packages, run:  install_required_packages() or Pkg.instantiate() \n\n" ) 

 

# using ApproximateGPs, Random, "CodeTracking",   "Setfield", "ParameterHandling" 
#   "AdvancedHMC", "DynamicHMC", "DistributionsAD", "Bijectors", "Libtask", "ReverseDiff", 
    # "Symbolics", "Logging", 


DEBUG = Ref{Any}()

# add this inside of a function to track vars
# Main.DEBUG[] = y,p,t
