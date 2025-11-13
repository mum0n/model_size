

model_size_presence_absence = function( p, todo="load", 
  num.threads="1:1:1", improve.fit=FALSE, theta0 = NULL) {
     
    p$selection$biologicals_using_snowcrab_filter_class = p$bioclass
   
    fn = file.path(  
        p$modeldir, 
        p$bioclass,
        paste("fit_", p$bioclass, ".rdz", sep="") 
    )

    fit = NULL
    if ( "load" %in% todo ) {
        if (file.exists(fn)){
            message( "Loading model fit: ", fn)
            fit = read_write_fast(fn)
            if (!is.null(fit)) return(fit)
        }
        message("\n\nModel fit not found ... regenerating")
    }
    
    message( "\n---\n---\n---\n---\n")
    
    message("Creating model fit: ", fn, "\n")
         
    M = model_size_data_carstm( p=p )  
    
    # summary on link scale
    MS = summarize_observations(
        obs = M[tag=="observations", pa ],
        offsets = M[tag=="observations", data_offset ], 
        family="binomial"
    )  
    
    ## these must be in global namespace for INLA
    
    sppoly <<- attributes(M)$sppoly

    H <<- inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] ) 

    if (is.null(theta0)) {
      if (exists( "theta", p )) theta0 = p[["theta"]][[ p$bioclass ]]
    }

    message("Trying default 'compact_tauc' inla mode: ")
    message( "\n---\n---\n" )

    fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        verbose=TRUE,
        inla.mode="compact",
        control.inla = list( strategy = "gaussian", int.strategy="eb" ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    ), silent=TRUE )

    message("Trying default 'classic' inla mode: ")
    message( "\n---\n---\n" )

    fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        inla.mode="classic",
        control.inla = list( strategy = "gaussian", int.strategy="eb" ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, config=TRUE, return.marginals.predictor=TRUE  ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    ), silent=TRUE )


    if (inherits(fit, "try-error" )) {
      message( "\n---\n---\n" )
      message("Trying default 'compact' inla mode with more robust settings: ")
      message( "\n---\n---\n" )
  
      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        inla.mode="compact",
        control.inla = list( strategy = "gaussian", fast=TRUE, h=0.01, int.strategy="eb", force.diagonal=TRUE,  cmin=0.0001 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE, save.memory=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
      ), silent=TRUE )
   }

   if (inherits(fit, "try-error" )) {
      message( "\n---\n---\n" )
      message( "Trying default 'compact' inla mode with using starting modes: ")
      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        inla.mode="compact",
        control.inla = list( strategy = "auto", fast=FALSE, int.strategy="auto" ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE, smtp="tauc" ),
        num.threads=num.threads
      ), silent=TRUE )
   }



    if (inherits(fit, "try-error" )) {
      message( "\n---\n---\n" )
      message( "Trying the 'classic' inla mode, more stable and robust but RAM intensive: ")
      message( "\n---\n---\n" )

      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        debug = TRUE,
        inla.mode="classic",
        control.inla = list( strategy="laplace", h=0.01, force.diagonal=TRUE, cmin=0.001 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
      ), silent=TRUE )

    }
 
    if (inherits(fit, "try-error" )) {
      message( "\n---\n---\n" )
      message( "Model fit error! Adjust options or check data. Giving up ...")
      message( "\n---\n---\n" )
      return(fit)
    }
 
 
    if (improve.fit) {
        message( "\n---\n---\n" )
        message( "Improving model fit: ", fn,)
        message( "\n---\n---\n" )

        fit = inla.hyperpar( fit, verbose=TRUE, restart=TRUE )  #  dz = 0.75, diff.logdens = 15,
    }

    message( "\n---\n---\n" )
    message( "Model theta estimate (use as modes for restarting):")
    message( paste0( round(fit$mode$theta, 4), collapse=", ") )
    message( "\n---\n---\n" )
    
    fit$modelinfo = list(
        bioclass = p$bioclass, 
        H = H,
        MS = MS,
        theta0 = theta0,
        fn = fn
    )
    
    message( "\n---\n---\n" )
    message( "Saving model fit: ", fn, "\n")
    message( "\n---\n---\n" )

    read_write_fast( data=fit, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)
}

