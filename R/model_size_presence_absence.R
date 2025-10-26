

model_size_presence_absence = function( p, todo="load", num.threads="2:1", improve.fit=FALSE) {
     
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
    }
    
    message( "\n\n\n---\n---\nCreating model fit: ", fn, "\n---\n---\n\n\n")
         
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


    theta0 = p$theta( p$bioclass )
    
    message("\n---\n---\nTrying default compact inla mode: \n---\n---\n")
    fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    ), silent=TRUE )


    if (inherits(fit, "try-error" )) {
      message("\n---\n---\nTrying default compact inla mode with more robust settings: \n---\n---\n")
      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        control.inla = list( strategy = "gaussian", control.vb = list(enable = TRUE), reordering="metis", diagonal=1e6, cmin=0.0001 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
      ), silent=TRUE )
   }

   if (inherits(fit, "try-error" )) {
      message("\n---\n---\nTrying default compact inla mode with generic starting modes: \n---\n---\n")
      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        num.threads=num.threads
      ), silent=TRUE )
   }



    if (inherits(fit, "try-error" )) {
      message("\n---\n---\nTrying the more stable 'classic' inla mode: \n---\n---\n")
      fit = try( inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        safe = FALSE,
        verbose=TRUE,
        inla.mode="classic",
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.inla = list( strategy="simplified.laplace", control.vb = list(enable = TRUE), reordering="metis", diagonal=1e5, cmin=0.0 ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
      ), silent=TRUE )

    }
 
    if (inherits(fit, "try-error" )) {
        message("\n---\n---\n\nModel fit error! \n\n---\n---\n")
    }
 
 
    if (improve.fit) {
        message( "\n---\n---\nImproving model fit: ", fn, "\n---\n---\n")
        fit = inla.hyperpar( fit, verbose=TRUE, restart=TRUE )  #  dz = 0.75, diff.logdens = 15,
    }

    message( "Model theta (modes) estimate:  \n---\n---\n", paste0( round(fit$mode$theta, 4), collapse=", "), "\n---\n---\n" )
    
    fit$modelinfo = list(
        bioclass = p$bioclass, 
        H = H,
        MS = MS,
        theta0 = theta0,
        fn = fn
    )
    
    message( "\n\n\n---\n---\nSaving model fit: ", fn, "\n", "\n---\n---\n\n\n" )
    read_write_fast( data=fit, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)
}

