

model_size_presence_absence = function( p, todo=c("fit", "fit_preds"), 
  num.threads="1:1:1", improve.fit=FALSE, theta0 = NULL, restart=TRUE,
  inla.mode="compact", verbose=TRUE,
  control.inla = list( strategy="gaussian", int.strategy="eb",  fast=TRUE, h=0.01, force.diagonal=TRUE,  cmin=0.0001 ),
  control.compute = list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE, save.memory=TRUE )
  ) {
     
    p$selection$biologicals_using_snowcrab_filter_class = p$bioclass
   
    fn = file.path(  
        p$modeldir, 
        p$bioclass,
        paste("fit_", p$bioclass, ".rdz", sep="") 
    )

    fn_preds = file.path(  
        p$modeldir, 
        p$bioclass,
        paste("fit_preds_", p$bioclass, ".rdz", sep="") 
    )

    fit = NULL
    fit_preds = NULL
    
    fit_theta = NULL
    fit_control_compute = NULL

    if ( "load" %in% todo ) {
        if (file.exists(fn)){
            message( "Loading model fit: ", fn)
            fit = read_write_fast(fn)
            if (!is.null(fit)) return(fit)
        }
        message("\n\nModel fit not found ... will regenerate")
    } 
    
    if ( "load_preds" %in% todo ) {
        if (file.exists(fn_preds)){
            message( "Loading model fit_preds: ", fn_preds)
            fit_preds = read_write_fast(fn_preds)
            if (!is.null(fit_preds)) return(fit_preds)
        }
        message("\n\nModel fit_preds not found ... will regenerate")
    
    }
    
    if ( "fit" %in% todo ) {
          
        message( "\n---\n---\n---\n---\n")
        message("Creating model fit: ", fn, "\n")
            
        M = model_size_data_carstm( p=p )  

        # NOTE: we fit on observations only .. (due to problem size and optimizer difficulties)
        # the predictions for the pred surface must be done separately 
        M = M[ tag=="observations", ]

        # summary on link scale
        MS = summarize_observations(
            obs = M[, pa ],
            offsets = M[, data_offset ], 
            family="binomial"
        )  
        
        ## these must be in global namespace for INLA
        
        sppoly <<- attributes(M)$sppoly

        H <<- inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] ) 

        if (is.null(theta0)) {
          if (exists( "theta", p )) theta0 = p[["theta"]][[ p$bioclass ]]
        }

        message("Trying defaults with inla mode: ",  inla.mode )
        message( "\n---\n---\n" )

        fit = try( inla( 
            data=M, 
            formula=p$formula, 
            family="binomial", 
            verbose=verbose,
            inla.mode=inla.mode,
            control.inla = control.inla,
            control.predictor = list(compute = TRUE,  link = 1), 
            control.compute=control.compute,
            control.mode= list( theta= theta0, restart=restart ),
            num.threads=num.threads
        ), silent=TRUE )

    
        if (inherits(fit, "try-error" )) {
        if (inla.mode != "classic" ) {
          message( "\n---\n---\n" )
          message("Trying more robust settings with inla mode: classic ")
          message( "\n---\n---\n" )
      
          fit = try( inla( 
            data=M, 
            formula=p$formula, 
            family="binomial", 
            safe = TRUE,
            verbose=verbose,
            inla.mode="classic",
            control.inla = list( strategy = "gaussian", fast=TRUE, h=0.01, int.strategy="eb", force.diagonal=TRUE,  cmin=0.0001 ),
            control.predictor = list(compute = TRUE,  link = 1), 
            control.compute = control.compute,
            control.mode= list( theta= theta0, restart=restart ),
            num.threads=num.threads
          ), silent=TRUE )
      }
      }

      if (inherits(fit, "try-error" )) {
          message( "\n---\n---\n" )
          message( "Trying with random start and inla mode: ", inla.mode )
          fit = try( inla( 
            data=M, 
            formula=p$formula, 
            family="binomial", 
            safe = FALSE,
            verbose=verbose,
            inla.mode=inla.mode,
            control.inla = list( strategy = "auto", fast=FALSE, int.strategy="auto" ),
            control.predictor = list(compute = TRUE,  link = 1), 
            control.compute = control.compute,
            num.threads=num.threads
          ), silent=TRUE )
        }



        if (inherits(fit, "try-error" )) {
        if (inla.mode != "classic" ) {
          message( "\n---\n---\n" )
          message( "Trying more stable and robust settings and debugging mode ... with inla mode: classic ")
          message( "\n---\n---\n" )

          fit = try( inla( 
            data=M, 
            formula=p$formula, 
            family="binomial", 
            safe = FALSE,
            verbose=verbose,
            debug = TRUE,
            inla.mode="classic",
            control.inla = list( strategy="laplace", h=0.01, force.diagonal=TRUE, cmin=0.001 ),
            control.predictor = list(compute = TRUE,  link = 1), 
            control.compute= control.compute,
            num.threads=num.threads
          ), silent=TRUE )
        }
        }
    
        if (inherits(fit, "try-error" )) {
          message( "\n---\n---\n" )
          message( "Model fit not converging! Adjust control.inla() options or check data. Giving up ...")
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

        fit_theta = fit$mode$theta
        fit_control_compute = fit$control.compute 
        fit = NULL; gc()

    }


    if ( "fit_preds" %in% todo ) {
        message( "\n---\n---\n---\n---\n")
        message("Creating model fit: ", fn, "\n")
            
        M = model_size_data_carstm( p=p )  
        sppoly <<- attributes(M)$sppoly

        H <<- inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] ) 

        if (is.null(fit_theta)) {
          fit = read_write_fast( fn=fn )
          fit_theta = fit$mode$theta
          fit_control_compute = fit$control.compute 
          fit = NULL
          gc()
        }

        # here we work on the prediction surface
        M = M[ tag=="predictions", ]

        fit_preds = inla( 
            data=M, 
            formula=p$formula, 
            family="binomial", 
            safe = FALSE,
            verbose=verbose,
            debug = FALSE,
            inla.mode="compact",
            control.predictor = list(compute = TRUE,  link = 1), 
            control.compute= fit_control_compute,
            control.mode = list( theta= fit_theta, restart=FALSE )
            num.threads=num.threads
        )

        fit_preds$modelinfo = list(
            bioclass = p$bioclass, 
            H = H,
            MS = MS,
            theta0 = theta0,
            fn = fn_preds
        )

        
        message( "\n---\n---\n" )
        message( "Saving model fit_preds: ", fn_preds, "\n")
        message( "\n---\n---\n" )

        read_write_fast( data=fit_preds, fn=fn_preds )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave
    }

    return("Done")
}

