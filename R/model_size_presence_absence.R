

model_size_presence_absence = function( p, theta0=NULL, todo="load", num.threads="4:3", improve.fit=FALSE) {
     
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
    
    message( "Creating model fit: ", fn)
         
    M = model_size_data_carstm( p=p )  

    if (p$zero_flag == "unit_zeros") {
        k = which(M$tag == "observations" & M$crabno==0)
        if (length(k) > 0) M$data_offset[k] = 1    
    }
    
    # summary on link scale
    MS = summarize_observations(
        obs = M[tag=="observations", pa ],
        offsets = M[tag=="observations", data_offset ], 
        family="binomial"
    )  
    
    ## these must be in global namespace for INLA
    
    sppoly <<- attributes(M)$sppoly

    H <<- inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] ) 

    fit = inla( 
        data=M, 
        formula=p$formula, 
        family="binomial", 
        verbose=TRUE, 
        # control.inla = list( strategy="adaptive", int.strategy="eb", h=0.05 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    )

    if (improve.fit) {
        message( "Improving model fit: ", fn)
        fit = inla.hyperpar( fit, verbose=TRUE, restart=TRUE )  #  dz = 0.75, diff.logdens = 15,
    }

    fit$modelinfo = list(
        bioclass = p$bioclass, 
        H = H,
        MS = MS,
        theta0 = theta0,
        fn = fn
    )
    
    message( "Saving model fit: ", fn)
    read_write_fast( data=fit, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)

}

