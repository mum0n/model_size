

model_size_presence_absence = function( sexid, span, modeldir, formula=NULL, theta0=NULL, 
    action="load", num.threads="4:3", improve.fit=FALSE) {
     
    fn = file.path( 
        p$modeldir, 
        sexid,
        paste("fit_", sexid, "_", paste0(span, collapse="_"), ".rdz", sep="") 
    )
    

    fit = NULL
    if (action != "redo" ) {
        if (file.exists(fn)){
            message( "Loading model fit: ", fn)
            fit = read_write_fast(fn)
            if (!is.null(fit)) return(fit)
        }
    }
    
    message( "Creating model fit: ", fn)
         
    M = model_size_data_carstm( p=p, sexid=sexid )  

    # summary on link scale
    MS = summarize_observations(
        obs = M[tag=="observations", pa ],
        offsets = M[tag=="observations", data_offset ], 
        family="binomial"
    )  
    
    ## these must be in global namespace for INLA
    sppoly <<- attributes(M)$sppoly
    H <<- inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] )  
    # model pa solutions close to final:

    # female:  
# maxld= -207571.2890 fn=2264 theta= -3.4225 -0.1573 1.0692 1.4218 2.7243 0.0157 0.2088 -1.8864 -3.0328 0.5815 2.4277 -2.1472 2.8986 -1.5248 -0.7769 0.9287 [10.65, 3.222]

    # male: 
# maxld= -341990.7511 fn=769 theta= -3.3990 -1.1159 1.2345 1.7680 1.5166 0.0099 0.9860 -0.6753 -0.3107 4.9739 3.4903 -1.9714 4.3707 -1.0976 -0.2291 0.8286 [10.49, 3.884]
# maxld= -341990.6576 fn=881 theta= -3.4036 -1.1141 1.2347 1.7686 1.5155 0.0099 0.9877 -0.6803 -0.3101 4.9744 3.4918 -1.9614 4.3707 -1.0946 -0.2417 0.8324 [10.50, 3.887]
# maxld= -341989.7669 fn=2065 theta= -3.7692 -0.9831 1.2247 1.7921 1.4443 0.0036 1.1232 -1.0129 -0.2673 5.0086 3.6189 -1.9536 4.3731 -1.0818 -0.3058 0.8301 [11.08, 3.870]
    if (is.null(theta0)) {
        theta0 = switch(
            sexid,
            female = c(-3.4225,-0.1573,1.0692,1.4218,2.7243,0.0157,0.2088,-1.8864,
                        -3.0328,0.5815,2.4277,-2.1472,2.8986,-1.5248,-0.7769,0.9287 ),
            male   = c(-3.3897, -1.1157, 1.2352, 1.7682, 1.5176, 0.0099, 0.9847, -0.6698,
                        -0.3111, 4.9736, 3.4893, -1.9710, 4.3707, -1.0877, -0.2271, 0.8203)
        )
    }


   fit = inla( 
        data=M, 
        formula=formula, 
        family="binomial", 
        verbose=TRUE, 
        # control.inla = list( strategy="adaptive", int.strategy="eb", h=0.05 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    )

    if (improve.fit) fit = inla.hyperpar( fit, verbose=TRUE, restart=TRUE )  #  dz = 0.75, diff.logdens = 15,
    
    fit$modelinfo = list(
        sexid = sexid, 
        H = H,
        MS = MS,
        theta0 = theta0,
        fn = fn
    )
    
    read_write_fast( data=fit, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)

}

