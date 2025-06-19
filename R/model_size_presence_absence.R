

model_size_presence_absence = function(p, sexid, sppoly, action="load") {

    fn = file.path( 
        p$modeldir, 
        sexid,
        paste("fit_", sexid, "_", paste0(p$span, collapse="_"), ".rds", sep="") 
    )
    print(fn)

    fit = NULL
    if (action != "redo" ) {
        if (file.exists(fn)){
            fit = read_write_fast(fn)
            if (!is.null(fit)) return(fit)
        }
    }
    
    M = model_size_data_carstm(p=p, sexid=sexid, sppoly=sppoly, redo=FALSE, carstm_set_redo=FALSE )  
    
    # summary on link scale
    MS = summarize_observations(
        obs = M[tag=="observations", pa ],
        offsets = M[tag=="observations", data_offset ], 
        family="binomial"
    )  

    H = inla_hyperparameters(  reference_sd = MS[["sd"]], alpha=0.5, MS[["median.50%"]] )  

    # model pa
 
    fit = inla( 
        formula=p$formula, 
        data=M, 
        family="binomial", 
        verbose=TRUE, 
        # control.inla = list( strategy="adaptive", int.strategy="eb", h=0.05 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        num.threads="4:3"
    )

    read_write_fast( fit, file=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)

}

