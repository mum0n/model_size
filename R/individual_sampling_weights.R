

individual_sampling_weights = function(p, sexid, sppoly, action="load") {

    fnO = file.path(  
        p$modeldir, 
        sexid,
        paste("obs_", sexid, "_", paste0(p$span, collapse="_"), ".rds", sep="") 
    )
    
    print(fnO)

    if (action != "redo" ){
        O = NULL
        if (file.exists(fnO)) {
            O = read_write_fast(fnO)
            if (!is.null(O)) return(O)
        }
    }

    M = model_size_data_carstm(p=p, sexid=sexid, sppoly=sppoly, redo=FALSE, carstm_set_redo=FALSE )  

    iobs = which(M$tag == "observations")
    ipreds = which(M$tag == "predictions")
      
    O = M[iobs, ]

    fit = model_size_presence_absence( p=p, sexid=sexid, sppoly=sppoly, action="load" ) 

    O$fitted_mean = fit$summary.fitted.values[["mean"]][iobs]
    O$fitted_sd = fit$summary.fitted.values[["sd"]][iobs]

    P = M[ipreds, ] 
    P$prediction_mean = fit$summary.fitted.values[["mean"]][ipreds]
    P$prediction_sd = fit$summary.fitted.values[["mean"]][ipreds]
      
    fit = M = NULL; gc()
      
    P = P[, .(AUID, year, cyclic, cwd, mat, prediction_mean, prediction_sd)]

    O = O[ pa==1,]

    span = p$span
    span[1] = log(span[1])
    span[2] = log(span[2])

    O$cwd = discretize_data( O$logcw, span=span ) 
    O$mat = as.numeric(O$mat)

    O = P[ O, on=.(AUID, year, cyclic, cwd, mat)]

    O$relative_rate = O$fitted_mean / O$prediction_mean

    sa_vars = c("cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea", "cfa24_surfacearea", "cfa4x_surfacearea", "au_sa_km2")
    for (i in sa_vars) {
        attr(sppoly[,i], "units") = NULL
        class(sppoly[,i]) = NULL
    }

    pg = data.table(sppoly)[, c("AUID", sa_vars, "strata_to_keep" )]
    
    O = pg[ O, on="AUID" ]

    read_write_fast( O, file=fnO )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(O)

}

