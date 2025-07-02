
individual_sampling_weights = function(p, sppoly, action="load") {

    fnout = file.path( p$modeldir, "individual_sampling_weights.rds" ) 
 
    print(fnout)
    
    out = NULL

    if (action != "redo" ){
        if (file.exists(fnout)) {
            out = read_write_fast(fnout)
            if (!is.null(out)) return(out)
        }
    }

    sa_vars = c("cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea", "cfa24_surfacearea", "cfa4x_surfacearea", "au_sa_km2")
    for (i in sa_vars) {
        attr(sppoly[,i], "units") = NULL
        class(sppoly[,i]) = NULL
    }
    tokeep = c("AUID", sa_vars, "strata_to_keep" )
    pg = data.table(sppoly)[, ..tokeep]

    for (sexid in c("female", "male")) {
        M = model_size_data_carstm(p=p, sexid=sexid, redo=FALSE, carstm_set_redo=FALSE )  

        iobs = which(M$tag == "observations")
        ipreds = which(M$tag == "predictions")
        
        O = M[iobs, ]

        fit = model_size_presence_absence( p=p, sexid=sexid, action="load" ) 

        O$individual_prob_mean = fit$summary.fitted.values[["mean"]][iobs]
        O$individual_prob_sd = fit$summary.fitted.values[["sd"]][iobs]

        P = M[ipreds, ] 
        P$auid_prob_mean = fit$summary.fitted.values[["mean"]][ipreds]
        P$auid_prob_sd = fit$summary.fitted.values[["mean"]][ipreds]
         
        # S = inla.posterior.sample( 100, fit, add.names=FALSE, num.threads=4 )  ## too slow to use

        fit = M = NULL; gc()
        
        P = P[, .(AUID, year, cyclic, cwd, mat, auid_prob_mean, auid_prob_sd)]

        O = O[ pa==1,]

        span = p$span(sexid)
        span[1] = log(span[1])
        span[2] = log(span[2])

        O$cwd = discretize_data( O$logcw, span=span ) 
        O$mat = as.numeric(O$mat)
 

        O = O[, .(
            AUID, year, cyclic, cwd, mat, 
            z, substrate.grainsize, dyear, t, pca1, pca2, sid,
            totno, totwgt, meansize, data.source, gear, cw, mass, data_offset
        )]

        O$year = factor2number( O$year )
        
        P$year = factor2number( P$year )
        P$cwd = factor2number( P$cwd )
        P$mat = factor2number( P$mat )

        O = P[ O, on=.(AUID, year, cyclic, cwd, mat)]

        O$relative_rate = O$individual_prob_mean / O$auid_prob_mean
        out = rbind( out, pg[ O, on="AUID" ] )

        attr(out, paste(sexid, "summary", sep="_" ) ) = summary(fit)
    }

    attr(out, "span" ) = p$span
    attr(out, "formula" ) = p$formula 

    read_write_fast( out, file=fnout )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(out)

}

