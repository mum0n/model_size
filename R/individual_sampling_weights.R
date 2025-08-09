
individual_sampling_weights = function(p, action="load", nposteriors = 100, mc.cores = 2 ) {

    fnout = file.path( p$modeldir, "individual_sampling_weights.rdz" ) 
 
    message( "\nCreating/loading: ", fnout)
    
    out = NULL
    ORR = NULL

    if (action != "redo" ){
        if (file.exists(fnout)) {
            out = read_write_fast(fnout)
            if (!is.null(out)) return(out)
        }
    }
 
    for (sexid in c("female", "male")) {
 
        message( "\nLoading model fit for: ", sexid)

        M = model_size_data_carstm( p=p, sexid=sexid )  
        M$year = factor2number( M$year )

        p$selection$biologicals_using_snowcrab_filter_class = sexid
        p$carstm_model_label = sexid
 
        sppoly = areal_units(p=p)
 
        sa_vars = c("cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea", "cfa24_surfacearea", "cfa4x_surfacearea", "au_sa_km2")
        for (i in sa_vars) {
            attr(sppoly[,i], "units") = NULL
            class(sppoly[,i]) = NULL
        }
        tokeep = c("AUID", sa_vars, "strata_to_keep" )
        pg = data.table(sppoly)[, ..tokeep]
 
        iobs = which(M$tag == "observations")
        ipreds = which(M$tag == "predictions")
        
        O = M[iobs, ]
        P = M[ipreds, ] 
        M = NULL; gc()
  
        fit = model_size_presence_absence( sexid=sexid, span=p$span(sexid), modeldir=p$modeldir, action="load" ) 

        O$individual_prob_mean = fit$summary.fitted.values[["mean"]][iobs]
        O$individual_prob_sd   = fit$summary.fitted.values[["sd"]][iobs]

        P$auid_prob_mean = fit$summary.fitted.values[["mean"]][ipreds]
        P$auid_prob_sd   = fit$summary.fitted.values[["sd"]][ipreds]


        P = P[, .(AUID, year, cyclic, cwd, mat, auid_prob_mean, auid_prob_sd)]
        P = P[ , PID := do.call(paste, .SD), .SDcols = c("AUID", "year", "cyclic", "cwd", "mat")]
 
        S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores )
            
        fit = M = NULL; gc()

        for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info 

        fkk = inla_get_indices("Predictor", tag=tag, start=start, len=length, model="direct_match")
        fkk = unlist(fkk)
        ndat = length(fkk)
        Osamples = array(NA, dim=c( ndat,  nposteriors ) )
        for (i in 1:nposteriors) Osamples[,i] = S[[i]]$latent[fkk,]  
        
        S = NULL ;  gc()

        SO = inverse.logit(Osamples[iobs,])  # same order as "O"
        SP = inverse.logit(Osamples[ipreds,])  #same order as "P"
        
        Osamples = NULL; gc()
 
        O = O[, .(
            AUID, year, cyclic, cwd, mat, 
            z, substrate.grainsize, dyear, t, pca1, pca2, sid,
            totno, totwgt, meansize, data.source, gear, cw, mass, data_offset,
            individual_prob_mean, individual_prob_sd
        )]

        # O$year = factor2number( O$year )
        # P$year = factor2number( P$year ) 

        O = P[ O, on=.(AUID, year, cyclic, cwd, mat)]
        O$sex = sexid

        O$relative_rate = O$individual_prob_mean / O$auid_prob_mean
        out = rbind( out, pg[ O, on="AUID" ] )

        ip = match( O$PID, P$PID )
        
        O =  P = NULL
        gc()

        ORR = rbind(ORR, SO / SP[ip,] )  # same order as O

        SO = SP = NULL
        gc()

    }

    attr(out, paste(sexid, "samples", sep="_" ) ) = ORR
    attr(out, "span" ) = p$span
    attr(out, "formula" ) = p$formula 

    read_write_fast( out, fn=fnout )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(out)

}

