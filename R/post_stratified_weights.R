
post_stratified_weights = function(p, todo="load", nposteriors=5000, mc.cores=1 ) {
 
    span = p$span(p$bioclass)

    vn = paste("post_stratified_weights", "_", paste0(span, collapse="_"), ".rdz", sep="") 

    vns = paste("post_stratified_weights_samples", "_", paste0(span, collapse="_"), ".rdz", sep="") 
 
    vnsb = paste("post_stratified_weights_samples_bias_adjusted", "_", paste0(span, collapse="_"), ".rdz", sep="") 
 
    fnout = file.path( p$modeldir, p$bioclass, vn ) 
    fnout_samples = file.path( p$modeldir, p$bioclass, vns ) 
    fnout_samples_bias = file.path( p$modeldir, p$bioclass, vnsb )

    O = NULL

    if ( "load" %in% todo ){
        message( "\nLoading from file: ", fnout)
        if (file.exists(fnout)) {
            O = read_write_fast(fnout)
            if (!is.null(O)) {
                return(O)
            }
        }
    }

    if ( "samples" %in% todo ){
        message( "\nLoading samples from file: ", fnout_samples)
        if (file.exists(fnout_samples)) {
            O = read_write_fast(fnout_samples)
            if (!is.null(O)) {
                return(O)
            }
        }
    }

    if ( "bias_adjusted" %in% todo ){
        message( "\nLoading samples from file: ", fnout_samples_bias)
        if (file.exists(fnout_samples_bias)) {
            O = read_write_fast(fnout_samples_bias)
            if (!is.null(O)) {
                return(O)
            }
        }
    }


    message( "\nCreating post-stratification weights for: ", p$bioclass)

    M = model_size_data_carstm( p=p )  
    M$year = factor2number( M$year )

    iobs = which(M$tag == "observations")
    ipreds = which(M$tag == "predictions")
    
    O = M[iobs, ]   # observations (individual, i)
    P = M[ipreds, ] # predictions (at areal units, a)
    M = NULL; gc()

    fit = model_size_presence_absence( p=p, todo="load" ) 

    # fetch direct predictions of mean, sd
    O$individual_prob_mean = fit$summary.fitted.values[["mean"]][iobs]
    O$individual_prob_sd   = fit$summary.fitted.values[["sd"]][iobs]

    P$auid_prob_mean = fit$summary.fitted.values[["mean"]][ipreds]
    P$auid_prob_sd   = fit$summary.fitted.values[["sd"]][ipreds]

    P = P[, .(AUID, year, cyclic, cwd, mat, auid_prob_mean, auid_prob_sd)]
    P = P[ , PID := do.call(paste, .SD), .SDcols = c("AUID", "year", "cyclic", "cwd", "mat")]

    # obtain joint-posterior samples

    message( "\nExtracting samples of joint posteriors")

    S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores )
        
    fit = NULL; gc()

    message( "\nPost-processing/reformatting/merging")

    # observations (i) and point estimates
    O = O[, .(
        AUID, year, cyclic, cwd, mat, 
        z, substrate.grainsize, dyear, t, pca1, pca2, sid,
        totno, totwgt, meansize, data.source, gear, cw, mass, data_offset,
        individual_prob_mean, individual_prob_sd
    )]

    # O$year = factor2number( O$year )
    # P$year = factor2number( P$year ) 

    # add associated areal unit level predictions (a)
    O = P[ O, on=.(AUID, year, cyclic, cwd, mat)]
    ip = match( O$PID, P$PID )  # get correct row order of areal units to match observations
    P = NULL
    gc()

    # note: 
    # we will add surface area outside of this function 
    # as there exist multiple surface areas, 
    # depending upon region of interest (W_a), cfanorth, south etc
    # and so need to be determined at time of aggregation
    
    # means are already on user scale .. probability 
    # this is the overall point estimate ... still needs multiplication with SA 
    O$post_stratified_ratio = O$auid_prob_mean / O$individual_prob_mean   

    attr(O, "bioclass" ) = p$bioclass
    attr(O, "span" ) = p$span(p$bioclass)
    attr(O, "formula" ) = p$formula 

    message( "\nSaving", fnout )

    read_write_fast( O, fn=fnout )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    O = NULL
    gc()

    message( "\nPosterior samples")

    for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info 

    fkk = inla_get_indices("Predictor", tag=tag, start=start, len=length, model="direct_match")
    fkk = unlist(fkk)
    ndat = length(fkk)
    Osamples = array(NA, dim=c( ndat,  nposteriors ) )
    for (i in 1:nposteriors) {
        Osamples[,i] =  S[[i]]$latent[fkk,] 
    }
   
    Osamples = inverse.logit( Osamples[ipreds[ip],] ) / inverse.logit(Osamples[iobs,]) # same order as O, samples of the ratio of probabilities (a/i)
 
    attr(Osamples, "bioclass" ) = p$bioclass
    attr(Osamples, "span" ) = p$span(p$bioclass)
    attr(Osamples, "formula" ) = p$formula 

    message( "\nSaving", fnout_samples )

    read_write_fast( Osamples, fn=fnout_samples )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    Osamples = NULL; gc()
    
    fss = inla_get_indices("inla.group(cwd, method = \"quantile\", n = 13)", tag=tag, start=start, len=length, model="direct_match")
    fss = unlist(fss)
    nobs = length(fss)
    Osamples_bias = array(NA, dim=c( nobs,  nposteriors ) )
     
    for (i in 1:nposteriors) {
        Osamples_bias[,i] = S[[i]]$latent[fss,] 
    }
    
    Osamples_bias = inverse.logit(Osamples_bias )
    
    S = fss = NULL ;  gc()
 
    attr(Osamples_bias, "bioclass" ) = p$bioclass
    attr(Osamples_bias, "span" ) = p$span(p$bioclass)
    attr(Osamples_bias, "formula" ) = p$formula 

    message( "\nSaving", fnout_samples_bias )

    read_write_fast( Osamples_bias, fn=fnout_samples_bias )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave


    return(fnout)

}

