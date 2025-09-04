
post_stratified_weights = function(p, toget="load", nposteriors=5000, mc.cores=1 ) {
 
    span = p$span(p$bioclass)

    vn = paste("post_stratified_weights", "_", paste0(span, collapse="_"), ".rdz", sep="") 

    vns = paste("post_stratified_weights_samples", "_", paste0(span, collapse="_"), ".rdz", sep="") 
 
    fnout = file.path( p$modeldir, p$bioclass, vn ) 
    fnout_samples = file.path( p$modeldir, p$bioclass, vns ) 
    
    O = NULL

    if (toget == "load" ){
        message( "\nLoading from file: ", fnout)
        if (file.exists(fnout)) {
            O = read_write_fast(fnout)
            if (!is.null(O)) {
                return(O)
            }
        }
    }

    if (toget == "samples" ){
        message( "\nLoading samples from file: ", fnout_samples)
        if (file.exists(fnout_samples)) {
            O = read_write_fast(fnout_samples)
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

    fit = model_size_presence_absence( p=p, toget="load" ) 

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

    for (z in c("tag", "start", "length") ) assign(z, attributes(S)[[".contents"]][[z]] )  # index info 

    fkk = inla_get_indices("Predictor", tag=tag, start=start, len=length, model="direct_match")
    fkk = unlist(fkk)
    ndat = length(fkk)

    Osamples = array(NA, dim=c( ndat,  nposteriors ) )
    for (i in 1:nposteriors) {
        # note: "S$latent" is on INLA's internal scale (logit) and needs to be back-transformed
        Osamples[,i] = inverse.logit( S[[i]]$latent[fkk,] )
    }
    S = fkk = NULL ;  gc()

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


    # these are the joint posterior samples of the above
    # saved as an attribute ... note it had the same order as "O"

    O = Osamples[ipreds[ip],] / Osamples[iobs,] # same order as O, samples of the ratio of probabilities (a/i)
    
    ipreds = iobs = ip = Osamples = NULL; gc()
    
    attr(O, "bioclass" ) = p$bioclass
    attr(O, "span" ) = p$span(p$bioclass)
    attr(O, "formula" ) = p$formula 

    message( "\nSaving", fnout_samples )

    read_write_fast( O, fn=fnout_samples )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(O)

}

