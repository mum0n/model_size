
post_stratified_weights = function(p, todo="load", nposteriors=5000, mc.cores=1 ) {
 
    fnout = file.path( p$modeldir, p$bioclass, "post_stratified_weights.rdz" ) 
    fnout_samples_obs = file.path( p$modeldir, p$bioclass, "post_stratified_weights_samples_obs.rdz" ) 
    fnout_samples_preds = file.path( p$modeldir, p$bioclass, "post_stratified_weights_samples.rdz" ) 
    fnout_samples_preds2 = file.path( p$modeldir, p$bioclass, "post_stratified_weights_samples2.rdz" ) 
    fnout_samples_bias = file.path( p$modeldir, p$bioclass, "post_stratified_weights_samples_bias_adjusted.rdz" )

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

 

    if ( "observation_samples" %in% todo ){
        message( "\nLoading obserevation samples from file: ", fnout_samples_obs)
        if (file.exists(fnout_samples_obs)) {
            O = read_write_fast(fnout_samples_obs)
            if (!is.null(O)) {
                return(O)
            }
        }
    }
 

    if ( "prediction_samples" %in% todo ){
        message( "\nLoading prediction samples from file (at observations): ", fnout_samples_preds)
        if (file.exists(fnout_samples_preds)) {
            O = read_write_fast(fnout_samples_preds)
            if (!is.null(O)) {
                return(O)
            }
        }

    }
    if ( "prediction_samples2" %in% todo ){
        message( "\nLoading prediction samples from file (at prediction time-slice): ", fnout_samples_preds2)
        if (file.exists(fnout_samples_preds2)) {
            O = read_write_fast(fnout_samples_preds2)
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

    # operate upon fit first and remove from meory to reduce RAM demand 
    # operate upon fit first and remove from meory to reduce RAM demand 
    fit = model_size_presence_absence( p=p, todo="load" ) 
    fit_summ_mean = fit$summary.fitted.values[["mean"]] 
    fit_summ_sd = fit$summary.fitted.values[["sd"]] 
    S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores )
    fit = NULL; gc() # no longer needed


    M = model_size_data_carstm( p=p )  
    
    
    M$year = factor2number( M$year )

    # fetch direct predictions of mean, sd
    M$individual_prob_mean = fit_summ_mean
    M$individual_prob_sd   = fit_summ_sd 

    fit_summ_sd = fit_summ_mean = NULL 

    iobs = which(M$tag == "observations")
    ipreds = which(M$tag == "predictions")
    
    setDT(M)
    setDT(M)

    # observations (individual, i)
    O = M[iobs, .(
        kuid, AUID, sid, year, cyclic, cwd, mat, sex,
        z, substrate.grainsize, dyear, dyri, t, pca1, pca2, sid,
        cw, mass, data_offset,
        individual_prob_mean, individual_prob_sd
    )]
    
    # matches to prediction time slice ( 1 September )
    O$cyclic_pred = p$prediction_dyear_index
    kuid2 = c("AUID", "year", "cyclic_pred", "cwd", "mat", "sex")
    O = O[, kuid_pred := do.call(paste, .SD), .SDcols = kuid2 ]
    O$cyclic_pred = NULL

    # get correct row order of areal units to match observations
    ip  = ipreds[ match( O$kuid,  M$kuid[ipreds] ) ]    # predictions (at areal units, a) on M

    ip2 = ipreds[ match( O$kuid_pred, M$kuid[ipreds] )]  # prediction time-slice on M
 
    # add associated areal unit level predictions (a)
    O$auid_prob_mean = M$individual_prob_mean[ ip ]
    O$auid_prob_sd = M$individual_prob_sd[ ip ]
    O$post_stratified_ratio = O$auid_prob_mean / O$individual_prob_mean   
 
    O$auid_prob_mean2 = M$individual_prob_mean[ ip2 ]
    O$auid_prob_sd2 = M$individual_prob_sd[ ip2 ]
    O$post_stratified_ratio2 = O$auid_prob_mean2 / O$individual_prob_mean   

    M = NULL; gc()
 
    gc()

    # note: 
    # we will add surface area outside of this function 
    # as there exist multiple surface areas, 
    # depending upon region of interest (W_a), cfanorth, south etc
    # and so need to be determined at time of aggregation
    
    # means are already on user scale .. probability 
    # this is the overall point estimate ... still needs multiplication with SA 

    attr(O, "bioclass" ) = p$bioclass
    attr(O, "formula" ) = p$formula 

    message( "\nSaving", fnout )

    read_write_fast( O, fn=fnout )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    O = NULL
    gc()

    # obtain joint-posterior samples

    message( "\nExtracting samples of joint posteriors")


    for (z in c("tag", "start", "length") ) {
        assign(z, attributes(S)[[".contents"]][[z]] )  # index info 
    }

    fkk = inla_get_indices(
        "Predictor", 
        tag=tag, 
        start=start, 
        len=length, 
        model="direct_match"
    )
    
    fkk = unlist(fkk)
    ndat = length(fkk)
    Osamples = array(NA, dim=c( ndat,  nposteriors ) )
    for (i in 1:nposteriors) {
        Osamples[,i] =  S[[i]]$latent[fkk,] 
    }
    

    message( "\nSaving", fnout_samples_obs )

    # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave
    read_write_fast( Osamples[iobs,], fn=fnout_samples_obs )  # on logit scale

    message( "\nSaving", fnout_samples_preds )

    # same order as O, samples of the ratio of probabilities (a/i)
    read_write_fast( Osamples[ ip,], fn=fnout_samples_preds )  # on logit scale
  
    message( "\nSaving", fnout_samples_preds2 )

    read_write_fast( Osamples[ ip2,], fn=fnout_samples_preds2 )  # on logit scale


    fss = inla_get_indices(
        "inla.group(cwd, method = \"quantile\", n = 11)", 
        tag=tag, start=start, 
        len=length, 
        model="direct_match"
    )
    
    fss = unlist(fss)
    nobs = length(fss)
    Osamples = array(NA, dim=c( nobs,  nposteriors ) )
     
    for (i in 1:nposteriors) {
        Osamples[,i] = S[[i]]$latent[fss,] 
    }
    
    Osamples = inverse.logit(Osamples )
    
    S = fss = NULL ;  gc()
  
    message( "\nSaving", fnout_samples_bias )

    read_write_fast( Osamples, fn=fnout_samples_bias )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave


    return(fnout)

}

