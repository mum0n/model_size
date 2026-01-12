model_size_results = function(p, todo="post_stratified_weights", 
  nposteriors=1000, mc.cores=1, only_observations=TRUE,   
  working.directory = tempfile(pattern = "to_delete_", tmpdir =p$project_output_directory, fileext = "")  ) {

    message("Files are being stored in the output directory. Remember to delete them once finished.\n", working.directory, "\n" )

    # all classes (observations) 
    fn_post_stratified_weights = file.path( p$modeldir, "post_stratified_results.rdz" ) 
    fn_size_selectivity = file.path( p$modeldir, "size_selectivity.rdz" )
    fn_fixed_effects = file.path( p$modeldir, "fixed_effects.rdz" )

    # class-specific results
    if (exists("bioclass", p)) {
      fn_loc = file.path( p$modeldir, p$bioclass )
      fn_observation_ratios = file.path( fn_loc, "post_stratified_ratios.rdz" ) 
      fn_observation_samples = file.path( fn_loc, "post_stratified_weights_samples_obs.rdz" ) 
      fn_psratio_samples = file.path( fn_loc, "post_stratified_ratio_samples.rdz" )  
      fn_psratio_samples_obs = file.path( fn_loc, "post_stratified_weights_samples_pred_obs.rdz" ) 
      fn_size_selectivity_samples = file.path( fn_loc, "post_stratified_size_selectivity_samples.rdz" )
      fn_fixed_effects_samples =  file.path( fn_loc, "fixed_effects_samples.rdz" )
    }

    fno = switch( todo,
      observation_ratios = fn_observation_ratios,
      observation_samples = fn_observation_samples,
      prediction_samples = fn_psratio_samples,
      prediction_samples_obs = fn_psratio_samples_obs,
      size_selectivity_samples = fn_size_selectivity_samples,
      fixed_effects_samples = fn_fixed_effects_samples,
      NULL
    )

    out = NULL
    size_selectivity = list()
    fixed_effects = list()

    if (!is.null(fno)) {
      message( "\nLoading from file: ", fno)
      if (file.exists(fno)) out = read_write_fast(fno)
      if (!is.null(out)) return(out)
    }

    if ( todo == "post_stratified_weights" ) {
        post_stratified_weights = read_write_fast(fn_post_stratified_weights) # n=1106430   
        if (only_observations) {
          # remove non-observations: 
          # when crabno==0 they represent data points with no observations but predictions are still made .. 
          # to be conservative we ignore them  .. the associated probabilities can be used as a post-stratified "weight" to zero-values 
          post_stratified_weights = post_stratified_weights[ crabno > 0, ]  # n=532176   ]
        }
        return( post_stratified_weights )
    }
 
    
    if ( todo == "size_selectivity" ) {
        size_selectivity = read_write_fast(fn_size_selectivity)
        return( size_selectivity )
    }

    if (todo == "fixed_effects" ) {
        fixed_effects = read_write_fast(fn_fixed_effects)
        return( fixed_effects )
    }

    if ( todo == "post_stratified_weights_redo" ){
      # assemble all bioclasses 
      # surface area
      pg = areal_units(p=p, areal_units_directory=p$project_data_directory)

      sa_vars = c(
        "cfanorth_surfacearea", 
        "cfasouth_surfacearea", 
        "cfa23_surfacearea", 
        "cfa24_surfacearea", 
        "cfa4x_surfacearea", 
        "au_sa_km2"
      )
      
      for (i in sa_vars) {
          attr(pg[,i], "units") = NULL
          class(pg[,i]) = NULL
      }
      
      tokeep = c("AUID", sa_vars, "strata_to_keep" )
      pg = data.table(pg)[, ..tokeep]
      attr(pg[["au_sa_km2"]], "units") = NULL
      class(pg[["au_sa_km2"]]) = NULL
      
      names(pg) = c("AUID", "cfanorth_sa", "cfasouth_sa", "cfa23_sa", "cfa24_sa", "cfa4x_sa", "cfaall_sa", "strata_to_keep")

      for (bc in c("f.imm", "f.mat", "m.imm", "m.mat") ) {
        
        p$bioclass = bc

        fn_loc = file.path( p$modeldir, p$bioclass )
        fn_observation_ratios = file.path( fn_loc, "post_stratified_ratios.rdz" ) 
        fn_observation_samples = file.path( fn_loc, "post_stratified_weights_samples_obs.rdz" ) 
        fn_psratio_samples = file.path( fn_loc, "post_stratified_weights_samples_pred.rdz" ) 
        fn_psratio_samples_obs = file.path( fn_loc, "post_stratified_weights_samples_pred_obs.rdz" ) 
        fn_size_selectivity_samples = file.path( fn_loc, "post_stratified_size_selectivity_samples.rdz" )
        fn_fixed_effects_samples =  file.path( fn_loc, "fixed_effects_samples.rdz" )

        # operate upon "fit" first and remove from memory to reduce RAM demand 
        message("\nSampling and extracting for: ", bc)

        inla.setOption( working.directory=working.directory )

        fit = model_size_presence_absence( p=p, todo="load_preds" ) 

        S = inla.posterior.sample( nposteriors, fit, add.names=FALSE, num.threads=mc.cores )
        
        fit_summ_mean = fit$summary.fitted.values[["mean"]] 
        fit_summ_sd = fit$summary.fitted.values[["sd"]] 
        fixeff = fit$summary.fixed

        vn_rnd = names(fit$summary.random)
        vn_cwd = vn_rnd[ grep(".*cwd,", vn_rnd) ]
        if (length(vn_cwd) !=1 ) {
          message("There is a problem with variable name expectations")
          browser()
        }
        
        sizeselect = fit$summary.random[[vn_cwd]]
        # plot( (sizeselect[,2]) ~ exp(sizeselect[,1]))  # log odds ratio
        # plot(exp(sizeselect[,2])~ exp(sizeselect[,1]))  # odds ratio
        # plot(1/exp(sizeselect[,2])~ exp(sizeselect[,1])) # selectivity ratio
        # input data 
        M = fit$.args$data
        # M = model_size_data_carstm( p=p )  
        
        M$year = factor2number( M$year )

        # fetch direct predictions of mean, sd
        M$individual_prob_mean = fit_summ_mean
        M$individual_prob_sd   = fit_summ_sd 

        fit = NULL; gc() # no longer needed
        fit_summ_sd = fit_summ_mean = NULL 

        iobs = which(M$tag == "observations")
        ipreds = which(M$tag == "predictions")
        
        setDT(M)
  
        # observations (individual, i)
        O = M[iobs, .(
            kuid, AUID, sid, crabno, year, cyclic, cwd, mat, sex,
            z, substrate.grainsize, dyear, dyri, t, pca1, pca2,
            cw, mass, data_offset,
            individual_prob_mean, individual_prob_sd
        )]
        
        # matches to prediction time slice ( 1 September )
        O$cyclic_pred = p$prediction_dyear_index
        kuid2 = c("AUID", "year", "cyclic_pred", "cwd", "mat", "sex")
        O = O[, kuid_pred := do.call(paste, .SD), .SDcols = kuid2 ]
        O$cyclic_pred = NULL

        # get correct row order of areal units to match observations
        ip  = ipreds[ match( O$kuid,      M$kuid[ipreds] )] # predictions (at areal units, a) on M
        ipo = ipreds[ match( O$kuid_pred, M$kuid[ipreds] )] # prediction time-slice on M
    
    
        # add associated areal unit level predictions (a)
        O$auid_prob_mean = M$individual_prob_mean[ ip ]
        O$auid_prob_sd = M$individual_prob_sd[ ip ]
    
        O$auid_prob_mean2 = M$individual_prob_mean[ ipo ]
        O$auid_prob_sd2 = M$individual_prob_sd[ ipo ]
        O$bioclass = bc
        
        M=NULL; gc()

        # note: 
        # we will add surface area outside of this function 
        # as there exist multiple surface areas, 
        # depending upon region of interest (W_a), cfanorth, south etc
        # and so need to be determined at time of aggregation
        
        # means are already on user scale .. probability 
        # this is the overall point estimate ... still needs multiplication with SA 

        attr(O, "bioclass" ) = p$bioclass
        attr(O, "formula" ) = p$formula 
        attr(O, "fixed_effects") = fixeff
        attr(O, "size_selectivity") = sizeselect

        message( "\nSaving: ", fn_observation_ratios )
        read_write_fast( O, fn=fn_observation_ratios )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

        out = rbind( out, O )  # retain for aggregate save at the end

        O = NULL;         gc()

        # obtain joint-posterior samples

        message( "\nSaving samples of joint posteriors")

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
        fkk = NULL
        message( "\nSaving:  ", fn_observation_samples )
        read_write_fast( Osamples[iobs,], fn=fn_observation_samples )  # on logit scale


        # same order as O, samples of the ratio of probabilities (a/i)
        post_stratified_ratio  = Osamples[ ip,] / Osamples[iobs,]  
        
        
        message( "\nSaving:  ", fn_psratio_samples )
        read_write_fast( post_stratified_ratio, fn=fn_psratio_samples )  # on logit scale
        post_stratified_ratio = NULL


        # prediction time-slice on M
        post_stratified_ratio_obs = Osamples[ ipo,] / Osamples[iobs,]   
        message( "\nSaving:  ", fn_psratio_samples_obs )
        read_write_fast( post_stratified_ratio_obs, fn=fn_psratio_samples_obs )  # on logit scale
        
        post_stratified_ratio_obs = NULL
        Osamples = NULL 
        gc()

        fe = inla_get_indices(
            "(Intercept)", 
            tag=tag, 
            start=start, 
            len=length, 
            model="direct_match"
        )

        fe = unlist(fe)
        nobs = length(fe)
        Osamples = array(NA, dim=c( nobs,  nposteriors ) )
        
        for (i in 1:nposteriors) {
            Osamples[,i] = S[[i]]$latent[fe,]   # log odds ratio
        }
        fe = NULL
        
        message( "\nSaving:  ", fn_fixed_effects_samples )
        read_write_fast( Osamples, fn=fn_fixed_effects_samples )   # on logit scale
        Osamples = NULL

        fss = inla_get_indices(
            vn_cwd, 
            tag=tag, start=start, 
            len=length, 
            model="direct_match"
        )
        
        fss = unlist(fss)
        nobs = length(fss)
        Osamples = array(NA, dim=c( nobs,  nposteriors ) )
        
        for (i in 1:nposteriors) {
            Osamples[,i] = S[[i]]$latent[fss,]   # log odds ratio
        }
        
        fss = NULL
        S = NULL 

        message( "\nSaving:  ", fn_size_selectivity_samples )
        read_write_fast( Osamples, fn=fn_size_selectivity_samples )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave
        
        Osamples = NULL
        gc()
      
        size_selectivity[[bc]] = sizeselect
        fixed_effects[[bc]] = fixeff

      }

      out = pg[out, on="AUID"] # merge SA's
      
      # variable numbers of sets in each AUID needs to be accounted
      # counting the number of sets and dividing the effective wgts 
      # rescales the weigts to a per AUID basis ... done at the script level
      ii = out[ , .(n_stations=length(unique(sid))), by=.(AUID, year) ]  
      out = ii[ out, on=.(AUID, year) ]
      out[ is.na(n_stations), "n_stations"] = 1

      out$post_stratified_ratio_obs  = out$auid_prob_mean  / out$individual_prob_mean    # observation time
      out$post_stratified_ratio = out$auid_prob_mean2 / out$individual_prob_mean    # time shifted 

      read_write_fast( out, fn=fn_post_stratified_weights )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave
      read_write_fast( size_selectivity, fn=fn_size_selectivity)
      read_write_fast( fixed_effects, fn=fn_fixed_effects)

      return(out)
    }

}


