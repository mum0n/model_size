model_size_results = function(p, todo="load_results", nposteriors=5000, mc.cores=1, region=NULL) {

    # all classes (observations) 
    post_stratified_weights = file.path( p$modeldir, "post_stratified_results.rdz" ) 

    # class-specific results
    if (exists("bioclass", p)) {
      fn_loc = file.path( p$modeldir, p$bioclass )
      observation_ratios = file.path( fn_loc, "post_stratified_ratios.rdz" ) 
      observation_samples = file.path( fn_loc, "post_stratified_weights_samples_obs.rdz" ) 
      prediction_samples = file.path( fn_loc, "post_stratified_weights_samples.rdz" ) 
      prediction_samples2 = file.path( fn_loc, "post_stratified_weights_samples2.rdz" ) 
      bias_adjusted = file.path( fn_loc, "post_stratified_weights_samples_bias_adjusted.rdz" )
    }

    O = NULL
    size_selectivity = list()

    fno = switch( todo,
      post_stratified_weights = post_stratified_weights,
      observation_ratios = observation_ratios,
      observation_samples = observation_samples,
      prediction_samples = prediction_samples,
      prediction_samples2 = prediction_samples2,
      bias_adjusted = bias_adjusted,
      NULL
    )

    if (!is.null(fno)) {
      message( "\nLoading from file: ", fno)
      if (file.exists(fno)) O = read_write_fast(fno)
      if (!is.null(O)) return(O)
    }

    
    if ("size_selectivity" %in% todo) {
        O = model_size_results(p=p, todo="post_stratified_weights" ) 
        size_selectivity = attr(O, "size_selectivity")
        if (!is.null(region)) size_selectivity = size_selectivity[[region]]
        return( size_selectivity )
    }


    if ( "post_stratified_weights_redo" %in% todo ){
      # assemble all bioclasses 
      # surface area
      pg = areal_units(p=p)

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

      for (bc in c("f.imm", "f.mat", "m.imm", "m.mat") ){
          print(bc)
          p$bioclass = bc
          o = NULL
          o = model_size_results( p=p, todo="observation_ratios" )
          ss = size_selectivity[[bc]] = attr(o, "size_selectivity") 
          # plot( (ss[,2]) ~ exp(ss[,1]))   log odds ratio
          # plot(exp(ss[,2])~ exp(ss[,1]))  # odds ratio
          # plot(1/exp(ss[,2])~ exp(ss[,1])) # selectivity ratio

          oddsratio = spline(
              x = exp(ss[,1]), 
              y = exp(ss[,2]), 
              xout=o$cw 
          )
          o$odds_ratio = oddsratio$y
          oddsratio = NULL
          o = o[, .(
              individual_prob_mean,
              individual_prob_sd,
              auid_prob_mean,
              auid_prob_sd,
              auid_prob_mean2,
              auid_prob_sd2,
              odds_ratio
          )]

          M = NULL
          M = model_size_data_carstm( p=p )  
          M = M[tag=="observations",]

          O = rbind( O, cbind(M, o) )
          M = NULL
          o = NULL
      }

      O$post_stratified_ratio  = O$auid_prob_mean  / O$individual_prob_mean   
      O$post_stratified_ratio2 = O$auid_prob_mean2 / O$individual_prob_mean   

      attr(O, "pg") = pg

      attr(O, "size_selectivity") = size_selectivity

      read_write_fast( O, fn=post_stratified_weights )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

      return(O)
    }


    if ("observation_weights" %in% todo) {
        # add region SA to results
        
        if (is.null(region)) stop("The parameter 'region' needs to be sent.")
        O = model_size_results(p=p, todo="post_stratified_weights" ) 
        pg = attr(O, "pg") 

        region_sa = switch( region,
          cfaall = "au_sa_km2",
          cfanorth = "cfanorth_surfacearea",
          cfasouth = "cfasouth_surfacearea",
          cfa4x = "cfa4x_surfacearea",
          cfa23 = "cfa23_surfacearea",
          cfa24 = "cfa24_surfacearea"
        )
        vn_keep = c("AUID", region_sa )
        pg = pg[, ..vn_keep]
        colnames(pg) = c("AUID", "SA")

        pg = pg[ O[,AUID], on="AUID" ] # bring in SA in correct sequence

        # finally, this is the post-stratified weight $\omega_i$ for sub-domain of focus
        O$post_stratified_weight = O$post_stratified_ratio * pg$SA  

        return(O)
    }


    if ( "observation_ratios_redo" %in% todo ){

        if (!exists("bioclass", p)) stop("Need to specify bioclass")
        
        message( "\nCreating post-stratification weights for: ", p$bioclass)

        # operate upon fit first and remove from meory to reduce RAM demand 
        # operate upon fit first and remove from meory to reduce RAM demand 
        fit = model_size_presence_absence( p=p, todo="load" ) 
        fit_summ_mean = fit$summary.fitted.values[["mean"]] 
        fit_summ_sd = fit$summary.fitted.values[["sd"]] 

        size_selectivity = fit$summary.random$'inla.group(cwd, method = "quantile", n = 11)' 
        # plot( (size_selectivity[,2]) ~ exp(size_selectivity[,1]))   log odds ratio
        # plot(exp(size_selectivity[,2])~ exp(size_selectivity[,1]))  # odds ratio
        # plot(1/exp(size_selectivity[,2])~ exp(size_selectivity[,1])) # selectivity ratio

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
    
        O$auid_prob_mean2 = M$individual_prob_mean[ ip2 ]
        O$auid_prob_sd2 = M$individual_prob_sd[ ip2 ]
        
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
        attr(O, "size_selectivity") = size_selectivity

        message( "\nSaving", observation_ratios )

        read_write_fast( O, fn=observation_ratios )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

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
        

        message( "\nSaving", observation_samples )
    
        read_write_fast( Osamples[iobs,], fn=observation_samples )  # on logit scale


        message( "\nSaving", prediction_samples )
        # same order as O, samples of the ratio of probabilities (a/i)
        post_stratified_ratio  = Osamples[ ip,] / Osamples[iobs,]  
        read_write_fast( post_stratified_ratio, fn=prediction_samples )  # on logit scale
        post_stratified_ratio = NULL

        message( "\nSaving", prediction_samples2 )
        post_stratified_ratio2 = Osamples[ ip2,] / Osamples[iobs,]   
        read_write_fast( post_stratified_ratio2, fn=prediction_samples2 )  # on logit scale


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
            Osamples[,i] = S[[i]]$latent[fss,]   # log odds ratio
        }
        
        Osamples = exp(Osamples )
        
        S = fss = NULL ;  gc()
    
        message( "\nSaving", bias_adjusted )

        read_write_fast( Osamples, fn=bias_adjusted )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave
 
        return(bias_adjusted)
    }


}