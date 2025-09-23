post_stratified_results = function(p, todo="load", region=NULL) {

    fn = file.path( p$modeldir, "post_stratified_results.rdz" ) 

    if ( "load" %in% todo ){
        message( "\nLoading from file: ", fn)
        if (file.exists(fn)) {
            O = read_write_fast(fn)
            if (!is.null(O)) {
                return(O)
            }
        }
    }

    if ("size_selectivity" %in% todo) {
        
        O = post_stratified_results(p=p, todo="load" ) 
        size_selectivity = attr(O, "size_selectivity")
        if (!is.null(region)) {
          size_selectivity = size_selectivity[[region]]
        }  
        return( size_selectivity )
    }


    if ("area_weighted" %in% todo) {
        
        if (is.null(region)) stop("The parameter 'region' needs to be sent.")
        O = post_stratified_results(p=p, todo="load" ) 
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


    # default is to "redo" if todo is not captured above
    O = NULL
    size_selectivity = list()

    for (bc in c("f.imm", "f.mat", "m.imm", "m.mat") ){
        print(bc)
				p$bioclass = bc
        o = NULL
        o = post_stratified_predictions( p=p, todo="load" )
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

    # add surface area information as an attribute
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

    attr(O, "pg") = pg

    attr(O, "size_selectivity") = size_selectivity

    read_write_fast( O, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(O)
}
