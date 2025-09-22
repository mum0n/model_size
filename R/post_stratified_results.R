post_stratified_results = function(p, todo="load"){

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
    
    attr(O, "size_selectivity") = size_selectivity

    read_write_fast( O, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(M)
}
