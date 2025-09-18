post_stratified_results = function(p, todo="load"){

    fn = file.path( p$modeldir, "post_stratified_results.rdz" ) 

    if ( "load" %in% todo ){
        message( "\nLoading from file: ", fnout)
        if (file.exists(fnout)) {
            O = read_write_fast(fnout)
            if (!is.null(O)) {
                return(O)
            }
        }
    }

    p$bioclass = "all"
    M = model_size_data_carstm( p=p )  

    O = NULL
    size_selectivity = list()

    for (bc in c("f.imm", "f.mat", "m.imm", "m.mat") ){
        p$bioclass = bioclass
        o = NULL
        o = post_stratified_predictions( p=p, todo="load" )
        ss = size_selectivity[[bc]] = attr(o, "size_selectivity") 
            
        # plot( (ss[,2]) ~ exp(ss[,1]))   log odds ratio
        # plot(exp(ss[,2])~ exp(ss[,1]))  # odds ratio
        # plot(1/exp(ss[,2])~ exp(ss[,1])) # selectivity ratio

        o$odds_ratio = splinefun(
            x = exp(ss[,1]), 
            y = exp(ss[,2]), 
            xout=o$cw 
        )

        o = o[, .(
            kuid, 
            individual_prob_mean,
            individual_prob_sd,
            auid_prob_mean,
            auid_prob_sd,
            odds_ratio
        )]

        O = rbind( O, o )
    }

    M = O[ M, on=.(kuid) ]
    O = NULL

    M$post_stratified_ratio  = M$auid_prob_mean / M$individual_prob_mean   
    M$post_stratified_ratio2 = M$auid_prob_mean2 / M$individual_prob_mean   
    
    attr(M, "size_selectivity") = size_selectivity

    read_write_fast( M, fn=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(M)
}
