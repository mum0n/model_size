 
size_data_tabulated = function(
    p=p, 
    outdir=NULL, 
    toget="",
    regions=c("cfanorth", "cfasouth", "cfa4x"), 
    span = NULL,
    redo=FALSE,
    add_zeros=FALSE,
    pg=NULL,
    Y=NULL ) { 

    # tabulations, computed directly and via model-form


    if (is.null(outdir)) outdir = file.path(p$project.outputdir, "size_structure") 
    if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings =FALSE) 
      

    # note ranges in CW will be log transformed later
    if (is.null(span)) {
        span = function( sexid) {
            switch(sexid,
                male   = c( 5, 155, 40),
                female = c( 5, 95,  40)
            )
        }
    }

    # sex codes
    # male = 0
    # female = 1
    # sex.unknown = 2

    # # maturity codes
    # immature = 0
    # mature = 1
    # mat.unknown = 2
 
  
    if (is.null(Y)) {
        if (exists("yrs", p)) Y = p$yrs
    }


    if (toget=="base_data") {
        # georeferenced "crude" data
        savedir = file.path(outdir, "base_data")
        if (!dir.exists(savedir)) dir.create(savedir, recursive=TRUE, showWarnings =FALSE) 

        if (!redo) {
            M = NULL
            if (!is.vector(Y)) stop("Y should be a year vector")
            for (yr in as.character(Y)) {
                fn = file.path( savedir, paste("size_distributions_base_data_", yr, ".rdz", sep="" ))
                if (file.exists(fn)) {
                    m = NULL
                    m = aegis::read_write_fast(fn)
                    m$year = yr
                    M = rbind( M, m)
                }
            }
            return(M)
        }
 
        # add georeference to pg

        set = snowcrab.db( DS="set.clean")
        setDT(set)
        set$sid = paste(set$trip, set$set, sep="~")
        set = set[, .(sid, lon, lat)]
        set$space_id = NA
        Z = sf::st_as_sf( set[,.(lon, lat)], coords=c("lon", "lat") )
        st_crs(Z) = st_crs( projection_proj4string("lonlat_wgs84") )
        for (aoi in 1:nrow(pg)) {
            ks = which(!is.na( st_points_in_polygons(pts=Z, polys=pg[aoi, "AUID"], varname= "AUID" ) ))
            if (length(ks) > 0 ) set$space_id[ks] = pg$AUID[aoi]
        }
        set = set[, .(sid, space_id)]

        for (yr in Y) {
            S = copy(set)
            basedata = size_data_tabulated(p=p, toget="tabulated_data", outdir=outdir, Y=yr)
            basedata = S[ basedata, on=.(sid)]
            fn = file.path( savedir, paste("size_distributions_base_data_", yr, ".rdz", sep="" ))
            read_write_fast( data=basedata, fn=fn )
            print(fn)
        }

        return(size_data_tabulated(p=p, toget="base_data", span=span, Y=Y, outdir=outdir, redo=FALSE))
    }


    if (toget=="tabulated_data") {
        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        fn = file.path( outdir, "size_distributions_tabulated_data.rdz" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) {
                M =aegis::read_write_fast(fn)
                if (!is.null(Y)) M = M[ year %in% Y, ]
                if (add_zeros) {
                    # merge zeros here so we do not have to store the massive intermediary file
                    # CJ required to get zero counts dim(N) # 171624960    
                    M = M[ CJ( sex, mat, cwd, sid, unique=TRUE ), 
                         on=.( sex, mat, cwd, sid ) ]
                }
                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
                return(M)
            }
        }

        M = size_data_tabulated(p=p, toget="base_data", span=span, outdir=outdir )

        # aggregate by cwd 
        M = M[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, sex, mat, cwd, sid) ]
        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$cwd = as.factor(M$cwd)
        read_write_fast( data=M, fn=fn )
        # return this way to add zeros, if required
        return( size_data_tabulated(p=p, toget="tabulated_data", outdir=outdir, add_zeros=add_zeros, redo=FALSE ) )
    }



    if (toget=="tabulated_data_by_stage") {
        # NOTE: sampling event = "sid"
        # NOTE: size = "cwd"  
        fn = file.path( outdir, "size_distributions_tabulated_data_by_stage.rdz" )
        if (!redo) {
            M = NULL
            if (file.exists(fn)) {
                M =aegis::read_write_fast(fn)
                if (!is.null(Y)) M = M[ year %in% Y, ]
                if (add_zeros) {
                    # merge zeros here so we do not have to store the massive intermediary file
                    # CJ required to get zero counts dim(N) # 171624960    
                    M = M[ CJ( stage, sid, unique=TRUE ), 
                        on=.(  stage, sid ) ]
                }
                M[ !is.finite(N),   "N"] = 0
                M[ !is.finite(mass), "mass"] = 0 
                M[ !is.finite(sa), "sa"] = 1 #dummy value
                M$density = M$N / M$sa
                M[ !is.finite(density), "density"] = 0  
                return(M)
            }
        }

        M = size_data_tabulated(p=p, toget="base_data", span=span, outdir=outdir )
        # aggregate by cwd 

        mds = modal_distribution_data(p=p, toget="modal_groups", outdir=outdir, redo=FALSE )
        
        M$stage = filter.stage( M, mds ) 
        M = M[!is.na(stage),]

        M = M[,  .( N=.N, mass=mean(mass, na.rm=TRUE), sa=mean(sa, na.rm=TRUE) ),  
            by=.( region, year, stage, sid) ]
        M$year = as.factor(M$year)
        M$region = as.factor(M$region)
        M$stage = as.factor(M$stage)
        read_write_fast( data=M, fn=fn )
        # return this way to add zeros, if required
        return( size_data_tabulated(p=p, toget="tabulated_data_by_stage", outdir=outdir, add_zeros=add_zeros, redo=FALSE ) )
    }

 
    if (toget=="linear_model") {
        require(biglm)
        fn = file.path( outdir, "size_distributions_lm.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_data_tabulated(p=p, toget="tabulated_data", span=span, outdir=outdir, add_zeros=TRUE )
        setDT(M)
        fit = biglm::bigglm( density ~ region:year:mat:cwd:sex - 1, data=M, family=gaussian(link="identity") )
        O = summary(fit)$coefficients
        res = tstrsplit(rownames(O), ":")
        setDT(res)
        setnames(res, new=c("region", "year", "mat", "cwd", "sex"))
        res[,region:=as.numeric(gsub( "region", "", region))]
        res[,year:=as.factor(gsub( "year", "", year))]
        res[,mat:=as.factor(gsub( "mat", "", mat))]
        res[,cwd:=as.numeric(gsub( "cwd", "", cwd))]
        res[,sex:=as.numeric(gsub( "sex", "", sex))]
        O = data.table(O)
        colnames(O) = c("density",  "density_se", "t", "p")
        O = cbind( res, O )
        read_write_fast( data=O, fn=fn )
        return(O)
    }


    if (toget=="poisson_glm") {
        stop("This takes far too long to use ... ")
        fn = file.path( outdir, "size_distributions_poisson_glm.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_data_tabulated(p=p, toget="tabulated_data", span=span, outdir=outdir, add_zeros=TRUE )
        M$ID = as.factor( paste( M$region, M$year, M$sex, M$mat, M$cwd, sep="_") )
        # subset 
        ss = M[ region=="cfanorth" & sex=="0" & year %in% as.character(2015:2022), which=TRUE]
        # NOTE this is too large of a problem for glm
        fit = biglm::bigglm( N ~ ID - 1 +offset(log_sa), data=M[ss,], 
            family=poisson(link="log"), na.action="na.omit" )
        P = data.table( ID = names(coef(fit)), mean=coef(fit) )
        nm = matrix( unlist( strsplit(P$ID, "_")), ncol=5, byrow=TRUE)
        P = cbind( P, nm)
        names(P) = c("ID", "N", "region", "year", "sex", "mat", "cwd")
        P$region = gsub("^ID", "", P$region)
        O = list( fit=fit, P=P )
        read_write_fast( data=O, fn=fn )
        return(O)
    }


    if (toget=="poisson_inla") {
        stop("This takes far too long to use ... ")
        require(INLA)
        fn = file.path( outdir, "size_distributions_inla_poisson.rdz" )
        O = NULL
        if (!redo) {
            if (file.exists(fn)) O =aegis::read_write_fast(fn)
            return(O)
        }
        M = size_data_tabulated(p=p, toget="tabulated_data", span=span, outdir=outdir, add_zeros=TRUE )
        M$tag ="o"
        P = CJ( 
            N = NA,
            log_sa = 0,   # log(1) ... 1km^2
            region = regions,
            year = p$yrs,  
            sex = c("0", "1"), 
            mat = c("0", "1"),
            cwd = levels(M$cwd),
            tag= "p"
        )
        pn = names(P)
        M = copy( rbind( M[, ..pn ], P ) ) # a deep copy and not a reference
        # M$ID = paste( M$region, M$year, M$sex, M$mat, M$cwd, sep="_")
        # subset 
        ss = M[ region=="cfanorth" & sex=="0" & year %in% as.character(2015:2022), which=TRUE]
        fit = inla( N ~ year + mat + cwd + year:mat:cwd - 1 + offset(log_sa), data=M[ss,], 
            family="poisson", verbose=TRUE)
        iP = M[tag=="p", which=TRUE]
        P$N = fit$summary.fitted.values$mean[iP]
        P$Nsd = fit$summary.fitted.values$sd[iP]
        P$Nlb = fit$summary.fitted.values$"0.025quant"[iP]
        P$Nub = fit$summary.fitted.values$"0.975quant"[iP]
        O = list( fit=fit, P=P )
        read_write_fast( data=O, fn=fn )
        return(O)
    }

}

 

 
