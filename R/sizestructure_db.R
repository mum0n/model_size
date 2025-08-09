sizestructure_db = function( p=NULL, DS=NULL, datasource="snowcrab",   
    rawdatadir=file.path(homedir, "projects", "model_size", "outputs", "size_structure", "modes_kernel_mixture_models_set"), 
    regions=c("cfanorth", "cfasouth", "cfa4x"), redo=FALSE, ... ) {
 
    if ( DS=="carstm_inputs") {
        if (datasource=="snowcrab") {
            km = NULL
            for (yr in p$yrs) {
                fn = file.path( rawdatadir, paste("kmm_parameter_summaries_", yr, ".rdz", sep="") )
                if (file.exists(fn)) {
                    km = rbind(km, aegis::read_write_fast(fn) )
                }
            }
            setDT(km)
            pg = areal_units( p=p , areal_units_directory=p$project.outputdir )  # reload
            M = bio.snowcrab::snowcrab.db(p=p, DS="carstm_inputs",  sppoly=pg, redo=FALSE )
            setDT(M)

            set = M[ tag=="observations", ]
            set$sid = gsub("\\.", "~", set$id) 

            # # sex codes
            # male = 0
            # female = 1
            # sex.unknown = 2

            # # maturity codes
            # immature = 0
            # mature = 1
            # mat.unknown = 2

            det = bio.snowcrab::snowcrab.db(p=p, DS="det.georeferenced" )
            setDT(det)
            det$sid = paste( det$trip, det$set, sep="~")
            dsex = det$sex
            det$sex = NA
            det[ dsex==0, "sex" ] = "m"
            det[ dsex==1, "sex" ] = "f"

            dmat = det$mat
            det$mat = NA
            det[ dmat==0, "mat" ] = "i"
            det[ dmat==1, "mat" ] = "m"

            det$yr = as.character(det$yr)

            sk = list()
            for (vn in unique(km$stage)) {
                sk[[vn]] = km[ stage==vn, ][ set, on="sid"]  # same order as m observations
                # merge creates na's
                sk[[vn]]$stage = vn
           
                sksex = sort(unique(sk[[vn]]$sex)) 
                skmat = sort(unique(sk[[vn]]$mat))
                skinstar = sort(unique(sk[[vn]]$instar))
           
                sk[[vn]]$sex =  sksex # sort removes na
                sk[[vn]]$mat = skmat
                sk[[vn]]$instar = skinstar
                sk[[vn]][ is.na(sigmasq_mean), "sigmasq_mean"] = 0 
#                sk[[vn]][ is.na(sigmasq_sd), "sigmasq_sd"] = 0
                sk[[vn]][ is.na(alpha_mean), "alpha_mean"] = 0 
#                sk[[vn]][ is.na(alpha_sd), "alpha_sd"] = 0
                inds = det[ 
                    sid %in% sk[[vn]][["sid"]] & 
                    mat==skmat & 
                    sex==sksex , ]
                ninds = inds[, .(Nkmm=.N), by=c("sid", "yr")]
                sk[[vn]] = ninds[ sk[[vn]], on=c("sid", "yr") ]
            }

            return( list(M=M[tag=="predictions",], km=km, pg=pg, set=set, sk=sk ) )
        }

    } 
    
    if ( DS=="arealunits") {
        fn = file.path(p$project.outputdir,  "areal_units_size_structure.rdz")
        pg = NULL
        if (!redo) {
            pg = read_write_fast( fn )
            return(pg)
        }
        pg = areal_units( p = p, areal_units_directory=p$project.outputdir )
        read_write_fast( data=pg, fn=fn)
        return(pg)
    }

    if ( DS=="get_rawdata") {
    
    }
}
