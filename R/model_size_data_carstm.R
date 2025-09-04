
model_size_data_carstm = function(p, redo=c(""), tokeep="necessary") { 
 
 
  if (is.null(p$bioclass)) stop("bioclass needs to be defined in p")

  p$selection$biologicals_using_snowcrab_filter_class = p$bioclass
  p$carstm_model_label = p$bioclass

  span = p$span(p$bioclass)

  fn = file.path( 
    p$project.outputdir, 
    paste( "size_distributions_tabulated_data_zeros_", p$bioclass, "_", paste0( span, collapse="_"), ".rdz", sep="" )  
  )


  Z = NULL 
  if ( ! "size_data" %in% redo ) {
    if (file.exists(fn)) {
      Z = read_write_fast( fn )
      if (tokeep=="necessary") {
        i = which(Z$tag == "observations")
        Z[ , uid := do.call(paste, .SD), .SDcols = c("AUID", "year", "cyclic", "cwd", "mat")]
        tokeep = unique( Z[i,uid])
        Z = Z[ uid %in% tokeep ,]
      }
      return(Z)
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
   # choose one:

  ### -------->>>  check cyclic matches dyear dyri
 
  # set level data from snow crab surveys 
  sppoly=areal_units( p=p )
  setcm = snowcrab.db( 
    p=p, 
    DS="carstm_inputs", 
    sppoly=sppoly, 
    redo=ifelse( "carstm_inputs" %in% redo, TRUE, FALSE )  
  )
  
  setDT(setcm)
  setnames(setcm, "id", "sid") 
  setnames(setcm, "mass", "mass_allspecies") 
  setnames(setcm, "len", "len_allspecies") 
  # length(unique(setcm$sid)) # 9375
  
  pred = setcm[tag=="predictions",] 
  setcm = setcm[tag=="observations",] 
  
  detcm = survey_db( DS="det", p=p  ) # size information, no, cm, kg
  setDT(detcm)
  setnames(detcm, "id", "sid") 
  setnames(detcm, "id2", "did") 

  isc = which( detcm$sid %in% unique( setcm$sid) )
  detcm = detcm[ isc, ]
  isc = NULL

  if (exists("selection", p)) {
    if (exists("biologicals", p$selection)) {  # filter biologicals
      isc = filter_data( detcm, p$selection$biologicals )
      if (length(isc) > 0) detcm = detcm[isc,]
      isc = NULL
    }

    if (exists("biologicals_using_snowcrab_filter_class", p$selection)) {  # filter biologicals using snow crab short-form ID
      warning( "Filtering using snow crab 'types' requires more data than is carried by survey_db. \n 
        .. Adding data directly from snowcrab.db .. this also means dropping other sources of data \n")
      det_sc = bio.snowcrab::snowcrab.db( DS ="det.georeferenced" )
      det_sc$spec = 2526
      det_sc$spec_bio =  taxonomy.recode( from="spec", to="parsimonious", tolookup=det_sc$spec ) # snow crab using groundfish codes
      det_sc$individual = paste( det_sc$trip, det_sc$set, det_sc$spec_bio, det_sc$crabno, sep=".")
      det_sc$filter.class = p$selection$biologicals_using_snowcrab_filter_class
      isc = bio.snowcrab::filter.class( x=det_sc, type=p$selection$biologicals_using_snowcrab_filter_class )
      if (length(isc) > 0) det_sc = det_sc[isc, c("individual", "filter.class") ]
      isc = NULL
      setDT( det_sc )
      detcm$individual = paste( detcm$did, detcm$individual, sep=".")
      detcm = det_sc[ detcm, on="individual" ]          
      detcm = detcm[ !is.na(detcm$filter.class), ]
    }
  }

  # --- NOTE det was not always determined and so totals from det mass != totals from cat nor set for all years
  # cf_det is the weight to make it sum up to the correct total catch (vs any subsamples) and tow length, etc
  
  detcm = detcm[ sex %in% c("0", "1") , ]  
  detcm = detcm[ mat %in% c("0", "1") , ]
  detcm$cw = detcm$len * 10  # mm -> cm
  detcm$logcw = log(detcm$cw)

  detcm = detcm[, .(sid, mat, cw, logcw, mass, cf_det_no)]  # mass in kg

  # trim a few strange data points
  o = lm( log(mass) ~ logcw , detcm)
  todrop = which(abs(o$residuals) > 0.5)
  if (length(todrop)>0) detcm = detcm[-todrop,]

  # internally on log 
  lspan = span
  lspan[1]  = log(span[1]) 
  lspan[2]  = log(span[2]) 

  todrop = detcm[ logcw < lspan[1], which=TRUE ] 
  if (length(todrop)>0) detcm = detcm[-todrop,]

  todrop = detcm[ logcw > lspan[2], which=TRUE ] 
  if (length(todrop)>0) detcm = detcm[-todrop,]
  
  # additional QA/QC 

  if (p$bioclass %in% c("female", "f.imm", "f.mat" )) {

    todrop = detcm[ logcw > log(80) & mat=="0", which=TRUE ] 
    if (length(todrop)>0) detcm = detcm[-todrop,]

    todrop = detcm[ logcw < log(35) & mat=="1", which=TRUE ] 
    if (length(todrop)>0) detcm = detcm[-todrop,]
  
  }

  if (p$bioclass %in% c("male", "m.imm", "m.mat" )) {

    todrop = detcm[ logcw > log(135) & mat=="0", which=TRUE ] 
    if (length(todrop)>0) detcm = detcm[-todrop,]

    todrop = detcm[ logcw < log(49) & mat=="1", which=TRUE ] 
    if (length(todrop)>0) detcm = detcm[-todrop,]

  }

  detcm$cwd = discretize_data( detcm$logcw, span=lspan  )  # this will truncate sizes

  detcm = detcm[ is.finite(cwd) ,]
  detcm$N = 1
  mats = unique(detcm$mat)

  # length(unique(det$sid)) # 8258

  Z = setcm[ detcm, on=.(sid)] # 569657
  Z$zid = paste( Z$sid, Z$mat, Z$cwd, sep="_" )

  # zeros at sampling locations: 
  #   .. CJ required to get zero counts dim(N) # 1237368 
  #   .. this is an approximationas the number of zeros are also controlled by size resoltuion/span
  Z0 = CJ( 
    mat=c("0", "1"), 
    cwd = discretize_data(span=lspan),  # midpoints
    sid=setcm$sid, 
    unique=TRUE 
  )
  Z0$mass = NA
  Z0$logcw = Z0$cwd  # fill with midpoints 
  Z0$cw = exp(Z0$logcw)  #  midpoint on original scale ... and a dummy to allow merges with observations
  Z0$N = 0

  Z0 = setcm[,.(sid, data_offset)][ Z0, on=.(sid)] 
  
  Z0$cf_det_no = 1 / Z0$data_offset  # this is actual cf_set_no ... as no is 0, assume set swept area as multiplier
  Z0$data_offset = NULL

  Z0 = setcm[ Z0, on=.( sid ) ]  # n=2690336
  Z0$zid = paste( Z0$sid, Z0$mat, Z0$cwd, sep="_" )

  withdata = unique(Z$zid)
  toremove = unique(Z0[ zid %in% withdata, which=TRUE])
  Z0 = Z0[-toremove, ]  

  Z = rbind(Z, Z0) # 1683904
  # length(unique(Z$sid)) # 9374
  # sum(Z$N, na.rm=T) # n=539181
  Z$zid = NULL
  Z0 = NULL 
  gc()

  Z = Z[ year %in%  p$yrs, ] # n=494 541
  
  # prediction surface
  pred$pid = paste( pred$AUID, pred$year, pred$dyri, sep="_")

  # to capture zeros 
  P0 = CJ( 
    mat=mats, 
    cwd=discretize_data( span=lspan ), 
    pid=pred$pid, 
    mass=NA,
    N =NA,
    cf_det_no = 1,
    unique=TRUE 
  )

  P0$logcw = P0$cwd

  pred = P0[ pred, on=.(pid), allow.cartesian=TRUE]
  pred$cw = exp(pred$logcw)
  setnames(pred, "pid", "sid" )
  
  Zvn = names(Z)

  Z = rbind(Z, pred[, ..Zvn ])

  # a final round of data trimming based upon size, sex, maturity and time of year

  if (p$bioclass %in% c("female", "f.imm", "f.mat" )) {
      
    todrop = Z[ logcw > log(95), which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]

    # todrop = Z[ logcw > log(80) & mat=="0", which=TRUE ] 
    # if (length(todrop)>0) Z = Z[-todrop,]

    # todrop = Z[ logcw < log(35) & mat=="1", which=TRUE ] 
    # if (length(todrop)>0) Z = Z[-todrop,]
  
  }


  if (p$bioclass %in% c("male", "m.imm", "m.mat" )) {

    todrop = Z[ logcw > log(155), which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]

    # todrop = Z[ logcw > log(135) & mat=="0", which=TRUE ] 
    # if (length(todrop)>0) Z = Z[-todrop,]

    # todrop = Z[ logcw < log(49) & mat=="1", which=TRUE ] 
    # if (length(todrop)>0) Z = Z[-todrop,]
  }
  
  # match prediction range to observation range for season (cyclic)
  data_dyears = range(Z[tag=="observations", "cyclic"] ) 
  # todrop = which(Z$cyclic < data_dyears[1] | Z$cyclic > data_dyears[2] )
  # if (length(todrop)>0) Z = Z[-todrop,]

  Z$cyclic = match( Z$dyri, p$cyclic_levels ) 
  Z$cyclic_space = Z$cyclic # copy cyclic for space - cyclic component .. for groups, must be numeric index


  Z$year = as.factor(Z$year)
  Z$region = as.factor(Z$region)
  Z$cwd = as.numeric(as.character(Z$cwd))
  
  Z$data_offset = 1/Z$cf_det_no

  # override -- "pa" was from set level determination, here we are using individual-level bernoulli 0/1 process
  setnames(Z, "pa", "pa_set")
  setnames(Z, "N",  "pa")  

  Z$mat =  as.numeric(as.character(Z$mat))
  Z$mat_group = Z$mat  # copy needed for INLA (to use a a group)

  Z$cwd2 = Z$cwd
  Z$time_cw = Z$time


  attr(Z, "bioclass") = p$bioclass
  attr(Z, "span") = span  
  attr(Z, "data_dyears") = data_dyears
  attr(Z, "yrs") = p$yrs
  attr(Z, "sppoly") = sppoly
  
  read_write_fast( data=Z, fn=fn )  # save full prediction surface .. subset on return (above)


  return("complete")
}
