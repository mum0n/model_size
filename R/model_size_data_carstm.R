
model_size_data_carstm = function(p, redo=c(""), zero_weight="unit" ) { 

  # p$selection$biologicals_using_snowcrab_filter_class = p$bioclass

  fn = file.path( 
    p$modeldir,  
    paste( "size_distributions_tabulated_data_zeros.rdz", sep="" )  
  )
  


  Z = NULL 
  if ( ! "size_data" %in% redo ) {
    if (file.exists(fn)) {
      Z = read_write_fast( fn )

      if (!exists("bioclass", p)) return(Z)  
      if (p$bioclass == "all") return(Z)  

      i = filter.class( Z, p$bioclass )
      if (length(i) > 0) Z = Z[i, ]

      if (zero_weight=="unit") {
        k = which(Z$tag=="observations" & Z$crabno==0 )
        if (length(k)>0) Z$data_offset[k] = 1  # override swept area with unit of prediction (1km^2) # this increases the importance of zero-values
      }

      Z$mat = as.character(Z$mat)
      Z$sex = as.character(Z$sex)
      
      # subset data 
      # essentially drop all predict times except those required to reduce computations
      ukuid_obs  = unique(Z[ tag=="observations", kuid])  # add all observations
      ukuid_pred = unique(Z[ cyclic==p$prediction_dyear_index, kuid ])  # and add time-slice for prediction period (1 sept)
      
      ukuid = unique( c(ukuid_obs, ukuid_pred) )
      
      Z = Z[ kuid %in% ukuid ,]

      # additional copies of data for INLA, refer to model formula
      Z$space2 = Z$space
      Z$cwd2 = Z$cwd

      Z$space_time = Z$space
      Z$time_space = Z$time


      Z$cwd3 = Z$cwd
      Z$time3 = Z$time
      # Z$time_cw = Z$time

      # Z$cyclic_space = Z$cyclic # copy cyclic for space - cyclic component .. for groups, must be numeric index


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
 
  # if (is.null(p$bioclass)) stop("bioclass needs to be defined in p")

  p$selection$biologicals_using_snowcrab_filter_class = "all"
  p$carstm_model_label = "all"
 
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
  setcm$pa = NULL

  # length(unique(setcm$sid)) # 9375
  
  pred = setcm[tag=="predictions",] 
  setcm = setcm[tag=="observations",] 
  
  detcm = survey_db( DS="det", p=p  ) # size information, no, cm, kg
  setDT(detcm)
  setnames(detcm, "id", "sid") 
  setnames(detcm, "id2", "did") 

  # dim(setcm) [1] 9374   38
  # dim(detcm) [1] 3452897      19

  # table(detcm$data.source)
  # groundfish   snowcrab 
  # 2669078      579052 

  isc = which( detcm$sid %in% unique( setcm$sid) )
  detcm = detcm[ isc, ]
  isc = NULL

  # dim(detcm)  # [1] 543396     19

  if (exists("selection", p)) {
    if (exists("biologicals", p$selection)) {  # filter biologicals on size
      isc = filter_data( detcm, p$selection$biologicals )
      if (length(isc) > 0) detcm = detcm[isc,]
      isc = NULL
    }

    if (exists("biologicals_using_snowcrab_filter_class", p$selection)) {  # filter biologicals using snow crab short-form ID
      # warning( "Filtering using snow crab 'types' requires more data than is carried by survey_db. \n 
      #  .. Adding data directly from snowcrab.db .. this also means dropping other sources of data \n")
      det_sc = bio.snowcrab::snowcrab.db( DS ="det.georeferenced" )
      det_sc$spec = 2526
      det_sc$spec_bio =  taxonomy.recode( from="spec", to="parsimonious", tolookup=det_sc$spec ) # snow crab using groundfish codes
      det_sc$individual = paste( det_sc$trip, det_sc$set, det_sc$spec_bio, det_sc$crabno, sep=".")
      det_sc$filter.class = p$selection$biologicals_using_snowcrab_filter_class
      isc = bio.snowcrab::filter.class( x=det_sc, type=p$selection$biologicals_using_snowcrab_filter_class )
      if (length(isc) > 0) det_sc = det_sc[isc, c("individual", "crabno",  "filter.class") ]
      isc = NULL
      setDT( det_sc )
      detcm$individual = paste( detcm$did, detcm$individual, sep=".")
      detcm = det_sc[ detcm, on="individual" ]          
      detcm = detcm[ !is.na(detcm$filter.class), ]
    }
  }

  # --- NOTE det was not always determined and so totals from det mass != totals from cat nor set for all years
  # cf_det is the weight to make it sum up to the correct total catch (vs any subsamples) and tow length, etc
 
  # table(detcm$sex, detcm$mat) # 541537
  #        sex 0  sex 1
  # mat 0 240393 104653
  # mat 1  98691  97800

  detcm = detcm[ sex %in% c(0, 1) , ]  
  detcm = detcm[ mat %in% c(0, 1) , ]

  detcm$cw = detcm$len * 10  # mm -> cm
  detcm$logcw = log(detcm$cw)

  detcm = detcm[, .(sid, crabno, sex, mat, cw, logcw, mass, cf_det_no)]  # mass in kg

  
  # trim a few strange data points
  o = lm( log(mass) ~ as.factor(mat) + as.factor(sex) + logcw:as.factor(mat) + logcw:as.factor(sex), detcm)
  todrop = which(abs(o$residuals) > 0.5)
  if (length(todrop)>0) detcm = detcm[-todrop,]

  # internally on log 
  lspan = c( log(p$size_range("all")), p$size_bandwidth) 

  lbrks = discretize_data(span=lspan)
  
  # dim(detcm) # [1] 541533


  # additional QA/QC 

  detcm$cwd = discretize_data( detcm$logcw, span=lspan  )  # this will truncate sizes

  detcm = detcm[ is.finite(cwd) ,]
  detcm$pa = 1
  
  # length(unique(detcm$sid)) # 7861

  Z = setcm[ detcm, on=.(sid)] # 541533
  Z$zid = paste( Z$sid, Z$sex, Z$mat, Z$cwd, sep="_" )

  # recover zero counts at sampling locations with no observations: 
  #   .. CJ required to get zero counts # 2812200 
  #   .. this is an approximation as the number of zeros are also controlled by size resoltuion/span
  Z0 = CJ( 
    sex=c(0, 1), 
    mat=c(0, 1), 
    cwd = lbrks,  # midpoints
    sid=setcm$sid, 
    unique=TRUE 
  )
  Z0$mass = NA
  Z0$logcw = Z0$cwd  # fill with midpoints 
  Z0$cw = exp(Z0$logcw)  #  midpoint on original scale ... and a dummy to allow merges with observations
  Z0$pa = 0
	
  # dim(Z0) # 1124880

  Z0 = setcm[,.(sid, data_offset)][ Z0, on=.(sid)] 

  Z0$cf_det_no = 1 / Z0$data_offset  # this is actual cf_set_no ... as no is 0, assume set swept area as multiplier
  Z0$data_offset = NULL

  Z0 = setcm[ Z0, on=.( sid ) ]  # n=2812200
  Z0$zid = paste( Z0$sid, Z0$sex, Z0$mat, Z0$cwd, sep="_" )

  withdata = unique(Z$zid)
  toremove = unique(Z0[ zid %in% withdata, which=TRUE])
  if (length(toremove) > 0) Z0 = Z0[-toremove, ]  
	Z0$crabno = 0

  Z = rbind(Z, Z0) # 1565641
  # length(unique(Z$sid)) # 9374
  # sum(Z$pa, na.rm=T) # 541533
  Z$zid = NULL
  Z0 = NULL 
  gc()

  Z = Z[ year %in%  p$yrs, ] # n=3189976
  
  # prediction surface
  pred$pid = paste( pred$AUID, pred$year, pred$dyri, sep="_")

  # to capture zeros 
  P0 = CJ( 
    sex=c(0, 1), 
    mat=c(0, 1), 
    cwd=lbrks, 
    pid=pred$pid, 
    mass=NA,
    pa =NA,
    cf_det_no = 1,
    unique=TRUE 
  )

  P0$logcw = P0$cwd

  pred = P0[ pred, on=.(pid), allow.cartesian=TRUE]  # 29090880
  pred$cw = exp(pred$logcw)
  setnames(pred, "pid", "sid" )
	pred$crabno = 0
  
  vnkeep = c(
    "sid", "AUID", "crabno",  "tag", "data_offset", 
    "z", "substrate.grainsize", "t", "pca1", "pca2",
    "year", "dyear", "dyri", "space", "time",  
    "sex", "mat", "cw", "logcw", "cwd", "mass", "cf_det_no",  "pa" 
  )

  Z = rbind(Z[, ..vnkeep], pred[, ..vnkeep ])  # 30,656,521

  # size truncation to range of reliable data
  bc = "f.imm"
    sizerange = log(p$size_range( bc ))
    todrop = Z[ sex==1 & mat==0 & cwd < sizerange[1], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
    todrop = Z[ sex==1 & mat==0 & cwd > sizerange[2], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
  
  bc = "f.mat"
    sizerange = log(p$size_range( bc ))
    todrop = Z[ sex==1 & mat==1 & cwd < sizerange[1], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
    todrop = Z[ sex==1 & mat==1 & cwd > sizerange[2], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
  
  bc = "m.imm"
    sizerange = log(p$size_range( bc ))
    todrop = Z[ sex==0 & mat==0 & cwd < sizerange[1], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
    todrop = Z[ sex==0 & mat==0 & cwd > sizerange[2], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
  
  bc = "m.mat"
    sizerange = log(p$size_range( bc ))
    todrop = Z[ sex==0 & mat==1 & cwd < sizerange[1], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]
    todrop = Z[ sex==0 & mat==1 & cwd > sizerange[2], which=TRUE ] 
    if (length(todrop)>0) Z = Z[-todrop,]

  Z$cyclic = match( Z$dyri, p$cyclic_levels ) 
  Z$year = as.factor(Z$year)
  
  Z$data_offset = 1/Z$cf_det_no
  Z$cf_det_no = NULL 
  
  key_vars = c("AUID", "year", "cyclic", "cwd", "mat", "sex" )

  Z[ , kuid := do.call(paste, .SD), .SDcols = key_vars ]  

  # match prediction range to observation range for season (cyclic)
  # data_dyears = range(Z[tag=="observations", "cyclic"] ) 
  # todrop = which(Z$cyclic < data_dyears[1] | Z$cyclic > data_dyears[2] )
  # if (length(todrop)>0) Z = Z[-todrop,]

  # override -- "pa" was from set level determination, here we are using individual-level bernoulli 0/1 process
 
  # reduce file size .. most of these are unused or can be recovered from data

  attr(Z, "brks") = lbrks  
  attr(Z, "yrs") = p$yrs
  attr(Z, "sppoly") = sppoly
  
  # attr(Z, "data_dyears") = data_dyears
  # dim(Z) # 45701417       30
  # table(Z[tag=="observations", .(sex, mat) ] )  # includes zeros
  #            mat
  # sex      0      1
  #   0 845639 294837
  #   1 637057 287564

  # table(detcm$sex, detcm$mat)  # excluding zeros
  #        sex 0  sex 1
  # mat 0 240393 104653
  # mat 1  98691  97800

  message("Saving data file: ", fn )
  
  read_write_fast( data=Z, fn=fn )  # save full prediction surface .. subset on return (above)

  return("complete")
}
