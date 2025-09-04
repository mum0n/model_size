surface_area_estimates = function(p, pg, region="cfaall",
  sa_vars = c("cfanorth_surfacearea", "cfasouth_surfacearea", "cfa23_surfacearea", "cfa24_surfacearea", "cfa4x_surfacearea", "au_sa_km2")) {

  for (i in sa_vars) {
      attr(pg[,i], "units") = NULL
      class(pg[,i]) = NULL
  }
  tokeep = c("AUID", sa_vars, "strata_to_keep" )
  pg = data.table(pg)[, ..tokeep]

  attr(pg[["au_sa_km2"]], "units") = NULL
  class(pg[["au_sa_km2"]]) = NULL

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

  return(pg)
}
