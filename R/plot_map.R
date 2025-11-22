plot_map = function( H, spatial_domain="snowcrab", plottype="points" ) {

  library(ggplot2)
  p = spatial_parameters(p=list(spatial_domain=spatial_domain))
  plt_crs = st_crs(p$aegis_proj4string_planar_km)
  cst = coastline_db( p=p, project_to=plt_crs ) 
  isob = isobath_db( DS="isobath", depths=c(100, 200, 300), project_to=plt_crs)
  isob$level = as.factor( isob$level)

  if (plottype == "points" ) {   
    plt = ggplot() +
        geom_sf( data=cst, show.legend=FALSE ) +
        geom_sf( data=isob, aes( alpha=0.1, fill=level), lwd=0.1, show.legend=FALSE) +
        geom_point(data=H, aes(x=plon, y=plat, colour=pt), size=5) +
        coord_sf(xlim = c(270, 940 ), ylim = c(4780, 5200 )) +
        theme(legend.position="inside", legend.position.inside=c(0.08, 0.8), 
            axis.title.x=element_blank(),
            axis.title.y=element_blank())
  }
  
  return(plt)
}
