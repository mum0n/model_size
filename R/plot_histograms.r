
plot_histograms = function( M, REG, YR, 
        plotsavedir=NULL, 
        fn_prefix="size_freq_naive_sum", nbins=40, plottype="histogram", 
        gg_labs=labs(x="cwd", y="density; n/km^2"),
        show_error = FALSE
    ) {


  if (plottype =="histogram") {
    pmm = ggplot( M[ region %in% REG & year==YR & sex==0 & mat==1, ], aes(cw) ) +  
        geom_histogram(bins=nbins, fill="lightgreen", col="grey") + 
        xlim(c(0, 145))
    pfm = ggplot( M[ region %in% REG & year==YR & sex==1 & mat==1, ], aes(cw) ) +  
        geom_histogram(bins=nbins, fill="lightgreen", col="grey") + 
        xlim(c(0, 85))
    pmi = ggplot( M[ region %in% REG & year==YR & sex==0 & mat==0, ], aes(cw) ) +  
        geom_histogram(bins=nbins, fill="lightgreen", col="grey") + 
        xlim(c(0, 145))
    pfi = ggplot( M[ region %in% REG & year==YR & sex==1 & mat==0, ], aes(cw) ) +  
        geom_histogram(bins=nbins, fill="lightgreen", col="grey") + 
        xlim(c(0, 85))
  }


  if (plottype =="barchart") {
      pmm = ggplot( M[ region %in% REG & sex==0 & mat==1, ], aes(cwd, mn) ) + 
        geom_col(position="dodge2", fill="lightgreen", col="grey") + 
        gg_labs + 
        xlim(c(0, 145))
        if (show_error) pmm = pmm + geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2, colour=NA) 
      pfm = ggplot( M[ region %in% REG & sex==1 & mat==1, ], aes(cwd, mn) ) + 
        geom_col(position="dodge2", fill="lightgreen", col="grey") + 
        gg_labs + 
        xlim(c(0, 85))
        if (show_error) pfm = pfm + geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2, colour=NA) 
      pmi = ggplot( M[ region %in% REG & sex==0 & mat==0, ], aes(cwd, mn) ) + 
        geom_col(position="dodge2", fill="lightgreen", col="grey") + 
        gg_labs + 
        xlim(c(0, 145))
        if (show_error) pmi = pmi + geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2, colour=NA) 
      pfi = ggplot( M[ region %in% REG & sex==1 & mat==0, ], aes(cwd, mn) ) + 
        geom_col(position="dodge2", fill="lightgreen", col="grey") + 
        gg_labs + 
        xlim(c(0, 85))
        if (show_error) pfi = pfi + geom_ribbon(aes(ymin=lb, ymax=ub), alpha=0.2, colour=NA) 
  }


  labels = paste( 
    c( 
      "Male mature", 
      "Female mature: ", 
      "Male immature", 
      "Female immature"
    )) 

  plt = ggpubr::ggarrange( pmm, pfm, pmi, pfi, 
    ncol=2, nrow=2, labels=labels, align = "v", font.label=list(size=8) 
  )
  
  plt = ggpubr::annotate_figure(plt, 
    bottom = ggpubr::text_grob( paste( REG, YR, fn_prefix), 
        color = "slategray", hjust = 1, x = 1, face = "plain", size = 8) 
  )

  if (!file.exists(plotsavedir)) dir.create( plotsavedir, showWarnings=FALSE, recursive=TRUE )
  fn_plt = file.path( plotsavedir, paste( fn_prefix, "_", rnm, "_", YR, ".png", sep="") )
  ggsave(filename=fn_plt, plot=plt, dpi=288 )

  message("Plot saved to: ", fn_plt ) 

  return(plt)
}

