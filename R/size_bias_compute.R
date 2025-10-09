size_bias_compute = function(p) {

  size_selectivity = model_size_results( p=p, todo="size_selectivity"  ) # alternative
  fixed_effects = model_size_results( p=p, todo="fixed_effects"  )

  # reformat
  names_ss = c( "log_cwd", "mean", "sd", "lb", "median", "ub", "mode", "kld", "bioclass")
  names_fe = c( "b0_mean","b0_sd", "b0_lb", "b0_median", "b0_ub", "b0_mode", "b0_kld", "bioclass")

  out = NULL

  for (bc in c("f.imm", "f.mat", "m.imm", "m.mat") ) {
    ss = data.table( size_selectivity[[ bc ]], bioclass = bc )
    names(ss) = names_ss
  
    fe = data.table( fixed_effects[[bc]], bioclass=bc )
    names(fe) = names_fe

    ss = fe[ ss, on=.(bioclass)]
    out = rbind(out, ss)
  }

  out[, bias:=mean + b0_mean ]

  return(out)
}

