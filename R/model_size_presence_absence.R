

model_size_presence_absence = function(sexid, M, span, modeldir, formula=NULL, action="load", num.threads="4:3") {
     
    fn = file.path( 
        p$modeldir, 
        sexid,
        paste("fit_", sexid, "_", paste0(span, collapse="_"), ".rds", sep="") 
    )
    print(fn)

    fit = NULL
    if (action != "redo" ) {
        if (file.exists(fn)){
            fit = read_write_fast(fn)
            if (!is.null(fit)) return(fit)
        }
    }
 
    # model pa solutions close to final:

    # female:  
# maxld= -209786.4105 fn=4817 theta= -2.5601 -1.0683 1.0352 1.4241 2.7500 0.1777 0.3675 -1.8650 -2.9363 0.5284 2.3864 -2.1478 2.8746 -1.5170 -0.8136 0.9255 [8.24, 3.155]
 

    # male: 
# maxld= .. fn=3781 theta= -1.7892 -1.2806 1.3076 1.7057 1.2670 0.0053 1.8229 -0.6221 1.0522 2.4770 3.4527 -1.9354 4.1381 -1.0566 -0.4377 0.8401 [7.16, 12.138]

    theta0 = switch(
        sexid,
        female = c(-2.5601, -1.0683, 1.0352, 1.4241, 2.7500, 0.1777, 0.3675, -1.8650, -2.9363, 0.5284, 2.3864, -2.1478, 2.8746, -1.5170, -0.8136, 0.9255),
        male   = c(-1.789, -1.280, 1.307, 1.705, 1.267, 0.005, 1.822, -0.622, 1.052, 2.477, 3.452, -1.9354, 4.138, -1.056, -0.437, 0.840)
    )

    fit = inla( 
        data=M, 
        formula=formula, 
        family="binomial", 
        verbose=TRUE, 
        # control.inla = list( strategy="adaptive", int.strategy="eb", h=0.05 ),
        control.predictor = list(compute = TRUE,  link = 1), 
        control.compute=list( dic=TRUE, waic=TRUE, cpo=FALSE, config=TRUE, return.marginals.predictor=TRUE ),
        control.mode= list(theta= theta0),
        num.threads=num.threads
    )

    read_write_fast( fit, file=fn )  # read_write_fast is a wrapper for a number of save/reads ... default being qs::qsave

    return(fit)

}

