showall(x) = show(stdout, "text/plain", x)


function modes_naive() 
    # these imodes come from naive kernel density analyses of size structure by (in R:):
    # mds = size_distributions(p=p, toget="modal_groups", strata=strata, bw=bw, np=np, ldx=ldx, 
    #    sigdigits=sigdigits, redo=TRUE )
    fn = projectdir( "outputs", "size_structure",  "modes_summary.rdz" )
    mds = load( fn, convert=true)["mds"]
    # limit size range (CW)
    return mds
end
  

function instar_prediction_naive( logcw, sex )
    if sex=="f"
        # cw = exp(2.198848 + 0.315026 * (instar - 4) ) # left bounds
        instar = floor( 4 + ( logcw - 2.198848 ) /  0.315026 )
    end
    
    if sex=="m"
        # cw = exp(1.917564 + 0.298914 * (instar - 3) )  # left bounds
        instar = floor( 3 + ( logcw - 1.917564 ) /  0.298914 )
    end    

    return instar
end



function instar_inference( mds=nothing; savedir= projectdir( "outputs" ), outformat="png"  ) 

    fnout = joinpath(savedir, "mds.jl2" )

    if isnothing(mds)
        if isfile(fnout) 
            mds = jldopen(fnout, "r")["mds"] 
            return mds
        end
    end

    # infer from historical OBSERVATIONS 
    mds[!, :instar] .= instar_prediction_naive.( mds[!, :peakmode], mds[!, :sex] ) 
    sort!(mds, [:sex, :mat, :instar])
    
    fi = findall( (mds.sex.=="f") .& (mds.mat.=="i") )
    fm = findall( (mds.sex.=="f") .& (mds.mat.=="m") )
    mi = findall( (mds.sex.=="m") .& (mds.mat.=="i") )
    mm = findall( (mds.sex.=="m") .& (mds.mat.=="m") )
  
    if false
      # females: dominated by instar 10 mature
      #  infer mature 9 is missing 
      #  infer immature 10 is missing
      pl = plot()
      pl = scatter( pl, mds.instar[fi], mds.peakmode[fi] ) 
      pl = scatter( pl, mds.instar[fm], mds.peakmode[fm]) 
       
  
      # males:  no issues, mature are a little smaller
      pl = plot()
      pl = scatter( pl, mds.instar[mi], mds.peakmode[mi] ) 
      pl = scatter( pl, mds.instar[mm], mds.peakmode[mm]) 
    end
    
  
    # ensure all instars are there
  
    jfi = setdiff( 4:10, unique(mds.instar[fi]) )
    jfm = setdiff( 8:11, unique(mds.instar[fm]) )
    jmi = setdiff( 4:12, unique(mds.instar[mi]) )
    jmm = setdiff( 11:13, unique(mds.instar[mm]) )
  
    o = DataFrame()
    if length(jfi) > 0
      append!(o, DataFrame( sex="f", mat="i", instar=jfi ) )
    end
    
    if length(jfm) > 0
        append!(o, DataFrame( sex="f", mat="m", instar=jfm ) )
    end
    
    if length(jmi) > 0
        append!(o, DataFrame( sex="m", mat="i", instar=jmi ) )
    end
    
    if length(jmm) > 0
        append!(o, DataFrame( sex="m", mat="m", instar=jmm ) )
    end

    if size(o, 1) > 0
        mds = outerjoin( mds, o, on=names(o) )
    end

    i = findall(nonunique(mds[:,[:sex, :mat, :instar]]))
    if length(i) > 0
      delete!(mds, i) 
    end
  
    sort!(mds, [:sex, :mat, :instar])
 
    fi = findall( (mds.sex.=="f") .& (mds.mat.=="i") )
    fm = findall( (mds.sex.=="f") .& (mds.mat.=="m") )
    mi = findall( (mds.sex.=="m") .& (mds.mat.=="i") )
    mm = findall( (mds.sex.=="m") .& (mds.mat.=="m") )
  
    modfi = lm(@formula(peakmode ~ instar), mds[fi,:])
    modfm = lm(@formula(peakmode ~ instar), mds[fm,:])
    modmi = lm(@formula(peakmode ~ instar), mds[mi,:])
    modmm = lm(@formula(peakmode ~ instar), mds[mm,:])
      
    #────────────────────────────────────────────────────────────────────────
    #                Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
    #────────────────────────────────────────────────────────────────────────
    #(Intercept)  1.05187   0.0430406   24.44    <1e-04   0.932365   1.17136
    #instar       0.321978  0.00640426  50.28    <1e-06   0.304196   0.339759
  
    #(Intercept)  1.49157    0.216241    6.90    0.0917  -1.25602     4.23917
    #instar       0.268652   0.0221858  12.11    0.0525  -0.0132459   0.550549
  
    #(Intercept)  1.10376   0.0259636   42.51    <1e-08   1.04237    1.16515
    #instar       0.308079  0.00308857  99.75    <1e-11   0.300775   0.315382
  
    #(Intercept)  0.335281  0.114697     2.92    0.2098  -1.12208    1.79264
    #instar       0.363371  0.00953601  38.11    0.0167   0.242204   0.484537
    #────────────────────────────────────────────────────────────────────────
      
    mds.predicted .= 0.0 
    mds.predicted[fi] .= predict(modfi, mds[fi,:])
    mds.predicted[fm] .= predict(modfm, mds[fm,:])
    mds.predicted[mi] .= predict(modmi, mds[mi,:])
    mds.predicted[mm] .= predict(modmm, mds[mm,:])
    sort!(mds, [:sex, :mat, :instar])
 
    inst = string.( Int32.(mds.instar) ) 
    ilen = findall( length.( inst) .==1 )
    inst[ ilen] = string.("0", inst[ ilen] )

    mds.stage .= string.( mds.sex, "|", mds.mat, "|", inst ) 
     
    # plots
    labelsize = 12
    
    # GR.setarrowsize(1.25)
    
    # using CairoMakie 
    
    # Makie does not work with missing values
    o = mds[:, [:sex, :mat, :instar, :peakmode]]
    dropmissing!(o)
    ofi = findall( (o.sex.=="f") .& (o.mat.=="i") )
    ofm = findall( (o.sex.=="f") .& (o.mat.=="m") )
    omi = findall( (o.sex.=="m") .& (o.mat.=="i") )
    omm = findall( (o.sex.=="m") .& (o.mat.=="m") )
 
    CairoMakie.activate!(type = "png")
    
    f = Figure()
    ax = Axis( 
        f[1, 1],
        title = "Growth in females",
        xlabel = "Stage (instar)",
        ylabel = "Carapace width (mm)",
        xticks = Int32.(sort( unique( mds.instar[ vcat( fi, fm)] ) ) )
    )

    CairoMakie.scatter!( ax, o.instar[ofi], exp.(o.peakmode[ofi]); 
        label = "Observed immature",
        marker=:diamond, alpha=0.5, color=:orange, markersize=12
    )
     
    CairoMakie.scatter!( ax, o.instar[ofm], exp.(o.peakmode[ofm]); 
        label = "Observed mature",
        marker=:circle, alpha=0.5, color=:red, markersize=12
    )

    fi = findall( (mds.sex.=="f") .& (mds.mat.=="i") )
    fm = findall( (mds.sex.=="f") .& (mds.mat.=="m") )
     
    CairoMakie.scatter!( ax, mds.instar[fi], exp.(mds.predicted[fi]); 
        label = "Predicted immature",
        marker=:xcross, alpha=0.5, color=:blue, markersize=12
    )
     
    CairoMakie.scatter!( ax, mds.instar[fm], exp.(mds.predicted[fm]); 
        label = "Predicted mature",
        marker=:cross, alpha=0.5, color=:darkblue, markersize=12
    )
 
    # immature
    for i in 1:(length(fi)-1)
        j = fi[[i, i+1]] 
        arrowlines!(ax, mds.instar[j], exp.( mds.predicted[j]); arrowstyle="--|>", alpha=0.6, color=:darkorange, linewidth=1.5)
    end

    # mature
    for i in 1:(length(fm))
        k = only( findall( (Int32.(mds.instar).==(mds.instar[fm[i]]-1.0) ) .& (mds.mat.=="i") .& (mds.sex.=="f")) )
        j = [k, fm[i]]
        arrowlines!(ax, mds.instar[j], exp.( mds.predicted[j]); arrowstyle="..|>", alpha=0.6, color=:darkblue, linewidth=1.5  )
    end
   
    axislegend(position = :lt)
  
    save(joinpath( savedir, string("plot_growth_female.", outformat ) ) , f ) 
   

    f = Figure()
    ax = Axis( 
        f[1, 1],
        title = "Growth in males",
        xlabel = "Stage (instar)",
        ylabel = "Carapace width (mm)",
        xticks = Int32.(sort( unique( mds.instar[ vcat( mi, mm)] ) ) )
    )

    CairoMakie.scatter!( ax, o.instar[omi], exp.(o.peakmode[omi]); 
        label = "Observed immature",
        marker=:diamond, alpha=0.5, color=:orange, markersize=12
    )
     
    CairoMakie.scatter!( ax, o.instar[omm], exp.(o.peakmode[omm]); 
        label = "Observed mature",
        marker=:circle, alpha=0.5, color=:red, markersize=12
    )

    mi = findall( (mds.sex.=="m") .& (mds.mat.=="i") )
    mm = findall( (mds.sex.=="m") .& (mds.mat.=="m") )
    
    CairoMakie.scatter!( ax, mds.instar[mi], exp.(mds.predicted[mi]); 
        label = "Predicted immature",
        marker=:xcross, alpha=0.5, color=:blue, markersize=12
    )
     
    CairoMakie.scatter!( ax, mds.instar[mm], exp.(mds.predicted[mm]); 
        label = "Predicted mature",
        marker=:cross, alpha=0.5, color=:darkblue, markersize=12
    )
 
    # immature
    for i in 1:(length(mi)-1)
        j = mi[[i, i+1]] 
        arrowlines!(ax, mds.instar[j], exp.( mds.predicted[j]); arrowstyle="--|>", alpha=0.6, color=:darkorange, linewidth=1.5)
    end

    # mature
    for i in 1:(length(mm))
        instar0 = Int32(mds.instar[mm[i]]-1.0)
        k = only( findall( (Int32.(mds.instar).==instar0) .& (mds.mat.=="i") .& (mds.sex.=="m")) )
        j = [k, mm[i]]
        arrowlines!(ax, mds.instar[j], exp.( mds.predicted[j]); arrowstyle="..|>", alpha=0.6, color=:darkblue, linewidth=1.5  )
    end
   
    axislegend(position = :lt)
  
    save( joinpath( savedir, string("plot_growth_male.", outformat ) ), f ) 
   
    jldsave(fnout; mds)

    return mds
end
  


function size_structured_data(; sigdigits=3 ) 
    # data created in R using : Y = size_distributions( p=p, toget="base_data" ) 
    fn = projectdir( "outputs", "size_structure",  "base_data_julia.rdz" )
    o = load( fn, convert=true)["out"]
    
    Y  = o["Y"] 
    Y.logcw = log.(Y.cw)
    Y.logsa = log.(Y.sa)
    Y.wt = 1.0 ./ Y.sa
    # data_resolution is ~1 to 2 mm (observations error)
    # bw =  round(mean(log.(22:142)-log.(20:140)), digits=2 )   # approx SD is interval
  
    # discretize time (weekly)
    Y[!,:ti] = Y[!,:year] .+ round.(trunc.(Y[!,:julian] ./ 365. .* 52 ) ./ 52, digits=sigdigits) # weekly 
  
    # dimensionality of problem
    nb = o["nb"]
    aus = o["au"]
    
    # recode sex and mat
     
    Y = Y[ in(["0", "1"]).(Y.sex) .& (.!ismissing.(Y.sex)),:]
    Y = Y[ in(["0", "1"]).(Y.mat) .& (.!ismissing.(Y.mat)),:]
    Y = Y[ .!ismissing.(Y.cw),:]
    Y = Y[ .!ismissing.(Y.wt),:]

    
    # sex codes
        # male = 0
        # female = 1
        # sex.unknown = 2
  
    # # maturity codes
        # immature = 0
        # mature = 1
        # mat.unknown = 2
    
    sexmap = Dict("0" => "m", "1" => "f")
    matmap = Dict("0" => "i", "1" => "m")
  
    Y.sex = [sexmap[i] for i in Y[:, "sex"]]
    Y.mat = [matmap[i] for i in Y[:, "mat"]]
    
    Y.cw = identity.( Y.cw )
    Y.logcw = identity.( Y.logcw )
    Y.wt = identity.( Y.wt )

    return Y, nb, aus
end


function kernel_density_estimated_modes(Y=nothing; 
        nb=nb, 
        aus=aus, 
        yrs=nothing, 
        xrange=[2, 250],  # wider to allow FFT-based Kernel density estimation
        ti_window = (-4,4), multiplier=20.0, nmin=3, 
        toget="summary" 
    ) 
 
    fnout_root = "kernel_density_modes"
    savedir = projectdir( "outputs", "size_structure", "kernel_density_modes" )
    mkpath(savedir) 
    
    if isnothing(yrs)
        if isnothing(Y) 
            yrs = Int.(sort( unique( Y.year ) ))
        else
            return( "yrs need to bespecified at a minimum" )
        end
    end


    res = DataFrame()
   
    if toget=="saved_results"
        for yr in yrs
            fnout = joinpath(savedir, string(fnout_root, "_", yr, ".jl2") )
            if isfile(fnout) 
                res = vcat( res, jldopen(fnout, "r")["res"] )
            end
        end
        return res 
    end
   
    np = 2048  # default number of points in kde
    xr = round.( log.(xrange), digits=2 ) 
    ldx = (diff(xr)/(np-1))[1]
    xvals = range( xr[1], xr[2], step = ldx )
  
    bw = ldx*multiplier * 0.9   # empirically reasonable (~0.05)
    pw = Int(floor( bw / ldx))   # length of peak window in one direction (index, not x-value)
   
    if toget=="summary"
        # summarize mode of modes
        kdem = kernel_density_estimated_modes(Y; yrs=yrs, nmin=nmin, toget="saved_results") 
        res = DataFrame()
        for sx in sexes 
            ks = findall( x -> x==sx, kdem.sex )  
            (length(ks) < nmin) && continue
            for ma in mats
                kmi = findall( x -> x==ma, kdem.mat[ks])  
                N =  length(kmi) 
                (N < nmin) && continue 
                km = ks[kmi]  # back on kdem
                U = kde( kdem.peakmode[km]; weights=kdem.Neffective[km], boundary=Tuple(xr), bandwidth=bw, kernel=Normal )  # normal kernal by default
                pks = findmaxima( U.density, pw )
                pks = peakproms( pks; min=0.001)  # compute prominences
                pks = peakwidths( pks; min=0.01)  # compute widths at 50% prominences
                # DataFrame(pks[Not(:data)])
                pl = plotpeaks( U.x, U.density; peaks=pks.indices, col ="green")
                out1 = repeat( DataFrame( sex=sx, mat=ma, Neffective=round( sum( kdem.Neffective[km]) ) ), length(pks.indices))
                out2 = DataFrame( peakmode=U.x[pks.indices], peakheight=U.density[pks.indices], peakwidth=pks.widths )
                out = hcat(out1, out2)
                append!(res, out)
               end # mat
        end # sex
    
        return res
    end
  
  
    for yr in yrs 
    
      res = DataFrame()
    
      for wk in 1:52 
        mti1 = yr .+ (wk .+ ti_window[1]) ./ 52
        mti2 = yr .+ (wk .+ ti_window[2]) ./ 52
        kt = findall( x -> ( mti1 <= x <= mti2), Y.ti)
  
        (length(kt) < 1)  && continue 
        auids = unique(Y.space_id[kt])
  
        for auid in auids 
          
          ii = findall(x -> x == auid, aus)
          (length(ii) != 1) && continue  # NA's  should not happen
          oo = reduce(vcat, nb[ii])
          nbs = aus[ oo ]  # nearest neighbours
          auloc = unique( vcat( nbs, auid ) )
          kain = in(auloc)
          kai = findall( x -> kain.(x), Y.space_id[kt])  
          (length(kai) < 1) && continue
          ka = kt[kai] # index back on original Y elements 
  
          for sx in sexes 
              
            ksi = findall( x -> x==sx, Y.sex[ka])  
            (length(ksi) < nmin) && continue
            ks = ka[ksi] # back on Y
             
            for ma in mats
              kmi = findall( x -> x==ma, Y.mat[ks])  
              N =  length(kmi) 
              (N < nmin) && continue 
              km = ks[kmi]  # back on Y
       
              U = kde( Y.logcw[km]; weights=Y.wt[km], boundary=Tuple(xr), bandwidth=bw, kernel=Normal )  # normal kernal by default
              
              pks = findmaxima( U.density, pw )
              pks = peakproms( pks; min=0.001)  # compute prominences
              pks = peakwidths( pks; min=0.1)  # compute widths at 50% prominences
  
              # DataFrame(pks[Not(:data)])
              # pl = plotpeaks( U.x, U.density; peaks=pks.indices, col ="green")
              # pl = plot(pl, U  )  # KD smooth
  
              out1 = repeat( DataFrame( sex=sx, mat=ma, au=auid, year=yr, wk=wk, Nsample=N, Neffective=round( sum( Y.wt[km]) ) ), length(pks.indices))
   
              out2 = DataFrame( peakmode=U.x[pks.indices], peakheight=U.density[pks.indices], peakwidth=pks.widths )
   
              out = hcat(out1, out2)
              # @show out
   
              append!(res, out)
   
            end # mat
          end # sex
  
        end  # au
      end # wk
  
      fnout = joinpath(savedir, string(fnout_root, "_", yr, ".jl2") )
      jldsave(fnout; res)
      print( string("\n year complete: ", yr, "\n") ) 
    
    end # year
    
    return kernel_density_estimated_modes(Y; yrs=yrs, toget="saved_results" ) 
  
end
   

  
Turing.@model function kmm(x, imodes, n_imodes, N, sd_imodes )
  # Kernel Mixture model with imodes as pre-specified components
  # alpha = concentration parameter of results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
  sigmasq ~ filldist(  truncated( InverseGamma(5.0, 0.1), 0.0, sd_imodes), n_imodes)  # variance prior 
  kernels = map( i -> Normal( imodes[i], sqrt(sigmasq[i]) ), 1:n_imodes ) 
  alpha ~ Dirichlet(n_imodes, 1.0)
  mixdist = Distributions.MixtureModel(kernels, alpha)
  x ~ filldist(mixdist, N)
end


Turing.@model function kmm_latent(x, imodes0, n_imodes, N, sd_imodes )
    # Kernel Mixture model with imodes as pre-specified components
    # alpha = concentration parameter of results in all sets of probabilities being equally likely, i.e., in this case the Dirichlet distribution of dimension k is equivalent to a uniform distribution over a k-1-dimensional simplex.  
    sigmasq ~ filldist( LogNormal( log(sd_imodes), 0.25), n_imodes)  # variance prior 
    imodes ~ arraydist( Normal.(imodes0, sd_imodes) )  # variance prior
    kernels = map( i -> Normal( imodes[i], sqrt(sigmasq[i]) ), 1:n_imodes ) 
    alpha ~ Dirichlet(n_imodes, 1.0)
    alpha_sum = sum(alpha)
    alpha_sum ~ Normal(1.0, 0.001*n_imodes )  # soft sum to 1 constraint (to make it a probability)
    mixdist = Distributions.MixtureModel(kernels, alpha)
    x ~ filldist(mixdist, N)
end

 
function kmm_samples( chain, N; imodes=nothing, modeltype="latent" ) 

    if modeltype == "latent"
        n_imodes = size( group(chain, "imodes"), 2 ) 
    elseif modeltype == "deterministic"
        n_imodes =  length(imodes)  # imodes must be passed
    end

    nch = size(chain)
    out = Any[];
    for a in 1:N
        g = rand(1:nch[1])  # sim
        l = rand(1:nch[3])  # chain
        sigmasq = [ chain[g, Symbol("sigmasq[$k]"), l] for k in 1:n_imodes]
        alpha = [ chain[g, Symbol("alpha[$k]"), l] for k in 1:n_imodes]
        if modeltype == "latent"
            imodes = [ chain[g, Symbol("imodes[$k]"), l] for k in 1:n_imodes]
        end
        kernels = map( i -> Normal(imodes[i], sqrt(sigmasq[i]) ), 1:n_imodes ) 
        sep = Any[]
        for b in 1:n_imodes
            if rand() <= alpha[b]
                alph = zeros(n_imodes) 
                alph[b] = 1
                mixdist = Distributions.MixtureModel(kernels, alph)
                push!(sep, rand(mixdist ) )
            else 
                push!(sep, missing)
            end
        end
        push!(out, sep)
    end  
    return reduce(hcat, out)'  
end


function kmm_sample( res, N=1, n_imodes=1 ) 
    nch = size(res)
    out = Any[];
    for a in 1:N
        g = rand(1:nch[1]) 
        sigmasq = [ res[g, Symbol("sigmasq[$k]") ] for k in 1:n_imodes]
        alpha = [ res[g, Symbol("alpha[$k]") ] for k in 1:n_imodes]
        imodes = [ res[g, Symbol("imodes[$k]") ] for k in 1:n_imodes]
        kernels = map( i -> Normal(imodes[i], sqrt(sigmasq[i]) ), 1:n_imodes ) 
        sep = Any[]
        for b in 1:n_imodes
            if rand() <= alpha[b]
                alph = zeros(n_imodes) 
                alph[b] = 1
                mixdist = Distributions.MixtureModel(kernels, alph)
                push!(sep, rand(mixdist ) )
            else 
                push!(sep, missing)
            end
        end
        push!(out, sep)
    end
    return reduce(hcat, out)'  
end;


function size_structured_data_subset(U, mds; 
    yrs=nothing, sexes=nothing, mats=nothing, regions=nothing, sids=nothing )

    if !isnothing(yrs) 
        U = U[ findall( in( Int64.(yrs)) , Int64.(U.year) ), :]
    end

    if !isnothing(sexes) 
        U = U[ findall( occursin(sexes), U.sex), :]
    end
    
    if !isnothing(mats) 
        U = U[ findall( occursin(mats), U.mat), :]
    end
    
    if !isnothing(regions) 
        if regions != "cfaall"
            U = U[ findall( occursin(regions), U.region), :]
        end
    end

    if !isnothing(sids) 
        U = U[ findall( occursin(sids), U.sid), :]
    end

    data = U[ :, :logcw ]  
    N = length(data)
 
    if !isnothing(sexes) 
        mds = mds[ findall( occursin(sexes), mds.sex), :]
    else
        print("no sexes specified, averaging .. make sure this is what you want to do")
    end

    if !isnothing(mats) 
        mds = mds[ findall( occursin(mats), mds.mat), :]
    else 
        print("no mats specified, averaging .. make sure this is what you want to do")
    end
    
    imodes = combine(groupby(mds, :instar), :predicted => mean) 
    
    n_imodes = size(imodes, 1)  
    # assume ~ 3SD so 1sd = xrez/3
    sd_imodes = round( maximum(diff(imodes[:,:predicted_mean])), digits=2 ) / 6.0 # approx magnitude of variability of size classes see :predicted_se
    
    return( data, N, imodes, n_imodes, sd_imodes )
end


function kmm_chain( U, mds; modeltype="latent", yrs=1999:2000, sexes=["m", "f"], mats=["i", "m"], regions=nothing,
    nmin=30, nsamples=1000, n_chains=4, sampler = Turing.NUTS(0.65), savedir=tempdir() ) 

    # create mcmc chain using deterministic imodes, for each year, sex, mat, region
    if modeltype == "latent"
        kmm_model = kmm_latent
    elseif modeltype == "deterministic"
        kmm_model = kmm
    end

    for yr in yrs
        print("year:", yr, "\n")
        UU = U[ findall( in( Int64.(yr)) , Int64.(U.year) ), :]
        
        savedir_yr = joinpath( savedir, string(yr)  )
        mkpath(savedir_yr) 
        
        size(UU,1) < nmin && continue

        if isnothing(regions)
            spaceids = unique( UU.sid ) 
        else
            spaceids = unique( UU.region ) 
        end

        for spaceid in spaceids
        for sex in sexes
        for mat in mats
            fnout = joinpath( savedir_yr, string( "kmm_", spaceid, "_", sex, "_", mat, "_", yr, ".jl2" ))
            isfile(fnout) && continue

            if !isnothing(regions)
                data, N, modes_mds, n_imodes, sd_imodes = size_structured_data_subset( UU, mds; yrs=yr, sexes=sex, mats=mat, regions=spaceid)
            else
                data, N, modes_mds, n_imodes, sd_imodes = size_structured_data_subset( UU, mds; sexes=sex, mats=mat, sids=spaceid)
            end
            
            N < nmin && continue 
            print(" ", N )

            imodes = modes_mds[:, :predicted_mean]
            M = kmm_model(data, imodes, n_imodes, N, sd_imodes )
            # chain = sample(M, sampler, nsamples; max_depth=8, init_ϵ=0.0125)  
            chain = sample(M, sampler, MCMCThreads(), nsamples, n_chains; max_depth=8, init_ϵ=0.01)  
            jldsave(fnout; chain)
            print(fnout)
        end
        end 
        end
        
    end
    
    print( savedir )
    return "Models completed and saved"
end

  

function kmm_chain_load( ; yr=nothing, sex=nothing, mat=nothing, region=nothing, sid=nothing, fnout=nothing, savedir=tempdir() ) 
    # get a specific chain from file
    if isnothing(fnout) 
        savedir_yr = joinpath( savedir, string(yr) )
        if !isnothing(region)
            fnout = joinpath(savedir_yr, string( "kmm_", region, "_", sex, "_", mat, "_", yr, ".jl2" ))
        end
        if !isnothing(sid)
            fnout = joinpath(savedir_yr, string( "kmm_", sid, "_", sex, "_", mat, "_", yr, ".jl2" ))
        end
    end
    # print(fnout)
    chain = missing
    if isfile(fnout) 
        chain = jldopen(fnout, "r")["chain"] 
    end
    return(chain)
end
 
 
function file_list( searchdir=".", patt=r".*" )
    fnall = readdir(searchdir)
    fns = [fn for fn in fnall if isfile(abspath(joinpath(searchdir, fn)))]
    fns = [m.match for m in match.(patt, fns) if m != nothing]
    return fns
end

function parse_filename(fn)
    fn = replace(fn, "kmm_" => "")
    fn = replace(fn, ".jl2" => "")
    fnx = split(fn, "_")
    return fnx
end


function kmm_summary( ; mds=Any[], yrs=1999:2000, sexes=["m", "f"], mats=["i", "m"], 
    regions="cfaall", toget="saved_results", savedir=tempdir(), modeltype="latent", save_r_data=false )
 
    out = DataFrame()
    fnout_root = "kmm_parameter_summaries"

    if toget=="saved_results"
        for yr in yrs
            fnout = joinpath(savedir, string(fnout_root, "_", yr, ".jl2") )
            if isfile(fnout) 
                append!( out, jldopen(fnout, "r")["out"] )
            end
        end
        return out
    end

    for yr in yrs
        fns = file_list( joinpath(savedir, string(yr) ),  r"^kmm.*\.jl2$" )

        length(fns) == 0 && continue

        fnout = joinpath(savedir, string(fnout_root, "_", yr, ".jl2") )
        out = DataFrame()
        for fn in fns
            chain = kmm_chain_load( fnout=joinpath(savedir, string(yr), fn) ) 
            ismissing(chain) && continue
            
            if modeltype == "latent"
                sid, sex, mat, yr = parse_filename(fn)
                n_imodes = size( group(chain, "imodes"), 2 ) 
                imodes =  [ reduce(hcat, chain[:, Symbol("imodes[$k]"), :])' for k in 1:n_imodes ] 
                sigmasq = [ reduce(hcat, chain[:, Symbol("sigmasq[$k]"), :])' for k in 1:n_imodes ] 
                alpha =   [ reduce(hcat, chain[:, Symbol("alpha[$k]"), :])'   for k in 1:n_imodes ] 
                o = DataFrame( 
                    sid=fill(sid, n_imodes), 
                    yr=fill(yr, n_imodes), 
                    sex=fill(sex, n_imodes), 
                    mat=fill(mat, n_imodes),
                    imodes = [ mean(imodes[k]) for k in 1:n_imodes ],
                    instar=mds[(mds.sex.==sex) .& (mds.mat.==mat), :instar],
                    stage=mds[(mds.sex.==sex) .& (mds.mat.==mat), :stage],
                    sigmasq_mean = [ mean(sigmasq[k]) for k in 1:n_imodes ],
                    sigmasq_sd = [ std(sigmasq[k]) for k in 1:n_imodes ], # ~ 30-40% of mean
                    alpha_mean = [ mean(alpha[k]) for k in 1:n_imodes ],
                    alpha_sd = [ std(alpha[k]) for k in 1:n_imodes ]  
                )
                #  fill(sid, nins)
            end

            if modeltype == "deterministic"
                region, sex, mat, yr = parse_filename(fn)
                imodes = mds[ (mds.sex.==sex) .& (mds.mat.==mat), :predicted ]
                n_imodes = length(imodes)
                sigmasq = [ reduce(hcat, chain[:, Symbol("sigmasq[$k]"), :])' for k in 1:n_imodes ] 
                alpha =   [ reduce(hcat, chain[:, Symbol("alpha[$k]"), :])'   for k in 1:n_imodes ] 
                o = DataFrame( region=region, yr=yr, sex=sex, mat=mat,
                    imodes =mds[(mds.sex.==sex) .& (mds.mat.==mat), :predicted],
                    instar=mds[(mds.sex.==sex) .& (mds.mat.==mat), :instar],
                    stage=mds[(mds.sex.==sex) .& (mds.mat.==mat), :stage],
                    sigmasq_mean = [ mean(sigmasq[k]) for k in 1:n_imodes ],
                    sigmasq_sd = [ std(sigmasq[k]) for k in 1:n_imodes ], # ~ 30-40% of mean
                    alpha_mean = [ mean(alpha[k]) for k in 1:n_imodes ],
                    alpha_sd = [ std(alpha[k]) for k in 1:n_imodes ]  
                )
            end
            append!(out, o) 
        end

        jldsave(fnout; out)
        println( fnout )
        if save_r_data
            # RCall.@rput out
            fnR = replace(fnout, "jl2"=>"rdz")
            cmdR = string( "read_write_fast( out, fn=", "'$fnR')" )
            # RCall.reval(cmdR)
        end
    end
    return kmm_summary( ; toget="saved_results", yrs=yrs, savedir=savedir )
end
    



function normalize_size_structure(U; sexes="0", mats="0", aus="7", yrs=2009, weeks=32, ti_window=[-1,1],
  bw=0.1,  nsmp=1024, np=256, xr=extrema(U), nmin=5, 
  outdir_nss=joinpath( project_directory, "size_structure", "posteriors_summaries" ), overwrite=false )   

  # deprcated: now done in R (bio.snowcrab::size_distributions.R)

  if false
      # debug: 
      sex="0"; mat="0"; au="103"; yr=1999; week=14
      sex="1"; mat="1"; au="304"; yr=1999; week=32
      sex="0"; mat="1"; au="318"; yr=1999; week=32
      
      nsmp=1024
  end


  for yr in yrs
      fn = string( "posterior_summaries_", yr, ".csv" )
      fnout  = joinpath( outdir_nss, fn )
      sample_the_data = true
      if isfile(fnout) 
          if !overwrite
              sample_the_data = false
          end
      end

      if sample_the_data

          out1 = DataFrame( sex="-1", mat="-1", au="-1", year=-1, week=-1, Nsample=-1, Neffective=-1 )
          out2 = DataFrame( Tables.table(zeros(np)') )

          for week in weeks
              mti = yr .+ round.( (week .+ ti_window) ./ 52, digits=3) # weekly   week + ti_window
              kt = findall( x -> x >= mti[1] && x <= mti[2], skipmissing(U.ti) )
            if length(kt) >= nmin
          for au in aus
              ka = intersect( kt, findall( x -> x==au , skipmissing(U.space_id) ) )
            if length(ka) >= nmin
          for sex in sexes
              ks = intersect( ka, findall( x -> x==sex, skipmissing(U.sex) ) )
            if length(ks) >= nmin
          for mat in mats
              n =  intersect( ks, findall( x -> x==mat, skipmissing(U.mat) ) )
              N =  length(n)
              tout = string("| sex: ", sex, "| mat: ", mat, "| au: ", au, "|year: ", yr, "| week: ", week, "| N: ", N ) 
              if N >= nmin 
                  @show tout
                  uu = kde( U.logcw[n], bandwidth=bw, boundary=xr, npoints=np, weights=U.wt[n] )
                  # CairoMakie.lines!( uu , color="orange" )
                  push!( out1, [sex , mat, au, yr, week, N, round( sum( U.wt[n]) ) ] )
                  push!( out2, uu.density' )
                  res = nothing
              end
          end # mat
            end # sex if
          end # sex
            end # au if
          end  # au
            end # time if
          end   # time  
                      
          CSV.write(fnout, hcat(out1, out2), writeheader=true)
          @show fnout 


      end
  end
end
   

function sample_posteriors(U; sexes="0", mats="0", aus="7", yrs=2009, weeks=32, 
  bw=0.1, xrez=bw*4.0, nsmp=1024, np=256, xr=extrema(U),
  ti_window=[-1,1], sampler=Turing.NUTS(0.65), nsamples=1000, n_chains=1, nmin=100, 
  savedir=joinpath( project_directory, "size_structure", "posteriors_summaries" ), overwrite=false )   

  if false
      # debug: 
      sex="0"; mat="0"; au="103"; yr=1999; week=14
      sex="1"; mat="1"; au="304"; yr=1999; week=32
      
      nsmp=1024
      xrez= round( maximum(diff(imodes)), digits=2 )
  end

  out = nothing

  for yr in yrs
      fn = string( "posterior_summaries_", yr, ".csv" )
      fnout  = joinpath( savedir, fn )
      sample_the_data = true
      if isfile(fnout) 
          if !overwrite
              sample_the_data = false
          end
      end

      if sample_the_data

          for week in weeks
              mti = yr .+ round.( (week .+ ti_window) ./ 52, digits=3) # weekly   week + ti_window
              kt = findall( x -> x >= mti[1] && x <= mti[2], skipmissing(U.ti) )
            if length(kt) >= nmin
          for au in aus
              ka = intersect( kt, findall( x -> x==au , skipmissing(U.space_id) ) )
            if length(ka) >= nmin
          for sex in sexes
              ks = intersect( ka, findall( x -> x==sex, skipmissing(U.sex) ) )
            if length(ks) >= nmin
          for mat in mats
              n =  intersect( ks, findall( x -> x==mat, skipmissing(U.mat) ) )
              N =  length(n)
              tout = string("| sex: ", sex, "| mat: ", mat, "| au: ", au, "|year: ", yr, "| week: ", week, "| N: ", N ) 
              if N >= nmin 
                  @show tout
                  uu = kde( U.logcw[n], bandwidth=bw, boundary=xr, npoints=np, weights=U.wt[n] ) # base2 .. large enough to cover x .. fft
                  # CairoMakie.lines!( uu , color="orange" )
                       
                  Neff = round( sum( U.wt[n]) )
                  vv = sample( uu.x, ProbabilityWeights(uu.density), nsmp )
                  M = kmm( vv, imodes, n_imodes, nsmp )
                  
                  chain = nothing
                  if n_chains==1
                      chain = sample(M, sampler, nsamples)
                  else
                      chain = sample(M, sampler, MCMCThreads(), nsamples, n_chains)
                  end

                      if false
                          showall(chain)
                          v = imodes
                          v = fill(4.14, n_imodes)
                          mp = predict( kmm(missing, v, n_imodes, nsmp  ), chain)
                          k = mp.value.data[1:1000,:,:]

                          # using CairoMakie

                          f = Figure()
                          ax = Axis(f[1, 1])

                          # xlim=(imodes[[1,n_imodes]]) 
                          for i in 1:N 
                              den = kde(vec(k[:,i,:]) )
                              CairoMakie.lines!( den, color="red" )
                          end
                          CairoMakie.lines!( kde(U.logcw[n]), color="blue" )
                          f

                      end

                  if !isnothing(chain) 
                      res = DataFrame( summarystats(chain) )
                      res[!,:sex]  .= sex
                      res[!,:mat]  .= mat
                      res[!,:au]   .= au
                      res[!,:year] .= yr
                      res[!,:week] .= week
                      res[!,:N]    .= Neff
                      out = vcat(out, res)
                      res = nothing
                  end

              end
          end # mat
            end # sex if
          end # sex
            end # au if
          end  # au
            end # time if
          end   # time  
          
          CSV.write(fnout, out, writeheader=true)
          @show fnout 
          out = nothing

      end
  end
end
 


  

