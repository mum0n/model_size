
# using predict for kmm

# source: https://discourse.julialang.org/t/using-mixturemodels-from-distributions-jl-in-turing-issues-with-posterior-prediction/71617

using Distributions, Turing

# Unidimensional Kernel Mixture model with K pre-specified components
# that cover the space from min_x to max_x
 

@model function KMM(x, min_x, max_x, k, σ, N=size(x, 1))
    linspan = range(min_x, stop=max_x, length=k)
    kernels = map(u -> Normal(u, σ), linspan)

    ω ~ Dirichlet(k, 1.0)
    mixdist = MixtureModel(kernels, ω)

    x ~ filldist(mixdist, N)
end

# Simulate data from a bimodal distribution
data = vcat(rand(Normal(-1, 0.5), 50), rand(Normal(1, 0.5), 50))

# Define a kernel mixture with 10 gaussian components, with means covering -2:2
model = KMM(data, -2.0, 2.0, 10, 0.5)

# Estimate weights
m1 = sample(model, NUTS(0.65), 1000)

pp_data = predict(KMM(missing, -2.0, 2.0, 10, 0.5, size(data, 1)), m1)


---


using Random, KernelDensity, Distributions
data = randn(1000)
h = KernelDensity.default_bandwidth(data)

newdata = [rand(Normal(rand(data) , h)) for _=1:1000]
std(newdata) # => 1.0014533279921076 (as expected)

---

using Distributions

mm = MixtureModel([Normal(-2.0, 1.2), Normal(), Normal(3.0, 2.5)], [0.1, 0.6, 0.3])
MixtureModel{Normal{Float64}}(K = 3)

components[1] (prior = 0.1000): Normal{Float64}(μ=-2.0, σ=1.2)
components[2] (prior = 0.6000): Normal{Float64}(μ=0.0, σ=1.0)
components[3] (prior = 0.3000): Normal{Float64}(μ=3.0, σ=2.5)

rand(mm)
1.882130062980293
using Plots

histogram(rand(mm, 100_000), normalize = true, xlabel = "Value", ylabel = "Frequency", label = "Mixture model")

