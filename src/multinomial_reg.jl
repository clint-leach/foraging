using Turing
using Plots
using StatsPlots
using PDMats
using LinearAlgebra

# Number of bouts
nb = 100

# Number of dives per bout
nd = 5

# Number of prey categories
np = 2

# Generate covariates at bout level
X = hcat(fill(1.0, nb), rand(Normal(0, 1), nb, 1))

# Response to driver
# β = rand(Normal(0, 0.5), 2, np)
β = [-1.0 -1.0; 0.0 -0.0]

# Generating fake data
Σp = ScalMat(2, 2.0)
ϵ = rand(MvNormal(Σp), nb)

μ = X * β + ϵ'
λ = sum(exp.(μ), dims = 2)
histogram(λ)

ϕ = exp.(μ) ./ λ
histogram(ϕ[:, 1], bins = 20)

N = [rand(Poisson(λ[i])) for i = 1:nb, j = 1:nd]
histogram(vec(N))

y = [rand(Multinomial(N[i, j], ϕ[i, :])) for i = 1:nb, j = 1:nd]

# Defining the model
@model function multinom_regression(X, N, y)

    # Dimensions
    p = size(X, 2)
    nb = size(N, 1)
    nd = size(N, 2)
    np = length(y[1, 1])

    # Priors
    β = Matrix{Float64}(undef, p, np)
    for k in 1:np
        β[:, k] ~ MvNormal(ScalMat(p, 1.0))
    end

    # Linear predictors
    μ = X * β
    λ = sum(exp.(μ), dims = 2)
    π = exp.(μ) ./ λ

    # Likelihood
    for i in 1:nb
        for j in 1:nd
            N[i, j] ~ Poisson(λ[i])
            y[i, j] ~ Multinomial(N[i, j], π[i, :])
        end
    end
end

iterations = 10000
ϵ = 0.05
τ = 10

chain = sample(multinom_regression(X, N, y), HMC(ϵ, τ), iterations)

plot(chain)
chain[:β]
