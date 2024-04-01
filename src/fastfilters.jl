function nodes_weights(N, kmin, kmax)
    transform = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(transform)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2
    w = (kmax - kmin) / 2 * w
    return x, w
end
