function T(x)
    return (3 / x^3) * (sin(x) - x * cos(x))#exp(-x^2)#
end

function Δ²(k, pk)
    return 4π*(k/(2π))^3*pk(k)
end

function dσ2dk(k, R, pk)
    x = k * R
    T² = T(x)^2
    dσ2dk = T² * Δ²(k, pk) / k
    return dσ2dk
end

function sigma2quadgk(pk, R::Real; kmin = 1e-5, kmax = 20/R)
    σ², _ = quadgk(x -> dσ2dk(x, R, pk), kmin, kmax)
    return σ²
end
