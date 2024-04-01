function R_NL_z_computation_turbo!(result, w, Wx, Δ²kz)
    @turbo for i ∈ axes(Wx,1), j ∈ axes(Δ²kz,1)
        Cmn = zero(eltype(w))
        for k ∈ axes(w,1)
            Cmn += w[k]*Wx[i,k]*Δ²kz[j,k]
        end
        result[i,j] = Cmn
    end
end

function get_R_NL_z!(R_NL_z, z_list, R_test_interp, pippami)
    for (zidx, myz) in enumerate(z_list)

        pippo = AkimaInterpolation(log10.(pippami[:,zidx]), log10.(R_test_interp), extrapolate=true)
        myf(x) = @. 10^pippo(log10(x))-1.686
        nlsol = nlsolve(myf, [ 0.0002])
        R_NL_z[zidx] = nlsol.zero[1]
    end
    return nothing
end

function get_R_NL_z(z_list, R_test_interp, w1, Δ²kz, Wx)
    R_NL_z = zeros(length(z_list))
    pippami = zeros(length(Wx[:,1]), length(Δ²kz[:,1]))
    R_NL_z_computation_turbo!(pippami, w1, Wx, Δ²kz)
    get_R_NL_z!(R_NL_z, z_list, R_test_interp, pippami)
    return R_NL_z
end
