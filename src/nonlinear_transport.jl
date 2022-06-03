module NonlinearTransport

using LinearAlgebra
using ..HopTB
using ..HopTB.Parallel: parallel_sum

################################################################################
##  Second order intrinsic nonlinear conductivity
################################################################################

function get_2nd_intrinsic_conductivity_k(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    bandidx::Integer,
    k::AbstractVector{<:Real};
    is_doubly_degenerate::Bool=false
)
    (is_doubly_degenerate && iseven(bandidx)) && error("wrong bandidx.")
    factor = -1.9626145625979834 # -2 e^2 1.0e6 / ħ / (2π)^3
    v = [getvelocity(tm, dir, k)[bandidx, bandidx] for dir in 1:3]
    Es = geteig(tm, k).values
    Aα = getA(tm, α, k)
    Aβ = getA(tm, β, k)
    Aγ = getA(tm, γ, k)
    if is_doubly_degenerate
        degenerate_band_indices = [bandidx, bandidx + 1]
    else
        degenerate_band_indices = [bandidx]
    end
    result = 0.0
    for n in degenerate_band_indices
        for m in 1:tm.norbits
            if abs(Es[n] - Es[m]) > HopTB.DEGEN_THRESH[1] && !(m in degenerate_band_indices)
                result += real((v[α] * Aβ[n, m] * Aγ[m, n] - v[β] * Aα[n, m] * Aγ[m, n]) / (Es[n] - Es[m]))
            end
        end
    end
    return factor * result / norm(v)
end


@doc raw"""
```julia
get_2nd_intrinsic_conductivity(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    fss::Vector{FermiSurface};
    is_doubly_degenerate::Bool=false
)
```

Calculate intrinsic nonlinear conductivity σ^{αβγ}. `fss` is the Fermi surfaces of the system.

If every band is doubly degenerate (for example, due to space-time inversion symmetry), `is_doubly_degenerate`
should be set to `true`. In this case, `fss` should only contain bands with odd band indices.

The function returns conductivity in μA/V^2.
"""
function get_2nd_intrinsic_conductivity(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    fss::Vector{FermiSurface};
    is_doubly_degenerate::Bool=false,
    batchsize::Int64=1
)
    result = 0.0
    for fs in fss
        result += parallel_sum(
            ik -> get_2nd_intrinsic_conductivity_k(
                tm, α, β, γ, fs.bandidx, fs.ks[:, ik], is_doubly_degenerate=is_doubly_degenerate
            ) * fs.weights[ik],
            1:size(fs.ks, 2),
            0.0;
            batchsize=batchsize
        )
    end
    return result
end


################################################################################
##  Second order drude weight
################################################################################


function get_2nd_drude_weight_k(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    bandidx::Integer,
    k::AbstractVector{<:Real};
    is_doubly_degenerate::Bool=false
)
    (is_doubly_degenerate && iseven(bandidx)) && error("wrong bandidx.")
    factor = -0.9813072812989917 # -e^2 1.0e6 / ħ / (2π)^3
    v = [getdEs(tm, dir, k)[bandidx] for dir in 1:3]
    return factor * v[α] * getdEs(tm, β, γ, k)[bandidx] / norm(v) * (is_doubly_degenerate ? 2.0 : 1.0)
end

@doc raw"""
```julia
get_2nd_drude_weight(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    fss::Vector{FermiSurface};
    is_doubly_degenerate::Bool=false,
    batchsize::Int64=1
)
```

Calculate second order Drude weight Λ^{αβγ}. `fss` is the Fermi surfaces of the system.

If every band is doubly degenerate (for example, due to space-time inversion symmetry), `is_doubly_degenerate`
should be set to `true`. In this case, `fss` should only contain bands with odd band indices.

The function returns conductivity in eV^2 μA / V^2.
"""
function get_2nd_drude_weight(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    γ::Integer,
    fss::Vector{FermiSurface};
    is_doubly_degenerate::Bool=false,
    batchsize::Int64=1
)
    result = 0.0
    for fs in fss
        result += parallel_sum(
            ik -> get_2nd_drude_weight_k(
                tm, α, β, γ, fs.bandidx, fs.ks[:, ik], is_doubly_degenerate=is_doubly_degenerate
            ) * fs.weights[ik],
            1:size(fs.ks, 2),
            0.0;
            batchsize=batchsize
        )
    end
    return result
end

end