module Topology

using LinearAlgebra
using ..Hop

export get_smooth_gauge


function parallel_transport(atm::AbstractTBModel, getU::Function,
    Ustart::Matrix{ComplexF64}, kpaths::AbstractMatrix{<:Real}, ndiv::Int64)
    kpts = Hop.Utilities.constructlinekpts(kpaths, ndiv)
    result = Ustart
    for i in 1:(size(kpts, 2) - 1)
        k1, k2 = kpts[:, i], kpts[:, i + 1]
        δk = k2 - k1
        S = getS(atm, k1)
        U = getU(k2)
        Awδk = sum(map(x->Hop.getAw(atm, x, k1), 1:3) .* (atm.rlat * δk))
        tmp = -im * Awδk + S
        result = U * (U') * (tmp') * result
    end
    return result
end


function get_wilson_spectrum(atm::AbstractTBModel, getU::Function,
    kpaths::AbstractMatrix{<:Real}, ndiv::Int64)
    kpts = Hop.Utilities.constructlinekpts(kpaths, ndiv)
    Ustart = getU(kpts[:, 1])
    W = Ustart' * getS(atm, kpts[:, end]) * parallel_transport(atm, getU, Ustart, kpaths, ndiv)
    tmp = eigvals(W)
    err = maximum(abs.(tmp) .- 1.0)
    if  err > 0.01
        @warn "W is not unitary!, error: $err"
    end
    return sort(imag(log.(tmp)))
end


function get_wilson_spectrum(atm::AbstractTBModel, bandind::Vector{Int64},
    kpaths::AbstractMatrix{<:Real}, ndiv::Int64)
    function getU(k)
        return geteig(atm, k).vectors[:, bandind]
    end
    return get_wilson_spectrum(atm, getU, kpaths, ndiv)
end


"""
```julia
get_smooth_gauge(tm::TBModel, getv::Function, v0::Vector{<:Number},
    kpts::Matrix{<:Real}; unwind::Bool=false)::Matrix{ComplexF64}
```

Parallel transport the Bloch function produced by `getv` at a line in BZ denoted by `kpts`.
The start point is `kpts[:, 1]` where the vector is assumed to be `v0`. If unwind is true, then
the first k point in `kpts` and the last k point in `kpts` are assumed to be
identical k point. In this case, the vector on the two k points are made
the same. kpts are assumed to be even spaced if unwind is true.

The output is a matrix of size (length(v0), length(kpts)).
"""
function get_smooth_gauge(tm::TBModel, getv::Function, v0::Vector{<:Number},
    kpts::Matrix{<:Real}; unwind::Bool=false)
    size(kpts, 1) == 3 || throw(ArgumentError("kpts is in wrong size."))
    size(kpts, 2) > 0 || throw(ArgumentError("At least one kpt is needed."))
    size(getv(kpts[:, 1])) == size(v0) || throw(ArgumentError("getv should return a vector with the same shape as v0."))
    nkpts = size(kpts, 2); vlen = length(v0)
    result = zeros(ComplexF64, vlen, nkpts)
    result[:, 1] = v0
    for i in 1:(nkpts-1)
        k1, k2 = kpts[:, i], kpts[:, i+1]
        δk = k2-k1
        S1 = getS(tm, k1); S2 = getS(tm, k2)
        v = getv(k2)
        Awδk = sum(map(x->Hop.getAw(tm, x, k1), 1:3).*(tm.rlat*δk))
        tmp = v*(v')*((-im*Awδk+S1)')*result[:, i]
        result[:, i+1] = tmp/√(tmp'*S2*tmp)
    end
    if unwind
        Δk = kpts[:, end]-kpts[:, 1]
        maximum(Δk-round.(Δk)) < 1.0e-3 || throw(ArgumentError("The start and end points do not seem to be the same."))
        max_ind = findmax(abs.(v0))[2]
        θ = angle(result[max_ind, end]/result[max_ind, 1])
        for i in 1:nkpts
            result[:, i] *= exp(-im*θ*(i-1)/(nkpts-1))
        end
    end
    return result
end

"""
```julia
get_smooth_gauge(tm::TBModel, getv::Function, nkmesh::Vector{Int64})
```

Find a smooth gauge in the BZ for the Bloch functions returned by `getv`.
Please make sure Chern number is 0 for any possible submanifold.
"""
function get_smooth_gauge(tm::TBModel, getv::Function, nkmesh::Vector{Int64})
    length(nkmesh) == 3 || throw(ArgumentError("nkmesh should have length 3."))
    v0 = getv([0.0, 0.0, 0.0])
    length(size(v0)) == 1 || throw(ArgumentError("getU should return a vector."))
    result = zeros(ComplexF64, length(v0), nkmesh...)
    nkmesh′ = nkmesh .+ 1
    kpts = zeros(Float64, 3, nkmesh′...)
    for ikx in 1:nkmesh′[1], iky in 1:nkmesh′[2], ikz in 1:nkmesh′[3]
        kpts[:, ikx, iky, ikz] = [(ikx-1)/nkmesh[1], (iky-1)/nkmesh[2], (ikz-1)/nkmesh[3]]
    end

    # first step
    result[:, :, 1, 1] = get_smooth_gauge(tm, getv, v0, kpts[:, :, 1, 1], unwind=true)[:, 1:end-1]

    # second step
    for i in 1:nkmesh[1]
        result[:, i, :, 1] = get_smooth_gauge(tm, getv, result[:, i, 1, 1], kpts[:, i, :, 1], unwind=true)[:, 1:end-1]
    end

    # third step
    for i in 1:nkmesh[1], j in 1:nkmesh[2]
        result[:, i, j, :] = get_smooth_gauge(tm, getv, result[:, i, j, 1], kpts[:, i, j, :], unwind=true)[:, 1:end-1]
    end

    return result
end

"""
```julia
get_smooth_gauge(tm::TBModel, bandind::Int64, nkmesh::Vector{Int64})
```

Find a smooth gauge in the BZ for the `bandind` band.
Please make sure Chern number is 0 for any possible submanifold.
"""
function get_smooth_gauge(tm::TBModel, bandind::Int64, nkmesh::Vector{Int64})
    function getv(k)
        return geteig(tm, k).vectors[:, bandind]
    end
    return get_smooth_gauge(tm, getv, nkmesh)
end

end
