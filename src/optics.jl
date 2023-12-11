module Optics

using Distributed, LinearAlgebra, SharedArrays
using ..HopTB
using ..HopTB.Utilities: constructmeshkpts, splitkpts
using ..HopTB.Parallel: ParallelFunction, claim!, stop!, parallel_do, parallel_sum
using ..HopTB.Meshes: UniformMesh

export getpermittivity, get_shift_cond, get_shg, cltberryconnection, cltshiftvector


@doc raw"""
```julia
get_generalized_dr(tm::AbstractTBModel, α::Int64, β::Int64, k::Vector{Float64})
    --> gdr::Matrix{ComplexF64}
```

Calculate generalized derivative of ``r``: ``r^α_{;β}``.

In the current gauge, generalized derivative of ``r`` coincides with derivative
of ``r``.
"""
function get_generalized_dr(tm::AbstractTBModel, α::Int64, β::Int64,
    k::Vector{Float64})::Matrix{ComplexF64}
    return getdr(tm, α, β, k)
end



################################################################################
##  Permittivity
################################################################################


@doc raw"""
```julia
get_permittivity_k!(
    χs::AbstractVector{ComplexF64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    ωs::AbstractArray{<:Real,1},
    μ::Float64,
    k::Vector{Float64};
    ϵ::Float64=0.1,
)
```

Calculate relative permittivity at `k` point and add the result to `χs`.
"""
function get_permittivity_k!(
    χs::AbstractVector{ComplexF64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    ωs::AbstractArray{<:Real,1},
    μ::Float64,
    k::Vector{Float64};
    ϵ::Float64=0.1,
)
    nωs = size(ωs, 1)
    rα = getA(tm, α, k)
    rβ = getA(tm, β, k)
    Es = geteig(tm, k).values
    for n in 1:tm.norbits, m in 1:tm.norbits
        En = Es[n]
        Em = Es[m]
        fn = (En < μ) ? 1 : 0
        fm = (Em < μ) ? 1 : 0
        if fn != fm
            tmp = 0.729494556 * (fn - fm) * (rα[n, m] * rβ[m, n])
            for iω in 1:nωs
                ω = ωs[iω] + im * ϵ
                χs[iω] += tmp / (Em - En - ω)
            end
        end
    end
    return nothing
end


@doc raw"""
```julia
getpermittivity(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    ωs::AbstractArray{<:Real,1},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
) --> Vector{ComplexF64} 
```

Calculate relative permittivity.

Permittivity is defined as

```math
χ^{αβ}=\frac{e^{2}}{\hbar}∫\frac{d^{3}\boldsymbol{k}}{(2π)^{3}}\sum_{n,m}
        \frac{f_{nm}r_{nm}^{α}r_{mn}^{β}}{ω_{mn}-ω-iϵ}.
```

This function returns ``χ^{αβ}/ϵ_0``.
"""
function getpermittivity(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    ωs::AbstractArray{<:Real,1},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
)::Vector{ComplexF64}
    nωs = size(ωs, 1)
    nks = prod(meshsize)
    mesh = UniformMesh(meshsize)
    χs = [SharedArray{ComplexF64}(nωs) for _ in 1:nprocs()]
    parallel_do(k -> get_permittivity_k!(χs[myid()], tm, α, β, ωs, μ, k; ϵ=ϵ), mesh, batchsize=batchsize)
    bzvol = abs(det(tm.rlat))
    return sum(χs) * bzvol / nks
end


################################################################################
##  Shift conductivity
################################################################################

@doc raw"""
```julia
get_shift_cond_k!(
    σs::AbstractVector{Float64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::Vector{Float64},
    μ::Float64,
    k::AbstractVector{Float64};
    ϵ::Float64=0.1
)
```

Calculate shift conductivity at `k` point and add the result to `σs`.
"""
function get_shift_cond_k!(
    σs::AbstractVector{Float64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::Vector{Float64},
    μ::Float64,
    k::AbstractVector{Float64};
    ϵ::Float64=0.1
)
    nωs = size(ωs, 1)
    Es = geteig(tm, k).values
    Aβ = getA(tm, β, k)
    Aγ = getA(tm, γ, k)
    gdrβα = get_generalized_dr(tm, β, α, k)
    gdrγα = get_generalized_dr(tm, γ, α, k)
    constant = -3.0828677458430857 * sqrt(1 / π) / ϵ # 10^6*(πe^3)/(ħ^2)/(2π)^3*ħ/e

    for n in 1:tm.norbits, m in 1:tm.norbits
        En = Es[n]
        Em = Es[m]
        fn = (En < μ) ? 1 : 0
        fm = (Em < μ) ? 1 : 0
        if fn != fm
            tmp = constant * (fn - fm) * imag(Aβ[m, n] * gdrγα[n, m] + Aγ[m, n] * gdrβα[n, m])
            for iω in 1:nωs
                ω = ωs[iω]
                σs[iω] += tmp * exp(-(En - Em - ω)^2 / ϵ^2)
            end
        end
    end

    return nothing 
end


function get_shift_cond_k(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::Vector{Float64},
    μ::Float64,
    k::AbstractVector{Float64};
    ϵ::Float64=0.1
)
    σs = zeros(Float64, length(ωs))
    get_shift_cond_k!(σs, tm, α, β, γ, ωs, μ, k; ϵ=ϵ)
    return σs
end


@doc raw"""
```julia
get_shift_cond(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64=β, ωs::Vector{Float64},
    μ::Float64, nkmesh::Vector{Int64}; ϵ::Float64=0.1, batchsize::Int64=1)
    --> Vector{Float64}
```

Calculate shift conductivity.

Shift conductivity is defined as

```math
σ^{αβγ}(ω)=
\frac{πe^3}{\hbar^2}∫\frac{d\boldsymbol{k}}{(2π)^3}\sum_{n,m}
f_{nm}I_{nm}^{αβγ}δ(ω_{nm}-ω),
```

where
```math
I_{nm}^{αβγ}=\Im[r_{mn}^β r_{nm;α}^γ+r_{mn}^γ r_{nm;α}^β].
```

This function returns shift conductivity in μA/V^2.
"""
function get_shift_cond(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::Vector{Float64},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
)
    nks = prod(meshsize)
    nωs = length(ωs)
    σs = [SharedArray{Float64}(nωs) for _ in 1:nprocs()]
    parallel_do(k -> get_shift_cond_k!(σs[myid()], tm, α, β, γ, ωs, μ, k; ϵ=ϵ),
        UniformMesh(meshsize), batchsize=batchsize)
    bzvol = abs(det(tm.rlat))
    return sum(σs) * bzvol / nks
end


get_shift_cond(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    ωs::Vector{Float64},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
) = get_shift_cond(tm, α, β, β, ωs, μ, meshsize; ϵ=ϵ, batchsize=batchsize)


################################################################################
##  Old stuff: collect something
################################################################################

@doc raw"""
Internal function for Berry connection parallel computing.
"""
function _cltberryconnection(atm::AbstractTBModel, α::Int64, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)
    berryconnection = zeros((atm.norbits, atm.norbits, nkpts))
    for ikpt in 1:nkpts
        k = kpts[:, ikpt]
        A = getA(atm, α, k)
        berryconnection[:, :, ikpt] .= abs2.(A)
    end
    return berryconnection
end


@doc raw"""
```julia
cltberryconnection(atm::AbstractTBModel, α::Int64, kpts::Matrix{Float64})
    --> berryconnection
```

Collect Berry connection for all k points.

The returned matrix berryconnection is
`berryconnection[n, m, ikpt] = |A[n, m]|^2` at `ikpt`.
"""
function cltberryconnection(atm::AbstractTBModel, α::Int64, kpts::AbstractMatrix{Float64})
    jobs = Vector{Future}()
    berryconnection = zeros((atm.norbits, atm.norbits, 0))
    kptslist = HopTB.Utilities.splitkpts(kpts, nworkers())
    for iw in 1:nworkers()
        job = @spawn _cltberryconnection(atm, α, kptslist[iw])
        append!(jobs, [job])
    end
    for iw in 1:nworkers()
        berryconnection = cat(berryconnection, HopTB.Utilities.safe_fetch(jobs[iw]), dims=(3,))
    end
    return berryconnection
end


"""
Internal function for shift vector parallel computing.
"""
function _cltshiftvector(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)
    shiftvector = zeros((atm.norbits, atm.norbits, nkpts))
    for ikpt in 1:nkpts
        k = kpts[:, ikpt]
        Aα = getA(atm, α, k)
        Aβα = getdr(atm, β, α, k)
        Aβ = getA(atm, β, k)
        for m in 1:atm.norbits
            for n in 1:atm.norbits
                if !(m == n)
                    shiftvector[n, m, ikpt] = imag(Aβα[m, n] / Aβ[m, n]) - real(Aα[m, m]) + real(Aα[n, n])
                end
            end
        end
    end
    return shiftvector
end


@doc raw"""
```julia
cltshiftvector(atm::AbstractTBModel, α::Int64, β::Int64, kpts::Matrix{Float64})
    --> shiftvector
```

Collect shift vector matrices for all k points.

`shiftvector[n, m, ikpt]` is the shift vector from band `m` to band `n` at `ikpt`.

Shift vector is defined to be
```math
R^{α,β}_{nm} = ∂_αϕ_{mn}^β-A_{mm}^α+A_{nn}^α.
```
"""
function cltshiftvector(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64})
    jobs = Vector{Future}()
    shiftvector = zeros((atm.norbits, atm.norbits, 0))
    kptslist = HopTB.Utilities.splitkpts(kpts, nworkers())
    for iw in 1:nworkers()
        job = @spawn _cltshiftvector(atm, α, β, kptslist[iw])
        append!(jobs, [job])
    end
    for iw in 1:nworkers()
        shiftvector = cat(shiftvector, HopTB.Utilities.safe_fetch(jobs[iw]), dims=(3,))
    end
    return shiftvector
end


################################################################################
##  Second harmonic generation
################################################################################

function _get_shg_inter_k_Es(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, k::Vector{Float64}, Es::Vector{Float64};
    ϵ::Float64=0.1)
    nωs = length(ωs)
    result = zeros(ComplexF64, nωs)
    Aα = getA(tm, α, k); Aβ = getA(tm, β, k); Aγ = getA(tm, γ, k)
    dEs = Es .- Es'
    fs = zeros(tm.norbits); fs[Es .< μ] .= 1.0; dfs = fs .- fs'
    ωs′ = ωs .+ im * ϵ
    tmp2 = zeros(ComplexF64, nωs, tm.norbits, tm.norbits)
    tmp3 = zeros(ComplexF64, nωs, tm.norbits, tm.norbits)
    for m in 1:tm.norbits, n in 1:tm.norbits, iω in 1:nωs
        ω = ωs′[iω]
        tmp2[iω, n, m] = 2 * dfs[n, m] / (dEs[m, n] - 2ω)
        tmp3[iω, n, m] = dfs[n, m] / (dEs[n, m] - ω)
    end
    for l in 1:tm.norbits, m in 1:tm.norbits, n in 1:tm.norbits
        if fs[n] == fs[m] == fs[l] continue end
        tmp1 = Aα[n, m] * (Aβ[m, l] * Aγ[l, n] + Aγ[m, l] * Aβ[l, n]) / (dEs[l, n] - dEs[m, l])
        for iω in 1:nωs
            result[iω] += tmp1 * (tmp2[iω, n, m] + tmp3[iω, l, n] + tmp3[iω, m, l])
        end
    end
    return result
end

function _get_shg_inter_k(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, k::Vector{Float64}; ϵ::Float64=0.1, scissor::Float64=0.0)
    Es = deepcopy(geteig(tm, k).values)
    Es[Es .> μ] .+= scissor
    return _get_shg_inter_k_Es(tm, α, β, γ, ωs, μ, k, Es, ϵ=ϵ)
end

function _get_shg_inter_k_B(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, B::Float64, Bdir::Int64, k::Vector{Float64};
    ϵ::Float64=0.1)
    Es = HopTB.Magnetism.get_field_modified_Es(tm, Bdir, B, k; double_degenerate=false)
    return _get_shg_inter_k_Es(tm, α, β, γ, ωs, μ, k, Es, ϵ=ϵ)
end

function _get_shg_mix_k_Es(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, k::Vector{Float64}, Es::Vector{Float64};
    ϵ::Float64=0.1)
    nωs = length(ωs)
    result = zeros(ComplexF64, nωs)
    Aα = getA(tm, α, k); Aβ = getA(tm, β, k); Aγ = getA(tm, γ, k)
    gdrβγ = get_generalized_dr(tm, β, γ, k); gdrγβ = get_generalized_dr(tm, γ, β, k);
    gdrαγ = get_generalized_dr(tm, α, γ, k); gdrαβ = get_generalized_dr(tm, α, β, k);
    gdrβα = get_generalized_dr(tm, β, α, k); gdrγα = get_generalized_dr(tm, γ, α, k);
    vβ = diag(getvelocity(tm, β, k)); vγ = diag(getvelocity(tm, γ, k));
    dEs = Es .- Es'
    fs = zeros(tm.norbits); fs[Es .< μ] .= 1.0; dfs = fs .- fs'
    ωs′ = ωs .+ im * ϵ
    for n in 1:tm.norbits, m in 1:tm.norbits
        if fs[n] == fs[m] continue end
        tmp1 = 2 * Aα[n, m] * (gdrβγ[m, n] + gdrγβ[m, n])
        tmp2 = gdrαγ[n, m] * Aβ[m, n] + gdrαβ[n, m] * Aγ[m, n]
        tmp3 = Aα[n, m] * (Aβ[m, n] * (vγ[m] - vγ[n]) + Aγ[m, n] * (vβ[m] - vβ[n]))
        tmp4 = gdrβα[n, m] * Aγ[m, n] + gdrγα[n, m] * Aβ[m, n]
        for iω in 1:nωs
            ω = ωs′[iω]
            result[iω] += im * dfs[n, m] * (
                tmp1 / (dEs[m, n] * (dEs[m, n] - 2ω)) +
                tmp2 / (dEs[m, n] * (dEs[m, n] - ω)) +
                tmp3 / dEs[m, n]^2 * (1 / (dEs[m, n] - ω) - 4 / (dEs[m, n] - 2ω)) -
                tmp4 / (2dEs[m, n] * (dEs[m, n] - ω))
            )
        end
    end
    return result
end

function _get_shg_mix_k(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, k::Vector{Float64}; ϵ::Float64=0.1, scissor::Float64=0.0)
    Es = deepcopy(geteig(tm, k).values)
    Es[Es .> μ] .+= scissor
    return _get_shg_mix_k_Es(tm, α, β, γ, ωs, μ, k, Es, ϵ=ϵ)
end

function _get_shg_mix_k_B(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, B::Float64, Bdir::Int64, k::Vector{Float64};
    ϵ::Float64=0.1)
    Es = HopTB.Magnetism.get_field_modified_Es(tm, Bdir, B, k; double_degenerate=false)
    return _get_shg_mix_k_Es(tm, α, β, γ, ωs, μ, k, Es, ϵ=ϵ)
end

function _get_shg_worker(k::Vector{Float64}, tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, ϵ::Float64, scissor::Float64)
    return _get_shg_inter_k(tm, α, β, γ, ωs, μ, k, ϵ=ϵ, scissor=scissor) + _get_shg_mix_k(tm, α, β, γ, ωs, μ, k, ϵ=ϵ, scissor=scissor)
end

function _get_shg_B_worker(k::Vector{Float64}, tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, B::Float64, Bdir::Int64, ϵ::Float64)
    return _get_shg_inter_k_B(tm, α, β, γ, ωs, μ, B, Bdir, k, ϵ=ϵ) + _get_shg_mix_k_B(tm, α, β, γ, ωs, μ, B, Bdir, k, ϵ=ϵ)
end

"""
```julia
get_shg(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64, ωs::Vector{Float64},
    μ::Float64, kmesh::Vector{Int64}; ϵ::Float64=0.1, scissor::Float64=0.0)
    --> σs::Vector{ComplexF64}
```

Calculate second harmonic generation.

The expression of second harmonic generation is in [Wang et al 2017].

The unit of σs is pm/V.
"""
function get_shg(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, nkmesh::Vector{Int64};
    ϵ::Float64=0.1, scissor::Float64=0.0)
    nks = prod(nkmesh); nωs = size(ωs, 1)
    ks = constructmeshkpts(nkmesh)
    pf = ParallelFunction(_get_shg_worker, tm, α, β, γ, ωs, μ, ϵ, scissor, len=10 * nworkers())
    σs = zeros(ComplexF64, nωs)
    @async for i in 1:nks
        pf(ks[:, i])
    end
    @sync @async for i in 1:nks
        σs += claim!(pf)
    end
    stop!(pf)
    bzvol = abs(det(tm.rlat))
    return σs * bzvol / nks * (-36.474728077745624) # -e/ϵ0*1e12/(2π)^3
end


"""
```julia
get_shg_B(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64, ωs::Vector{Float64},
    μ::Float64, B::Float64, Bdir::Int64, nkmesh::Vector{Int64}; ϵ::Float64=0.1)
    --> σs::Vector{ComplexF64}
```

Calculate second harmonic generation with Zeeman interaction.

The expression of second harmonic generation is in [Wang et al 2017].

The unit of σs is pm/V.

!!! warning
    This function has not been tested and I have not understand whether the expression
    in [Wang et al 2017] can be used for TRS-broken systems.
"""
function get_shg_B(tm::AbstractTBModel, α::Int64, β::Int64, γ::Int64,
    ωs::Vector{Float64}, μ::Float64, B::Float64, Bdir::Int64, nkmesh::Vector{Int64}; ϵ::Float64=0.1)
    nks = prod(nkmesh); nωs = size(ωs, 1)
    ks = constructmeshkpts(nkmesh)
    pf = ParallelFunction(_get_shg_B_worker, tm, α, β, γ, ωs, μ, B, Bdir, ϵ, len=10 * nworkers())
    σs = zeros(ComplexF64, nωs)

    @async for i in 1:nks
        pf(ks[:, i])
    end
    @sync @async for i in 1:nks
        σs += claim!(pf)
    end

    stop!(pf)
    bzvol = abs(det(tm.rlat))
    return σs * bzvol / nks * (-36.474728077745624) # -e/ϵ0*1e12/(2π)^3
end


################################################################################
##  Drude weight
################################################################################

@doc raw"""
```julia
get_Drude_weight_k(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    μ::Float64,
    k::AbstractVector{Float64};
    kBT::Float64=0.01
)  --> ComplexF64
```

This function calculates Drude weight contribution from point `k`.
"""
function get_Drude_weight_k(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    μ::Float64,
    k::AbstractVector{Float64};
    kBT::Float64=0.01
)
    result = 0.0im
    Es = geteig(tm, k).values
    vα = getvelocity(tm, α, k)
    vβ = getvelocity(tm, β, k)
    for n in 1:tm.norbits
        tmp = (Es[n] - μ) / kBT
        result += vα[n, n] * vβ[n, n] / (1 + cosh(tmp)) / (2kBT)
    end
    return result * 98.13072812989917 # e^2/hbar*10^8/(2π)^3
end


@doc raw"""
```julia
function get_Drude_weight(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    μ::Float64,
    meshsize::AbstractVector{Int64};
    kBT::Float64=0.01,
    batchsize::Int64=1
)
```

This function returns Drude weight in eV/(Ω⋅cm).
"""
function get_Drude_weight(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    μ::Float64,
    meshsize::AbstractVector{Int64};
    kBT::Float64=0.01,
    batchsize::Int64=1
)
    nks = prod(meshsize)
    mesh = UniformMesh(meshsize)
    result = parallel_sum(
        k -> get_Drude_weight_k(tm, α, β, μ, k; kBT=kBT),
        mesh,
        0.0im;
        batchsize=batchsize
    )
    bzvol = abs(det(tm.rlat))
    return result * bzvol / nks
end


@doc raw"""
```julia
get_injection_conductivity_k!(
    σs::AbstractVector{ComplexF64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::AbstractVector{Float64},
    μ::Float64,
    k::Vector{Float64};
    ϵ::Float64=0.1
)
```

Calculate injection conductivity at `k` point and add the result to `σs`.
"""
function get_injection_conductivity_k!(
    σs::AbstractVector{ComplexF64},
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::AbstractVector{Float64},
    μ::Float64,
    k::Vector{Float64};
    ϵ::Float64=0.1
)
    nωs = length(ωs)
    vα = getvelocity(tm, α, k)
    rβ = getA(tm, β, k)
    rγ = getA(tm, γ, k)
    Es = geteig(tm, k).values
    for n in 1:tm.norbits, m in 1:tm.norbits
        En = Es[n]
        Em = Es[m]
        if abs(En - Em) > only(HopTB.DEGEN_THRESH)
            fn = (En < μ) ? 1 : 0
            fm = (Em < μ) ? 1 : 0
            for iω in 1:nωs
                ω = ωs[iω]
                delta = exp(-(Em - En - ω)^2 / ϵ^2) / ϵ / √π
                σs[iω] += (vα[n, n] - vα[m, m]) * rβ[n, m] * rγ[m, n] * (fn - fm) * delta * 3.082867745843085 # π * e^2 / (ħ * (2π)^3)
            end
        end
    end
    return nothing
end

function get_injection_conductivity_k(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::AbstractVector{Float64},
    μ::Float64,
    k::Vector{Float64};
    ϵ::Float64=0.1
)
    σs = zeros(ComplexF64, length(ωs))
    get_injection_conductivity_k!(σs, tm, α, β, γ, ωs, μ, k; ϵ=ϵ)
    return σs
end

@doc raw"""
```julia
get_injection_conductivity(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::AbstractVector{Float64},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
)
```

This function returns injection conductivity in μA * eV / V^2.
"""
function get_injection_conductivity(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    ωs::AbstractVector{Float64},
    μ::Float64,
    meshsize::Vector{Int64};
    ϵ::Float64=0.1,
    batchsize::Int64=1
)
    nks = prod(meshsize)
    nωs = length(ωs)
    σs = [SharedArray{ComplexF64}(nωs) for _ in 1:nprocs()]
    parallel_do(k -> get_injection_conductivity_k!(σs[myid()], tm, α, β, γ, ωs, μ, k; ϵ=ϵ),
        UniformMesh(meshsize), batchsize=batchsize)
    bzvol = abs(det(tm.rlat))
    return sum(σs) * bzvol / nks
end

end
