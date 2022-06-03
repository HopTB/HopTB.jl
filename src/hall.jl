module Hall

using LinearAlgebra, Distributed
using ..Hop
using ..Hop.Utilities: fermidirac, constructmeshkpts, splitkpts
using ..Hop.Parallel: ParallelFunction, claim!, stop!, parallel_sum

export getahc


function _getahc(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64};
    Ts::Vector{Float64} = [0.0], μs::Vector{Float64} = [0.0])
    nkpts = size(kpts, 2)
    itgrd = zeros(ComplexF64, length(Ts), length(μs))
    for ik in 1:nkpts
        k = kpts[:, ik]
        egvals, egvecs = geteig(atm, k)
        order = [0, 0, 0]
        order[α] = 1
        Sbar_α = egvecs' * getdS(atm, Tuple(order), k) * egvecs
        order = [0, 0, 0]
        order[β] = 1
        Sbar_β = egvecs' * getdS(atm, Tuple(order), k) * egvecs
        Abar_α = egvecs' * Hop.getAw(atm, α, k) * egvecs
        Abar_β = egvecs' * Hop.getAw(atm, β, k) * egvecs
        Dα = Hop.getD(atm, α, k)
        Dβ = Hop.getD(atm, β, k)
        order = [0, 0, 0]
        order[α] = 1
        dAw_βα = Hop.getdAw(atm, β, Tuple(order), k)
        order = [0, 0, 0]
        order[β] = 1
        dAw_αβ = Hop.getdAw(atm, α, Tuple(order), k)
        Ωbar_αβ = egvecs' * (dAw_βα - dAw_αβ) * egvecs
        tmp1 = Sbar_α*Abar_β
        tmp2 = Sbar_β*Abar_α
        for iT in 1:length(Ts)
            for iμ in 1:length(μs)
                for n in 1:atm.norbits
                    f = fermidirac(Ts[iT], egvals[n]-μs[iμ])
                    itgrd[iT, iμ] +=  f*(Ωbar_αβ[n, n]-tmp1[n, n]+tmp2[n, n])
                end
                for n in 1:atm.norbits, m in 1:atm.norbits
                    fm = fermidirac(Ts[iT], egvals[m] - μs[iμ])
                    fn = fermidirac(Ts[iT], egvals[n] - μs[iμ])
                    itgrd[iT, iμ] += (fm - fn) * (im * Dα[n, m] * Dβ[m, n] + Dα[n, m] * Abar_β[m, n] - Dβ[n, m] * Abar_α[m, n])
                end
            end
        end
    end
    return real.(itgrd)
end


@doc raw"""
```julia
getahc(atm::AbstractTBModel, α::Int64, β::Int64, nkmesh::Vector{Int64};
    Ts::Vector{Float64}=[0.0], μs::Vector{Float64}=[0.0])::Matrix{Float64}
```

Calculate anomalous Hall conductivity $σ^{αβ}$.

Anomalous Hall conductivity is defined by
```math
σ^{αβ}=-\frac{e^2}{ħ}\int\frac{d\boldsymbol{k}}{(2pi)^3}f_nΩ_{nn}^{αβ}.
```

The returned matrix $σ^{αβ}[m, n]$ is AHC for temperature Ts[m] and
chemical potential μs[n].

The returned AHC is in unit (Ω⋅cm)^-1.
"""
function getahc(atm::AbstractTBModel, α::Int64, β::Int64, nkmesh::Vector{Int64};
    Ts::Vector{Float64} = [0.0], μs::Vector{Float64} = [0.0])
    @assert size(nkmesh, 1) == 3
    nkpts = prod(nkmesh)
    kpts = Hop.Utilities.constructmeshkpts(nkmesh)
    kptslist = Hop.Utilities.splitkpts(kpts, nworkers())

    jobs = Vector{Future}()
    for iw in 1:nworkers()
        job = @spawn _getahc(atm, α, β, kptslist[iw]; Ts = Ts, μs = μs)
        append!(jobs, [job])
    end


    σs = zeros(Float64, length(Ts), length(μs))
    for iw in 1:nworkers()
        σs += Hop.Utilities.safe_fetch(jobs[iw])
    end

    bzvol = abs(dot(cross(atm.rlat[:, 1], atm.rlat[:, 2]), atm.rlat[:, 3]))
    return σs * bzvol / nkpts * (-98.130728142) # -e**2/(hbar*(2pi)^3)*1.0e10/100
end


function _collect_berry_curvature(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)
    berry_curvature = zeros(atm.norbits, nkpts)
    for ik in 1:nkpts
        k = kpts[:, ik]
        egvals, egvecs = geteig(atm, k)
        order = [0, 0, 0]; order[α] = 1; Sbar_α = egvecs' * getdS(atm, Tuple(order), k) * egvecs
        order = [0, 0, 0]; order[β] = 1; Sbar_β = egvecs' * getdS(atm, Tuple(order), k) * egvecs
        Abar_α = egvecs' * Hop.getAw(atm, α, k) * egvecs
        Abar_β = egvecs' * Hop.getAw(atm, β, k) * egvecs
        Dα = Hop.getD(atm, α, k)
        Dβ = Hop.getD(atm, β, k)
        order = [0, 0, 0]; order[α] = 1; dAw_βα = Hop.getdAw(atm, β, Tuple(order), k)
        order = [0, 0, 0]; order[β] = 1; dAw_αβ = Hop.getdAw(atm, α, Tuple(order), k)
        Ωbar_αβ = egvecs' * (dAw_βα - dAw_αβ) * egvecs
        berry_curvature[:, ik] = real.(diag(Ωbar_αβ - Sbar_α * Abar_β + Sbar_β * Abar_α - im * Dα * Dβ + 
            im * Dβ * Dα - Dα * Abar_β + Abar_β * Dα + Dβ * Abar_α - Abar_α * Dβ))
    end
    return berry_curvature
end

"""
```julia
collect_berry_curvature(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64})::Matrix{Float64}
```

Collect berry curvature.

Standard units is used (eV and Å).

The returned matrix Ω[n, ik] is berry curvature for band n at ik point.
"""
function collect_berry_curvature(atm::AbstractTBModel, α::Int64, β::Int64, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)

    kptslist = Hop.Utilities.splitkpts(kpts, nworkers())
    jobs = Vector{Future}()
    for iw in 1:nworkers()
        job = @spawn _collect_berry_curvature(atm, α, β, kptslist[iw])
        append!(jobs, [job])
    end

    result = zeros((atm.norbits, 0))
    for iw in 1:nworkers()
        result = cat(result, Hop.Utilities.safe_fetch(jobs[iw]), dims = (2,))
    end
    return result
end


function shc_worker(kpts::Matrix{Float64}, tm::TBModel, α::Int64, β::Int64, γ::Int64,
    Ts::Vector{Float64}, μs::Vector{Float64}, ϵ::Float64)
    result = zeros(length(Ts), length(μs))
    nkpts = size(kpts, 2)
    for ik in 1:nkpts
        k = kpts[:, ik]
        vα = getvelocity(tm, α, k); vβ = getvelocity(tm, β, k)
        sγ = getspin(tm, γ, k)
        jαγ = (sγ*vα+vα*sγ)/2
        egvals, _ = geteig(tm, k)
        for (iT, T) in enumerate(Ts), (iμ, μ) in enumerate(μs)
            for n in 1:tm.norbits
                ϵn = egvals[n]
                fn = fermidirac(T, ϵn-μ)
                for m in 1:tm.norbits
                    ϵm = egvals[m]
                    result[iT, iμ] += -fn*imag(jαγ[n, m]*vβ[m, n]/((ϵn-ϵm)^2+ϵ^2))
                end
            end
        end
    end
    return result
end

@doc raw"""
```julia
getshc(tm::TBModel, α::Int64, β::Int64, γ::Int64, nkmesh::Vector{Int64};
    Ts::Vector{Float64}=[0.0], μs::Vector{Float64}=[0.0], ϵ::Float64=0.1)::Matrix{Float64}
```

Calculate spin Hall conductivity for different temperature (`Ts`, first dimension)
and chemical potential (`μs`, second dimension).

Spin Hall conductivity is defined as
```math
σ_{αβ}^{γ} = eħ\int\frac{d^3 \boldsymbol{k}}{(2π)^3}\sum_n f_n Ω^{γ}_{n,αβ},
```
where the spin Berry curvature is
```math
Ω_{n,αβ}^{γ} = -2 \text{Im} [\sum_{m≠n} \frac{⟨n|\hat{j}_α^γ|m⟩⟨m|\hat{v}_β|n⟩}{(ϵ_n-ϵ_m)^2+ϵ^2}]
```
and the spin current operator is
```math
\hat{j}_α^γ = \frac{1}{2} \{\hat{v}_a, \hat{s}_c\}.
```

Spin Hall conductivity from this function is in ħ/e (Ω*cm)^-1.
"""
function getshc(tm::TBModel, α::Int64, β::Int64, γ::Int64, nkmesh::Vector{Int64};
    Ts::Vector{Float64}=[0.0], μs::Vector{Float64}=[0.0], ϵ::Float64=0.1)
    size(nkmesh, 1) == 3 || error("nkmesh should be a 3-element vector.")
    nkpts = prod(nkmesh)
    kptslist = splitkpts(constructmeshkpts(nkmesh), nworkers())
    pf = ParallelFunction(shc_worker, tm, α, β, γ, Ts, μs, ϵ)

    for iw in 1:nworkers()
        pf(kptslist[iw])
    end
    σs = zeros(Float64, length(Ts), length(μs))
    for iw in 1:nworkers()
        σs += claim!(pf)
    end
    stop!(pf)
    bzvol = abs((tm.rlat[:, 1]×tm.rlat[:, 2])⋅tm.rlat[:, 3])
    return σs*bzvol/nkpts*98.130728142 # e**2/(hbar*(2π)^3)*1.0e10/100
end


################################################################################
##  Berry curvature dipole
################################################################################


@doc raw"""
```
get_berry_curvature_dipole(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    fss::Vector{FermiSurface}
)::Float64
```

Calculate
```math
\sum_n \int_{\text{FS}_n} \frac{d \sigma}{(2\pi)^3} \Omega_{n}^{\alpha \beta} \frac{v_n^{\gamma}}{|\boldsymbol{v}_n|}
```
which is related to the Berry curvature dipole contribution to the second order photocurrent.
"""
function get_berry_curvature_dipole(
    tm::AbstractTBModel,
    α::Int64,
    β::Int64,
    γ::Int64,
    fss::Vector{FermiSurface}
)
    result = 0.0
    for fs in fss
        result += parallel_sum(ik -> begin
            k = fs.ks[:, ik]
            Ω = get_berry_curvature(tm, α, β, k)[fs.bandidx]
            v = real([getvelocity(tm, i, k)[fs.bandidx, fs.bandidx] for i in 1:3])
            Ω * v[γ] * fs.weights[ik] / norm(v)
        end, 1:size(fs.ks, 2), 0.0)
    end
    return result / (2π)^3
end

end
