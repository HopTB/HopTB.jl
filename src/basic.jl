using LinearAlgebra
using .Memoize

export getH, getdH, getS, getdS, geteig, getdEs, getA, getdr, getvelocity, getspin, get_berry_curvature

const DEGEN_THRESH = [1.0e-4]

function set_degen_thresh(val)
    @warn "DEGEN_THRESH should be set before any calculations."
    DEGEN_THRESH[1] = val
    return nothing
end

const σ1 = [0 1; 1 0]
const σ2 = [0 -im; im 0]
const σ3 = [1 0; 0 -1]
const σs = [σ1, σ2, σ3]

@doc raw"""
```julia
getdH(tm::AbstractTBModel, order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate `order` derivative of Hamiltonian.
"""
@memoize k function getdH(tm::TBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dH = zeros(ComplexF64, tm.norbits, tm.norbits)
    for (R, hopping) in tm.hoppings
        Rc = tm.lat*R # R in Cartesian coordinate
        phase = exp(im*2π*(k⋅R))
        coeff = (im*Rc[1])^(order[1])*(im*Rc[2])^(order[2])*(im*Rc[3])^(order[3])*phase
        @. dH += coeff*hopping
    end
    return dH
end

@memoize k function getdH(
    sm::SharedTBModel,
    order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real}
)::Matrix{ComplexF64}
    nhalfRs = (size(sm.Rs, 2) + 1) ÷ 2
    norbits = sm.norbits
    coeffs = im^sum(order) * prod(sm.Rcs.^order; dims=1) .* exp.(im * 2π * (k' * sm.Rs))
    tmp = reshape(sm.H * (coeffs[1, 1:nhalfRs]), (norbits, norbits))
    return tmp' + tmp
end


function _getdS(nm::TBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dS = zeros(ComplexF64, nm.norbits, nm.norbits)
    for (R, overlap) in nm.overlaps
        Rc = nm.lat*R # R in Cartesian coordinate
        phase = exp(im*2π*(k⋅R))
        coeff = (im*Rc[1])^(order[1])* (im*Rc[2])^(order[2])*(im*Rc[3])^(order[3])*phase
        @. dS += coeff*overlap
    end
    return dS
end

@doc raw"""
```julia
getdS(tm::AbstractTBModel, order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate `order` derivative of overlap.
"""
@memoize k function getdS(tm::TBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return tm.isorthogonal ? (order == (0, 0, 0) ? I(tm.norbits) : 0I(tm.norbits)) : _getdS(tm, order, k)
end

function _getdS(
    sm::SharedTBModel,
    order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real}
)::Matrix{ComplexF64}
    nhalfRs = (size(sm.Rs, 2) + 1) ÷ 2
    norbits = sm.norbits
    coeffs = im^sum(order) * prod(sm.Rcs.^order; dims=1) .* exp.(im * 2π * (k' * sm.Rs))
    tmp = reshape(sm.S * (coeffs[1, 1:nhalfRs]), (norbits, norbits))
    return tmp' + tmp
end

@memoize k function getdS(sm::SharedTBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return sm.isorthogonal ? (order == (0, 0, 0) ? I(sm.norbits) : 0I(sm.norbits)) : _getdS(sm, order, k)
end


@doc raw"""
```julia
getdAw(tm::AbstractTBModel, α::Int64, order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate `order` derivative of ``i⟨u_n^{(W)}|∂_{k_α}u_m^{(W)}⟩``.
"""
@memoize k function getdAw(tm::TBModel, α::Int64, order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dAw = zeros(ComplexF64, tm.norbits, tm.norbits)
    for (R, pos) in tm.positions
        Rc = tm.lat*R
        phase = exp(im*2π*(k⋅R))
        coeff = (im*Rc[1])^(order[1])*(im*Rc[2])^(order[2])*(im*Rc[3])^(order[3])*phase
        @. dAw += coeff*pos[α]
    end
    return dAw'
end

@memoize k function getdAw(
    sm::SharedTBModel,
    α::Int64,
    order::Tuple{Int64,Int64,Int64},
    k::AbstractVector{<:Real}
)::Matrix{ComplexF64}
    nhalfRs = (size(sm.Rs, 2) + 1) ÷ 2
    norbits = sm.norbits
    coeffs = im^sum(order) * prod(sm.Rcs.^order; dims=1) .* exp.(im * 2π * (k' * sm.Rs))
    dAw = reshape(sm.r[α] * (coeffs[1, :]), (norbits, norbits))
    return dAw'
end


@doc raw"""
```julia
getH(tm::AbstractTBModel, k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate Hamiltonian at a reduced `k` point.
"""
function getH(tm::AbstractTBModel, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return getdH(tm, (0, 0, 0), k)
end


@doc raw"""
```julia
getS(tm::AbstractTBModel, k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate overlap matrix at a reduced `k` point.
"""
function getS(tm::AbstractTBModel, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return getdS(tm, (0, 0, 0), k)
end


"""
HermEig wraps eigenvalues and eigenvectors of a Hermitian eigenvalue problem.

# Fields
 - `values::Vector{Float64}`: eigenvalues
 - `vectors::Matrix{ComplexF64}`: eigenvectors stored in column
"""
struct HermEig
    values::Vector{Float64}
    vectors::Matrix{ComplexF64}
end

Base.iterate(S::HermEig) = (S.values, Val(:vectors))
Base.iterate(S::HermEig, ::Val{:vectors}) = (S.vectors, Val(:done))
Base.iterate(S::HermEig, ::Val{:done}) = nothing

@doc raw"""
```julia
geteig(tm::AbstractTBModel, k::AbstractVector{<:Real})::HermEig
```

Calculate eigenvalues and eigenvectors of `tm` at a reduced `k` point.
"""
@memoize k function geteig(tm::AbstractTBModel, k::AbstractVector{<:Real})::HermEig
    H = getH(tm, k)
    if tm.isorthogonal
        (Es, V) = eigen(Hermitian(H))
    else
        S = getS(tm, k)
        (Es, V) = eigen(Hermitian(H), Hermitian(S))
    end
    return HermEig(Es, V)
end


@doc raw"""
```julia
getAw(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
```

Calculate ``i⟨u_n^{(W)}|∂_{k_α}u_m^{(W)}⟩``.
"""
function getAw(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return getdAw(tm, α, (0, 0, 0), k)
end


@memoize k function getdHbar(tm::AbstractTBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dH = getdH(tm, order, k)
    V = geteig(tm, k).vectors
    return V'*dH*V
end


@memoize k function getdSbar(tm::AbstractTBModel, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dS = getdS(tm, order, k)
    V = geteig(tm, k).vectors
    return V'*dS*V
end

@memoize k function getdAwbar(tm::AbstractTBModel, α::Int64, order::Tuple{Int64,Int64,Int64}, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    dAw = getdAw(tm, α, order, k)
    V = geteig(tm, k).vectors
    return V'*dAw*V
end


function getAwbar(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    return getdAwbar(tm, α, (0, 0, 0), k)
end


function getorder(α::Int64)::Tuple{Int64, Int64, Int64}
    order = zeros(Int64, 3)
    order[α] += 1
    return Tuple(order)
end

function getorder(α::Int64, β::Int64)::Tuple{Int64, Int64, Int64}
    order = zeros(Int64, 3)
    order[α] += 1; order[β] += 1
    return Tuple(order)
end


@memoize k function getD(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    order = getorder(α)
    dHbar = getdHbar(tm, order, k); dSbar = getdSbar(tm, order, k)
    Es = geteig(tm, k).values
    D = zeros(ComplexF64, tm.norbits, tm.norbits)
    Awbar = getAwbar(tm, α, k)
    for m in 1:tm.norbits, n in 1:tm.norbits
        En = Es[n]; Em = Es[m]
        if abs(En-Em) > DEGEN_THRESH[1]
            D[n, m] = (dHbar[n, m]-Em*dSbar[n, m])/(Em-En)
        else
            D[n, m] = im*Awbar[n, m]
        end
    end
    return D
end


"""
```julia
getdEs(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})
-->dEs::Vector{Float64}
```

Calculate dE/dk for `tm` at `k` in the `α` direction.

Calculation method is provided in [Wang et al, 2019]. The relevant equation
is Eq. (13). Although in that equation, there is no energy different denominator,
it is still implicitly assumed that the band is nondegenerate. Therefore, dEs[n]
is only correct if n is nondegenerate or completely degenerate.

This function is memoized, which means the arguments and results of the function should
never be modified.
"""
@memoize k function getdEs(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Vector{Float64}
    order = getorder(α)
    dHbar = getdHbar(tm, order, k)
    dSbar = getdSbar(tm, order, k)
    Es = geteig(tm, k).values
    dEs = zeros(tm.norbits)
    for n in 1:tm.norbits
        dEs[n] = real(dHbar[n, n]-Es[n]*dSbar[n, n])
    end
    return dEs
end

@doc raw"""
```julia
getdEs(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    k::AbstractVector{<:Real}
) => dEs::Vector{Float64}
```

Calculate d^2 E / dkα dkβ for `tm` at `k`. `α` and `β` are Cartesian directions.

dEs[n] is only correct if n is nondegenerate or completely degenerate.

This function is memoized, which means the arguments and results of the function should
never be modified.
"""
@memoize k function getdEs(
    tm::AbstractTBModel,
    α::Integer,
    β::Integer,
    k::AbstractVector{<:Real}
)::Vector{Float64}
    Es, _ = geteig(tm, k)
    dHαbar = getdHbar(tm, getorder(α), k)
    dHβbar = getdHbar(tm, getorder(β), k)
    dSαbar = getdSbar(tm, getorder(α), k)
    dSβbar = getdSbar(tm, getorder(β), k)
    dHαβbar = getdHbar(tm, getorder(α, β), k)
    dSαβbar = getdSbar(tm, getorder(α, β), k)
    Dα = getD(tm, α, k)
    Dβ = getD(tm, β, k)
    dEs = zeros(tm.norbits)
    foo1 = dHαbar * Dβ + dHβbar * Dα
    foo2 = dSαbar * Dβ + dSβbar * Dα
    foo3α = getdEs(tm, α, k)
    foo3β = getdEs(tm, β, k)
    for n in 1:tm.norbits
        dEs[n] += real(dHαβbar[n, n] - dSαβbar[n, n] * Es[n])
        dEs[n] += real(foo1[n, n])
        dEs[n] -= real(foo2[n, n] * Es[n])
        dEs[n] -= real(dSαbar[n, n] * foo3β[n] + dSβbar[n, n] * foo3α[n])
        dEs[n] -= real(Dα[n, n] * foo3β[n] + Dβ[n, n] * foo3α[n])
    end
    return dEs
end


"""
```julia
getA(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})
-->A::Matrix{ComplexF64}
```

Calculate Berry connection ``i⟨u_n|∂_α|u_m⟩`` for `tm` at `k`.

Calculation method is provided in [Wang et al, 2019]. The relevant equation
is Eq. (14). Since Eq. (14) assumes band m is nondegenerate, A[n, m] is only
correct if m is nondegenerate or completely degenerate.
"""
@memoize k function getA(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    A = im*getD(tm, α, k)+getAwbar(tm, α, k)
    return A
end


"""
```julia
function getdr(tm::AbstractTBModel, α::Int64, β::Int64, k::AbstractVector{<:Real})
-->dr::Matrix{ComplexF64}
```

Compute ``∂_β r_α`` for `tm` at `k`. r[n, m] = A[n, m] if n != m and r[n, n] = 0.

dr is calculated by directly differentiating Eq. (14) of [Wang et al, 2019].
dr[n, m] is only correct when (i) both band n and band m are nondegenerate or
(ii) both band n and band m are completely degenerate but ``E_n≠E_m``.
"""
@memoize k function getdr(tm::AbstractTBModel, α::Int64, β::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    orderα = getorder(α); orderβ = getorder(β); orderαβ = getorder(α, β)
    dHαbar = getdHbar(tm, orderα, k); dHαβbar = getdHbar(tm, orderαβ, k)
    dSαbar = getdSbar(tm, orderα, k); dSαβbar = getdSbar(tm, orderαβ, k)
    Awαbar = getAwbar(tm, α, k); dAwαβbar = getdAwbar(tm, α, orderβ, k)
    Es = geteig(tm, k).values
    dEs = getdEs(tm, β, k)
    Dα = getD(tm, α, k); Dβ = getD(tm, β, k)
    dr = zeros(ComplexF64, tm.norbits, tm.norbits)
    tmpH = dHαbar*Dβ+dHαβbar+Dβ'*dHαbar; tmpS = dSαbar*Dβ+dSαβbar+Dβ'*dSαbar
    tmpA = Awαbar*Dβ+dAwαβbar+Dβ'*Awαbar
    for m in 1:tm.norbits, n in 1:tm.norbits
        En = Es[n]; Em = Es[m]; dEn = dEs[n]; dEm = dEs[m]
        if abs(En-Em) > DEGEN_THRESH[1]
            dr[n, m] += im*tmpH[n, m]/(Em-En)
            dr[n, m] -= im*dEm*dSαbar[n, m]/(Em-En)
            dr[n, m] -= im*Em*tmpS[n, m]/(Em-En)
            dr[n, m] -= im*(dEm-dEn)*Dα[n, m]/(Em-En)
            dr[n, m] += tmpA[n, m]
        end
    end
    return dr
end


@doc raw"""
```julia
getvelocity(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Number})
```

Calculate velocity matrix in the α direction.

Velocity matrix is calculated by the following expression
```math
v_{nm}^α = ∂_α ϵ_n δ_{nm} + i (ϵ_n-ϵ_m) A_{nm}^α.
```
Therefore, the velocity is actually ħ*velocity.

v[n, m] is only correct when band m is nondegenerate or completely degenerate.
"""
@memoize k function getvelocity(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})::Matrix{ComplexF64}
    v = zeros(ComplexF64, tm.norbits, tm.norbits)
    # diagonal elements
    dEs = getdEs(tm, α, k)
    for n in 1:tm.norbits
        v[n, n] = dEs[n]
    end
    # off-diagonal elements
    Aα = getA(tm, α, k)
    Es = geteig(tm, k).values
    for m in 1:tm.norbits, n in 1:tm.norbits
        if n != m
            v[n, m] = im*(Es[n]-Es[m])*Aα[n, m]
        end
    end
    return v
end


"""
```julia
getspin(tm::AbstractTBModel, α::Int64, k::AbstractVector{<:Real})
```

Calculate ⟨n|σα|m⟩ at `k` point.

If `tm` is a TBModel, the function checks whether `tm.isspinful` is true.
"""
@memoize k function getspin(
    tm::AbstractTBModel,
    α::Integer,
    k::AbstractVector{<:Real}
)::Matrix{ComplexF64}
    if tm isa TBModel && !(tm.isspinful === true)
        error("TBModel should be spinful.")
    end
    length(k) == 3 || error("k should be a 3-element vector.")
    α in 1:3 || error("α should be 1, 2 or 3.")
    nspinless = tm.norbits ÷ 2
    V = geteig(tm, k).vectors
    return V' * kron(σs[α], getS(tm, k)[1:nspinless, 1:nspinless]) * V
end


"""
```julia
get_berry_curvature(tm::AbstractTBModel, α::Int64, β::Int64, k::Vector{<:Real})::Vector{Float64}
```

Calculate Berry curvature Ω for `tm` at `k`.

Ω[n] is only correct if band n is nondegenerate or completely degenerate.
"""
@memoize k function get_berry_curvature(tm::AbstractTBModel, α::Int64, β::Int64, k::Vector{<:Real})
    _, V = geteig(tm, k)
    Sbar_α = V' * getdS(tm, getorder(α), k) * V
    Sbar_β = V' * getdS(tm, getorder(β), k) * V
    Abar_α = V' * HopTB.getAw(tm, α, k) * V
    Abar_β = V' * HopTB.getAw(tm, β, k) * V
    D_α = HopTB.getD(tm, α, k)
    D_β = HopTB.getD(tm, β, k)
    dAw_βα = HopTB.getdAw(tm, β, getorder(α), k)
    dAw_αβ = HopTB.getdAw(tm, α, getorder(β), k)
    Ωbar_αβ = V' * (dAw_βα - dAw_αβ) * V
    return real(diag(Ωbar_αβ - Sbar_α * Abar_β + Sbar_β * Abar_α - im * D_α * D_β +
        im * D_β * D_α - D_α * Abar_β + Abar_β * D_α + D_β * Abar_α - Abar_α * D_β))
end
