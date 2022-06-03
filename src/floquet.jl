module Floquet

using LinearAlgebra
using ..Hop
using HCubature
using SpecialFunctions
using Distributed


"""
```julia
mutable struct FTBModel{T<:Number}
    Ω::Float64
    tm::TBModel{T}
    num_independent_orbits::Int64
    forder::Int64
end
```

`FTBModel` is a struct representing a Floquet tight binding model, which is
is a tight binding model (whose number of orbits being `num_independent_orbits`)
with a periodic Hamiltonian (frequency being `Ω`).
A Floquet tight binding model can be mapped to a
tight binding model with an electric field, stored as `tm`. `tm` should
in principle have infinite number of orbits and we
cut off `tm` from `-forder` to `forder` in the frequency dimension.

The basis function for `tm` is |n⟩e^{-ilΩt}, where `n` labels tight binding
orbitals of the original model. The basis function above corresponds to
the `(l+forder)*num_independent_orbits+n` row/column in the matrix.
"""
mutable struct FTBModel{T<:Number}
    Ω::Float64
    tm::TBModel{T}
    num_independent_orbits::Int64
    forder::Int64
end


"""
```julia
sethopping!(ftm::FTBModel, R::AbstractVector{Int64}, n::Int64, m::Int64,
    p::Int64, hopping::Number)
```

Set ⟨⟨0; n e^{-ilΩt}|H(t)|R; m e^{-i(l-p)Ωt}⟩⟩ to `hopping` for all l.
Alternatively, `hopping` is ⟨n|H(p)|m⟩, where H(t) = ∑H(p)e^{-ipΩt}.

The Floquet Hamiltonian H_F=H-i∂t, in the basis functions |n⟩e^{-ilΩt}, is
(H_F)_{l, l-p}=⟨⟨0; n e^{-ilΩt}|H(t)|R; m e^{-i(l-p)Ωt}⟩⟩-δ_{p,0}lΩ.
This function automatically handles the above δ function.
"""
function Hop.sethopping!(ftm::FTBModel, R::AbstractVector{<:Integer}, n::Int64, m::Int64,
    p::Int64, hopping::Number)
    niorbs = ftm.num_independent_orbits
    δ = p==0 ? 1 : 0
    overl = R in keys(ftm.tm.overlaps) ? ftm.tm.overlaps[R][n, m] : 0
    forder = ftm.forder
    for l in max(-forder+p, -forder):min(forder, forder+p)
        sethopping!(ftm.tm, R, (l+forder)*niorbs+n, (l-p+forder)*niorbs+m,
            hopping-δ*l*ftm.Ω*overl)
    end
end


"""
```julia
addhopping!(ftm::FTBModel, R::AbstractVector{Int64}, n::Int64, m::Int64,
    p::Int64, hopping::Number)
```

Add `hopping` to ⟨⟨0; n e^{-ilΩt}|H(t)|R; m e^{-i(l-p)Ωt}⟩⟩.

This function does not add the onsite energy due to periodic driving.
"""
function Hop.addhopping!(ftm::FTBModel, R::AbstractVector{Int64}, n::Int64, m::Int64,
    p::Int64, hopping::Number)
    niorbs = ftm.num_independent_orbits
    forder = ftm.forder
    for l in max(-forder+p,-forder):min(forder,forder+p)
        addhopping!(ftm.tm, R, (l+forder)*niorbs+n, (l-p+forder)*niorbs+m, hopping)
    end
end


"""
```julia
function getFTBModel(tm::TBModel{T}, Ω::Float64, forder::Int64;
    keep_original_hoppings::Bool=false) where T<:Number
```

Construct a `FTBModel` from a `TBModel`.
"""
function getFTBModel(tm::TBModel{T}, Ω::Float64, forder::Int64;
    keep_original_hoppings::Bool=false) where T<:Number
    ntm = TBModel{T}(tm.norbits*(2forder+1), tm.lat, isorthogonal=tm.isorthogonal)
    if !tm.isorthogonal
        for R in keys(tm.overlaps), l in -forder:forder, m in 1:tm.norbits, n in 1:tm.norbits
            n′ = n+(l+forder)*tm.norbits
            m′ = m+(l+forder)*tm.norbits
            setoverlap!(ntm, R, n′, m′, tm.overlaps[R][n, m])
        end
    end
    for R in keys(tm.positions), α in 1:3, l in -forder:forder, m in 1:tm.norbits, n in 1:tm.norbits
        n′ = n+(l+forder)*tm.norbits
        m′ = m+(l+forder)*tm.norbits
        setposition!(ntm, R, n′, m′, α, tm.positions[R][α][n, m])
    end
    ftm = FTBModel{T}(Ω, ntm, tm.norbits, forder)
    if keep_original_hoppings
        for (R, hopping) in tm.hoppings
            for m in 1:tm.norbits, n in 1:tm.norbits
                sethopping!(ftm, R, n, m, 0, tm.hoppings[R][n, m])
            end
        end
    else
        for n in 1:tm.norbits
            sethopping!(ftm, [0, 0, 0], n, n, 0, 0.0)
        end
    end
    return ftm
end


Hop.geteig(ftm::FTBModel, k::Vector{<:Real}) = geteig(ftm.tm, k)

Hop.getH(ftm::FTBModel, k::Vector{<:Real}) = getH(ftm.tm, k)


"""
Calculate the nth Fourier components for exp(im*x*sin(t))exp(im*y*cos(t)).

exp(im*x*sin(t))exp(im*y*cos(t)) = ∑_p _get_fourier_components(x, y, p)*exp(-ipt)
exp(im*x*sin(t))exp(im*y*cos(t)) = ∑_p ∑_q i^q J_{-p-q}(x) J_q(y) e^{-ipt}
"""
function _get_fourier_components(x::Float64, y::Float64, p::Int64; cutoff::Int64=10)
    return sum([(1.0im)^q*besselj(-p-q, x)*besselj(q, y) for q in -cutoff:cutoff])
end


"""
A(t) = real(A*exp(iΩt))

⟨0n|H'(t)|Rm⟩ = ⟨0n|H|Rm⟩*exp(im*(A(t)⋅(R+rm-rn)))
where the negative charge of electron has been taken into account.
"""
function addlight(tm::TBModel{T}, A::Vector{ComplexF64}, Ω::Float64,
    forder::Int64) where T<:Number
    ftm = getFTBModel(tm, Ω, forder, keep_original_hoppings=false)
    Ar = real.(A)
    Ai = imag.(A)
    for (R, hopping) in tm.hoppings
        if R == Hop.R0 || !(R in keys(ftm.tm.hoppings))
            Rc = tm.lat*R
            for n = 1:tm.norbits
                rn = Hop.get_orbital_position(tm, n)
                for m = 1:tm.norbits
                    rm = Hop.get_orbital_position(tm, m)
                    rnm = rm-rn+Rc
                    x = -Ai⋅rnm
                    y = Ar⋅rnm
                    for p = -ftm.forder:ftm.forder
                        sethopping!(ftm, R, n, m, p, hopping[n, m]*_get_fourier_components(x, y, p))
                    end
                end
            end
        end
    end
    return ftm
end

end
