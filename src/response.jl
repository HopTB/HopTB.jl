module Response

"""
THIS MODULE ONLY WORKS FOR MODELS WITH ORTHOGONAL BASIS!
"""

using LinearAlgebra, Distributed
using ..HopTB.Utilities: integrate, constructmeshkpts, splitkpts
using ..HopTB.Parallel: ParallelFunction, claim!, stop!

"""
```julia
getresponse(getH::Function, getO::Function, getV::Function,
    getn::Function, k::AbstractVector{Float64}; Γ::Float64=0.1, atol::Float64=0.0,
    maxevals::Int64=typemax{Int64})
```

Calculate response at `k`.

The response is the linear change of ⟨O⟩ due to V added to the original Hamiltonian.
`getH`, `getO` and `getV` all take `k` as the only argument and
output corresponding matrices.

`Γ` is the broadening factor. If the every degree of the system is
coupled to a wide-band heat bath, Γ=i T^2 ν, where T is the hopping constant
between the system and the heat bath and ν is the density of states of the
heat bath. The heat bath has a occupation distribution `getn(ω)`.

The integration over ω would have an absolute tolerance `atol`, with maximal
evaluation number being `maxevals`.
"""
function getresponse(getH::Function, getO::Function, getV::Function,
    getn::Function, k::AbstractVector{Float64}; Γ::Float64=0.1, atol::Float64=0.0,
    maxevals::Int64=typemax(Int64))
    O = getO(k)
    V = getV(k)
    egvals, egvecs = eigen(Hermitian(getH(k)))
    f = let egvals=egvals, egvecs=egvecs, k=k, Γ=Γ, getn=getn, O=O, V=V
            function f(ω)
                GR = egvecs*Diagonal(1 ./(ω.-egvals.+im*Γ))*(egvecs')
                GA = GR'
                Σ = 2*im*Γ*getn(ω) # this is lesser self energy
                tmp = GR*V*GR*Σ*GA
                -im*(tr(O*(tmp-tmp')))/(2π)
            end
        end
    return integrate(f, atol=atol, maxevals=maxevals)[1]
end


"""
```julia
getresponse(getH::Function, getO::Function, getV::Function,
    getn::Function, nkmesh::Vector{Int64}; Γ::Float64=0.1, atol::Float64=0.0,
    maxevals::Int64=typemax(Int64), offset::Vector{Float64}=[0.0, 0.0, 0.0],
    k1::Vector{Float64}=[0.0, 0.0, 0.0], k2::Vector{Float64}=[1.0, 1.0, 1.0])
```

Calculate response.

The response is the linear change of ⟨O⟩ due to V added to the original Hamiltonian.
`getH`, `getO` and `getV` all take `k` as the only argument and
output corresponding matrices.

`Γ` is the broadening factor. If the every degree of the system is
coupled to a wide-band heat bath, Γ=i T^2 ν, where T is the hopping constant
between the system and the heat bath and ν is the density of states of the
heat bath. The heat bath has a occupation distribution `getn(ω)`.

The response function is integrated in a rectangle with corners `k1` and `k2`.
The mesh is `nkmesh` and every k point is offsetted by `offset`.

The integration over ω would have an absolute tolerance `atol`, with maximal
evaluation number being `maxevals` for each k point.

For a crystal, since this function does not know the volume of the Brillouin
zone, it simply sums over k points and divide the result by the number of k
points. The most common situation is calculating the density of ⟨O⟩, in which case
the result should be multiplied by the volume of the Brillouin zone and divided
by (2π)^3.
"""
function getresponse(getH::Function, getO::Function, getV::Function,
    getn::Function, nkmesh::Vector{Int64}; Γ::Float64=0.1, atol::Float64=0.0,
    maxevals::Int64=typemax(Int64), offset::Vector{Float64}=[0.0, 0.0, 0.0],
    k1::Vector{Float64}=[0.0, 0.0, 0.0], k2::Vector{Float64}=[1.0, 1.0, 1.0])

    length(nkmesh) == 3 || error("nkmesh should be a 3-element vector.")
    nks = prod(nkmesh)
    ks = constructmeshkpts(nkmesh; offset=offset, k1=k1, k2=k2)
    result = 0.0
    for i in 1:nks
        k = ks[:, i]
        result += getresponse(getH, getO, getV, getn, k;
            Γ=Γ, atol=atol, maxevals=maxevals)
    end
    return result/nks
end


"""
```julia
get_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, k::Vector{Float64}; Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))
```

Calculate Floquet response at `k`.

The response is the linear change of ⟨O⟩ due to V added to the original Hamiltonian.
`getH`, `getO` and `getV` all take `k` as the only argument and
output corresponding matrices.

We here assume operator O and V are time independent, such that their
size should be number of orbital. However, H should be Floquet Hamiltonian
and should have size (2*forder+1)*norbits.

⟨O⟩ might be time dependent but periodic. `l` denotes that the output is
the coefficient of ⟨O⟩ before e^{-i l Ω t}.

`Γ` is the broadening factor. If the every degree of the system is
coupled to a wide-band heat bath, Γ=i T^2 ν, where T is the hopping constant
between the system and the heat bath and ν is the density of states of the
heat bath. The heat bath has a occupation distribution `getn(ω)`.

The integration over ω would have an absolute tolerance `atol`, with maximal
evaluation number being `maxevals`.
"""
function get_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, k::Vector{Float64}; Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))

    norbits = size(getO([0.0, 0.0, 0.0]), 1)
    tmp = getH([0.0, 0.0, 0.0])
    forder = (size(tmp, 1)÷norbits-1)÷2
    Ω = real(tmp[1, 1]-tmp[norbits+1, norbits+1])
    O = getO(k)
    V = kron(I(2*forder+1), getV(k))
    index1 = (l+forder)*norbits+1:(l+forder+1)*norbits
    index2 = forder*norbits+1:(forder+1)*norbits
    egvals, egvecs = eigen(Hermitian(getH(k)))
    Vbar = egvecs'*V*egvecs
    Obar = egvecs'[:, index2]*O*egvecs[index1, :]
    f = let egvals=egvals, egvecs=egvecs, Γ=Γ, getn=getn, forder=forder,
            norbits=norbits, Ω=Ω, Obar=Obar, Vbar=Vbar
            function f(ω)
                GR = Diagonal(1 ./(ω.-egvals.+im*Γ))
                GA = GR'
                # This is Floquet lesser self energy
                Σ = 2*im*Γ*kron(Diagonal([getn(ω+n*Ω) for n in -forder:forder]),
                    I(norbits))
                Σbar = egvecs'*Σ*egvecs
                tmp = GR*Vbar*GR*Σbar*GA
                -im*tr(Obar*(tmp-tmp'))/(2π)
            end
        end
    return integrate(f, atol=atol, maxevals=maxevals)[1]
end


function Floquet_response_worker(ks::Matrix{Float64}, getH::Function,
    getO::Function, getV::Function, getn::Function, l::Int64, Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))

    nks = size(ks, 2)
    result = 0.0im
    for i in 1:nks
        k = ks[:, i]
        result += get_Floquet_response(getH, getO, getV, getn, l, k;
            Γ=Γ, atol=atol, maxevals=maxevals)
    end
    return result
end


"""
```julia
get_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, nkmesh::Vector{Int64}; Γ::Float64=0.1,
    offset::Vector{Float64}=[0.0, 0.0, 0.0], k1::Vector{Float64}=[0.0, 0.0, 0.0],
    k2::Vector{Float64}=[1.0, 1.0, 1.0], atol::Float64=0.0, maxevals::Int64=typemax(Int64))
```

Calculate Floquet response.

The response is the linear change of ⟨O⟩ due to V added to the original Hamiltonian.
`getH`, `getO` and `getV` all take `k` as the only argument and
output corresponding matrices.

We here assume operator O and V are time independent, such that their
size should be number of orbital. However, H should be Floquet Hamiltonian
and should have size (2*forder+1)*norbits.

⟨O⟩ might be time dependent but periodic. `l` denotes that the output is
the coefficient of ⟨O⟩ before e^{-i l Ω t}.

`Γ` is the broadening factor. If the every degree of the system is
coupled to a wide-band heat bath, Γ=i T^2 ν, where T is the hopping constant
between the system and the heat bath and ν is the density of states of the
heat bath. The heat bath has a occupation distribution `getn(ω)`.

The response function is integrated in a rectangle with corners `k1` and `k2`.
The mesh is `nkmesh` and every k point is offsetted by `offset`.

The integration over ω would have an absolute tolerance `atol`, with maximal
evaluation number being `maxevals` for each k point.

For a crystal, since this function does not know the volume of the Brillouin
zone, it simply sums over k points and divide the result by the number of k
points. The most common situation is calculating the density of ⟨O⟩, in which case
the result should be multiplied by the volume of the Brillouin zone and divided
by (2π)^3.
"""
function get_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, nkmesh::Vector{Int64}; Γ::Float64=0.1,
    offset::Vector{Float64}=[0.0, 0.0, 0.0], k1::Vector{Float64}=[0.0, 0.0, 0.0],
    k2::Vector{Float64}=[1.0, 1.0, 1.0], atol::Float64=0.0, maxevals::Int64=typemax(Int64))

    length(nkmesh) == 3 || error("nkmesh should be a 3-element vector.")
    nks = prod(nkmesh)
    ks = constructmeshkpts(nkmesh; offset=offset, k1=k1, k2=k2)
    klist = splitkpts(ks, 2*nworkers())
    pf = ParallelFunction(Floquet_response_worker, getH, getO, getV, getn, l, Γ,
        atol, maxevals, len=2*nworkers())
    for i in 1:2*nworkers()
        pf(klist[i])
    end
    result = 0.0im
    for i in 1:2*nworkers()
        result += claim!(pf)
    end
    stop!(pf)
    return result/nks
end


function clt_Floquet_response_worker(ks::Matrix{Float64}, getH::Function,
    getO::Function, getV::Function, getn::Function, l::Int64, Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))

    nks = size(ks, 2)
    result = zeros(ComplexF64, nks)
    for i in 1:nks
        k = ks[:, i]
        result[i] = get_Floquet_response(getH, getO, getV, getn, l, k;
            Γ=Γ, atol=atol, maxevals=maxevals)
    end
    return (ks, result)
end


"""
```julia
clt_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, ks::Matrix{Float64}; Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))
```

Collect k-resolved Floquet response in `ks`. `ks` is in reduced coordinates and
stored in columns.

Most options are identical to `get_Floquet_response`.
"""
function clt_Floquet_response(getH::Function, getO::Function, getV::Function,
    getn::Function, l::Int64, ks::Matrix{Float64}; Γ::Float64=0.1,
    atol::Float64=0.0, maxevals::Int64=typemax(Int64))

    klist = splitkpts(ks, 2*nworkers())
    pf = ParallelFunction(clt_Floquet_response_worker, getH, getO, getV, getn, l, Γ,
        atol, maxevals, len=2*nworkers())
    for i in 1:2*nworkers()
        pf(klist[i])
    end
    ks = zeros((3, 0))
    χs = zeros(ComplexF64, 0)
    for i in 1:2*nworkers()
        tmp = claim!(pf)
        ks = [ks tmp[1]]
        χs = [χs; tmp[2]]
    end
    stop!(pf)
    return ks, χs
end


end
