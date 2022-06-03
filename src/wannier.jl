module Wannier

using LinearAlgebra, HCubature
using ..Hop

export getwf, interpolatewf


@doc raw"""
```julia
inner(vs::Matrix{ComplexF64}, lfs::Dict{Vector{Int64}, Matrix{T}},
    k::Vector{<:Real}) where T<:Number
```

Calculate ⟨ψ_nk|g_m⟩, where g is a localized wave function.
`lfs` is stored in format {R: ⟨Rn|gm⟩}.
"""
function inner(vs::Matrix{ComplexF64}, lfs::Dict{Vector{Int64}, Matrix{T}},
    k::Vector{<:Real}) where T<:Number
    nlfs = 0
    for (R, coeff) in lfs
        nlfs = size(coeff, 2)
        break
    end
    result = zeros(ComplexF64, (size(vs, 2), nlfs))
    for (R, coeff) in lfs
        result += exp(-im*2pi*(k⋅R))*(vs')*coeff
    end
    return result
end


@doc raw"""
```julia
project_orthogonalize(vs::Matrix{ComplexF64},
    lfs::Dict{Vector{Int64}, Matrix{T}}, k::Vector{<:Real};
    tol::Float64=0.01) where T<:Number
```

Get a new Bloch basis by project the Bloch wave functions onto localized
functions.

The returned matrix stores ⟨ψ_m^(W)|φ_n⟩ in column, where φ_n is the new Bloch
basis.
"""
function project_orthogonalize(vs::Matrix{ComplexF64},
    lfs::Dict{Vector{Int64}, Matrix{T}}, k::Vector{<:Real};
    tol::Float64=0.1) where T<:Number
    A = inner(vs, lfs, k)
    S = A'*A
    @assert real(det(S)) > tol "det(S) < $tol at $k"
    return vs*A*inv(sqrt(S))
end


"""
```
getwf(getvs::Function, Rs::AbstractMatrix{Int64}; rtol::Float64=√eps(),
    atol::Float64=0.0, constant_components::Vector{Int64}=zeros(Int64, 0),
    ndiv::Int64=1)
```

Calculate Wannier functions in the home unit cell.

`getvs` should be a function that takes a three-component vector k as argument and
yields a matrix (orbital_index, band_index). Wannier functions will only be computed
at `Rs`. The other parameters are directly sent to `Hop.Wannier.pintegrate`.
"""
function getwf(getvs::Function, Rs::AbstractMatrix{Int64}; rtol::Float64=√eps(),
    atol::Float64=0.0, constant_components::Vector{Int64}=zeros(Int64, 0),
    ndiv::Int64=1)

    norbits, nwfs = size(getvs(zeros(3)))
    nRs = size(Rs, 2)

    function getitrd(k)
        vs = getvs(k)
        tmp = zeros(ComplexF64, norbits, nwfs, nRs)
        for iR in 1:nRs
            tmp[:, :, iR] = exp(im*2π*(k⋅Rs[:, iR]))*vs
        end
        return tmp
    end

    wfstmp, _ = Hop.Utilities.pintegrate(getitrd, zeros(3), ones(3),
        rtol=rtol, atol=atol, constant_components=constant_components, ndiv=ndiv)

    wfs = Dict{Vector{Int64}, Matrix{ComplexF64}}()

    for iR in 1:nRs
        wfs[Rs[:, iR]] = wfstmp[:, :, iR]
    end

    return wfs
end

"""
```
getwf(tm::TBModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, Rs::AbstractMatrix{Int64}; tol::Float64=0.1, rtol::Float64=√eps(),
    atol::Float64=0.0, constant_components::Vector{Int64}=zeros(Int64, 0), ndiv::Int64=1) where T<:Number
```

Calculate Wannier functions with trial Wannier functions in `twfs`.

The project orthogonalization procedure has a tolerance `tol`. Other keyword
arguments are directly sent to `Hop.Wannier.pintegrate`.
"""
function getwf(tm::TBModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, Rs::AbstractMatrix{Int64}; tol::Float64=0.1, rtol::Float64=√eps(),
    atol::Float64=0.0, constant_components::Vector{Int64}=zeros(Int64, 0), ndiv::Int64=1) where T<:Number
    for i in bands
        i in 1:tm.norbits || error("bands not in the range of 1..tm.norbits.")
    end
    function getvs(k)
        egvecs = geteig(tm, k).vectors[:, bands]
        return project_orthogonalize(egvecs, twfs, k, tol=tol)
    end
    return getwf(getvs, Rs, rtol=rtol, atol=atol, constant_components=constant_components, ndiv=ndiv)
end

"""
```julia
getwf(getvs::Function, nkmesh::Vector{Int64}, nrmesh::Vector{Int64})
```

Calculate home unit cell Wannier functions with k mesh `nkmesh`.

Wannier functions are only computed at unit cells denoted by `nrmesh`.
"""
function getwf(getvs::Function, nkmesh::Vector{Int64}, nrmesh::Vector{Int64})
    @assert size(nkmesh, 1) == 3
    @assert size(nrmesh, 1) == 3
    nbasis, nwfs = size(getvs([0.0, 0.0, 0.0]))
    # generate k points
    nkpts = prod(nkmesh)
    kpts = zeros(3, nkpts)
    i = 1
    for ikx in 1:nkmesh[1], iky in 1:nkmesh[2], ikz in 1:nkmesh[3]
        kpts[:, i] = [(ikx-1)/nkmesh[1], (iky-1)/nkmesh[2], (ikz-1)/nkmesh[3]]
        i += 1
    end
    # generate r points
    nrpts = prod(2*nrmesh.+1)
    wf = Dict{Vector{Int64}, Matrix{ComplexF64}}()
    for Rx in (-nrmesh[1]):nrmesh[1], Ry in (-nrmesh[2]):nrmesh[2], Rz in (-nrmesh[3]):nrmesh[3]
        wf[[Rx, Ry, Rz]] = zeros(ComplexF64, (nbasis, nwfs))
    end
    # perform integration
    for ik in 1:nkpts
        vs = getvs(kpts[:, ik])
        for R in keys(wf)
            wf[R] += exp(im*2pi*(kpts[:, ik]⋅R))*vs/nkpts
        end
    end

    return wf
end


@doc raw"""
```julia
getwf(tm::TBModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, nkmesh::Vector{Int64},
    nrmesh::Vector{Int64}; tol::Float64=0.01) where T<:Number
```

Calculate Wannier functions for bands indexed by `bands`.

Both returned Wannier functions and `lfs` are stored in format {R: ⟨Rn|gm⟩}, where g
are localized (Wannier) functions.
"""
function getwf(tm::TBModel, twfs::Dict{Vector{Int64}, Matrix{T}},
    bands::Vector{Int64}, nkmesh::Vector{Int64},
    nrmesh::Vector{Int64}; tol::Float64=0.1) where T<:Number
    for i in bands
        @assert i in 1:tm.norbits
    end
    function getvs(k)
        egvecs = geteig(tm, k).vectors[:, bands]
        return project_orthogonalize(egvecs, twfs, k, tol=tol)
    end
    return getwf(getvs, nkmesh, nrmesh)
end


@doc raw"""
```julia
interpolatewf(gettm, tmparams::Matrix{<:Number},
    twfs::Dict{Vector{Int64}, Matrix{T}}, bands::Vector{Int64},
    nkmesh::Vector{Int64}, nrmesh::Vector{Int64}, ndiv::Int64) where T<:Number
```

Calculate Wannier functions by interpolation between tight binding models.

`gettm` is a function taking parameters and returning tight binding models.
With the parameters stored in column in tmparams, the Wannier function is obtained
by interpolation. More precisely, the `twfs` are used for initial guess and Wannier
functions from last step is used as trial Wannier function for the next step.
"""
function interpolatewf(gettm, tmparams::Matrix{<:Number},
    twfs::Dict{Vector{Int64}, Matrix{T}}, bands::Vector{Int64},
    nkmesh::Vector{Int64}, nrmesh::Vector{Int64}, ndiv::Int64) where T<:Number
    tm = gettm(tmparams[:, 1]...)
    wf = getwf(tm, twfs, bands, nkmesh, nrmesh)
    for ipath in 1:(size(tmparams, 2)-1)
        for itm in 1:ndiv
            tmparam = tmparams[:, ipath]+(tmparams[:, ipath+1]-tmparams[:, ipath])*itm/ndiv
            tm = gettm(tmparam...)
            wf = getwf(tm, wf, bands, nkmesh, nrmesh)
        end
    end
    return wf
end

end
