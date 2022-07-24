using SharedArrays, StaticArrays

export AbstractTBModel, TBModel, SharedTBModel
export sethopping!, addhopping!, setoverlap!, setposition!, set_orbital_types!, set_is_canonical_ordered!
export gethopping, changebasis


const R0 = [0, 0, 0]
const IMAG_TOL = 1.0e-4
const HERM_TOL = 1.0e-4


@doc raw"""
Abstract type for tight binding models.

AbstractTBModel should at least implement 3 functions in basic.jl: `getdH`, `getdS`, `getdAw`.
"""
abstract type AbstractTBModel{T<:Number} end


@doc raw"""
`TBModel{T<Number}` is a data type representing tight binding models.

# Fields
 - `norbits::Int64`: number of orbits in the model
 - `lat::SMatrix{3,3,Float64,9}`: primitive lattice vectors stored in columns
 - `rlat::SMatrix{3,3,Float64,9}`: primitive reciprocal lattice vectors stored in columns
 - `hoppings::Dict{SVector{3,Int16},Matrix{T}}`: R -> ⟨0n|H|Rm⟩, where n is the
   first index and m is the second index.
 - `overlaps::Dict{SVector{3,Int16},Matrix{T}}`: R -> ⟨0n|Rm⟩
 - `positions::Dict{SVector{3,Int16},SVector{3,Matrix{T}}}`: R -> [⟨0n|rx|Rm⟩, ⟨0n|ry|Rm⟩, ⟨0n|rz|Rm⟩]
 - `isorthogonal`: whether the different orbitals of the model is orthogonal to each
   other.
 - `nsites::Union{Missing,Int16}`: It is possible that the orbitals can be grouped into several sites. Every group
   of orbitals should share the same position within the group, stored in `site_positions`.
 - `site_norbits::Union{Missing,Vector{Int64}}`: number of orbitals for each site.
 - `site_positions::Union{Missing,Matrix{Float64}}`: position of sites stored in column
 - `orbital_types::Union{Missing,Vector{Vector{Int16}}}`: It is possible that
   the orbitals have definite symmetry representations. Example: `orbital_types = [[0], [0, 1]]`
   denotes that there are two sites, the first of which has 1 s orbital and the second
   of which has 1 s orbital and 1 p orbital (4 orbitals in total not counting spin).
   The orbitals of the TBModel should appear in a consistent order denoted by `orbital_types`.
 - `isspinful::Union{Missing,Bool}`: If a TBModel is spinful, the first half of the orbitals
   should be spin up and the second half should be spin down. The two halves of the orbitals
   should be in one-to-one correspondence.
 - `is_canonical_ordered::Union{Missing,Bool}`: If the orbital types are know, and
   `is_canonical_ordered` is true, then the orbitals are guaranteed to appear
   in a canonical order defined by decreasing Lz value. For example, canonical
   order of p orbitals should -px-i*py, pz, px-i*py. See the wikipedia for
   other orbitals.

# Missing data

Not all fields of TBModels are required to have a valid value. For example,
it is possible the site information is missing for some models. Here are some
rules:
 - `nsites`, `site_norbits`, `site_positions` should all either be missing or valid.
 - `orbital_types` can only be valid if `nsites`, `site_norbits` and `site_positions`
   are valid.
 - `is_canonical_ordered` can only be valid if `orbital_types` is valid.

# Consistency check

Any functions that directly modifies a TBModel should always maintain the following
consistencies (if relevant fields are not missing):
 - `hoppings`, `overlaps` and `positions` matrices should always be Hermitian.
 - `overlap[[0, 0, 0]]` should always be positive definite.
 - `diag(position[[0, 0, 0]][α])/diag(overlap[[0, 0, 0]])` should be consistent with
   `site_positions`.
 - number of orbits should be consistent.
"""
mutable struct TBModel{T} <: AbstractTBModel{T}
    norbits::Int64
    lat::SMatrix{3,3,Float64,9}
    rlat::SMatrix{3,3,Float64,9}
    hoppings::Dict{SVector{3,Int16},Matrix{T}}
    overlaps::Dict{SVector{3,Int16},Matrix{T}}
    positions::Dict{SVector{3,Int16},SVector{3,Matrix{T}}}
    isorthogonal::Bool
    nsites::Union{Missing,Int64}
    site_norbits::Union{Missing,Vector{Int16}}
    site_positions::Union{Missing,Matrix{Float64}}
    orbital_types::Union{Missing,Vector{Vector{Int16}}}
    isspinful::Union{Missing,Bool}
    is_canonical_ordered::Union{Missing,Bool}
end


@doc raw"""
`SharedTBModel` is a data type that encodes all the information of a `TBModel`
into several SharedArrays.

`SharedTBModel` is generally more efficient than `TBModel`. A typical workflow
is to construct `TBModel` and then convert it into a `SharedTBModel` for
large scale calculations.

# Fields

 - `norbits::Int64`: number of orbits
 - `isorthogonal::Bool`: whether orbits of the model are orthonormal
 - `lat::SMatrix{3,3,Float64,9}`: primitive lattice vectors stored in columns
 - `rlat::SMatrix{3,3,Float64,9}`: primitive reciprocal lattice vectors stored in columns
 - `Rs::Matrix{Int16}`: R vectors stored in columns for hopping, overlap and position matrices.
   The number of R vectors is odd. The first R vector is [0, 0, 0]. For the rest of the matrix
   `Rs[:, 2:end]`, the first half and the second half is guaranteed to be in one-to-one
   correspondence R <-> -R.
 - `Rcs::Matrix{Float64}`: R vectors in Cartesian coordinates.
 - `H::SharedArray{T,2}`: `H = reshape(Hmatrix, (norbits * norbits, :))`.
   `Hmatrix[:, :, iR]` is ⟨0n|H|Rm⟩ where R is `Rs[:, iR]`. Since Hamiltonian is
   Hermitian, only first half R vectors in `Rs` are stored. In addition, `Hmatrix[:, :, 1]` is different:
   `Hmatrix[:, :, 1]` is ⟨0n|H|0m⟩/2.
 - `S::Union{Nothing,SharedArray{T,3}}`: `S = reshape(Smatrix, (norbits * norbits, :))`.
   `Smatrix[:, :, iR]` is ⟨0n|Rm⟩ where R is `Rs[:, iR]`. Since overlap matrix is
   Hermitian, only first half R vectors in `Rs` are stored. In addition, `Smatrix[:, :, 1]` is different:
   `Smatrix[:, :, 1]` is ⟨0n|0m⟩/2. 
 - `r:SVector{3,SharedArray{T,2}}`: `S[α] = reshape(Smatrices[α], (norbits * norbits, :))`.
   `Smatrices[α][:, :, iR]` is ⟨0n|rα|Rm⟩ where R is `Rs[:, iR]`.
"""
struct SharedTBModel{T} <: AbstractTBModel{T}
    norbits::Int64
    isorthogonal::Bool
    lat::SMatrix{3,3,Float64,9}
    rlat::SMatrix{3,3,Float64,9}
    Rs::Matrix{Int16}
    Rcs::Matrix{Float64}
    H::SharedArray{T,2}
    S::Union{Nothing,SharedArray{T,2}}
    r::SVector{3,SharedArray{T,2}}
end


function _get_halfRs(tm::TBModel{T}) where T <: Number
    Rs = unique(vcat(collect.(keys.([tm.hoppings, tm.overlaps, tm.positions]))...))
    nRs = length(Rs)
    @assert isodd(nRs)

    nhalfRs = (nRs + 1) ÷ 2
    halfRs = Vector{SVector{3,Int16}}()
    push!(halfRs, R0)
    for R in Rs
        if !(-R in halfRs)
            push!(halfRs, R)
        end
    end
    @assert length(halfRs) == nhalfRs
    return halfRs
end

@doc raw"""
```julia
SharedTBModel(tm::TBModel{T}) where T <: Number
```
"""
function SharedTBModel(tm::TBModel{T}) where T <: Number
    halfRs = _get_halfRs(tm)
    nhalfRs = length(halfRs)
    nRs = 2nhalfRs - 1
    norbits = tm.norbits

    Rs = zeros(Int16, 3, nRs)
    Hmatrix = zeros(T, norbits, norbits, nhalfRs)
    Smatrix = tm.isorthogonal ? nothing : zeros(T, norbits, norbits, nhalfRs)
    rmatrices = [zeros(T, norbits, norbits, nRs) for _ in 1:3]

    kH, kS, kr = keys.([tm.hoppings, tm.overlaps, tm.positions])
    # R ≠ [0, 0, 0]
    for (iR, R) in enumerate(halfRs[2:end])
       iR′ = iR + 1
       Rs[:, iR′] = R
       Rs[:, nhalfRs + iR] = -R
       if R in kH
           Hmatrix[:, :, iR′] = tm.hoppings[R]
       end
       if !tm.isorthogonal && R in kS
           Smatrix[:, :, iR′] = tm.overlaps[R]
       end
       if R in kr
           for i = 1:3
               rmatrices[i][:, :, iR′] = tm.positions[R][i]
               rmatrices[i][:, :, nhalfRs + iR] = tm.positions[-R][i]
           end
       end
    end

    # R = [0, 0, 0]
    Rs[:, 1] = R0
    Hmatrix[:, :, 1] = tm.hoppings[R0] / 2
    if !tm.isorthogonal
        Smatrix[:, :, 1] = tm.overlaps[R0] / 2
    end
    for i = 1:3
        rmatrices[i][:, :, 1] = tm.positions[R0][i]
    end

    return SharedTBModel{T}(
        norbits,
        tm.isorthogonal,
        tm.lat,
        tm.rlat,
        Rs,
        tm.lat * Rs,
        reshape(Hmatrix, (norbits * norbits, :)),
        tm.isorthogonal ? nothing : reshape(Smatrix, (norbits * norbits, :)),
        [reshape(rmatrix, (norbits * norbits, :)) for rmatrix in rmatrices]
    )
end


function Base.show(io::IO, tm::TBModel)
    print(io, "TBModel: $(tm.norbits) orbitals")
end


function Base.show(io::IO, sm::SharedTBModel)
    print(io, "SharedTBModel: $(sm.norbits) orbitals")
end


@doc raw"""
```julia
TBModel{T=ComplexF64}(norbits::Int64, lat::AbstractMatrix{Float64};
    isorthogonal::Bool=true) where T <: Number
```

Construct a TBModel.

Overlap matrix is automatically set as
identity for R = [0, 0, 0] no matter the model is orthogonal or not.
All extra information is missing. 
"""
function TBModel{T}(norbits::Int64, lat::AbstractMatrix{Float64};
    isorthogonal::Bool=true) where T <: Number
    size(lat) == (3, 3) || error("Size of lat is not correct.")
    rlat = 2π * inv(lat)'
    tm = TBModel{T}(norbits, lat, rlat, Dict(), Dict(), Dict(), isorthogonal,
        missing, missing, missing, missing, missing, missing)
    tm.overlaps[R0] = Matrix{T}(I, norbits, norbits)
    return tm
end


@doc raw"""
```julia
TBModel{T=ComplexF64}(lat::AbstractMatrix{Float64}, orbital_positions::AbstractMatrix{Float64};
    isorthogonal::Bool = true) where T <: Number
```

Construct a TBModel.

All extra information is missing. Overlap matrix is automatically set as
identity for R = [0, 0, 0] no matter the model is orthogonal or not. Diagonal
elements of the position matrix for R = [0, 0, 0] is set to the `orbital_positions`.
"""
function TBModel{T}(lat::AbstractMatrix{Float64}, orbital_positions::AbstractMatrix{Float64};
    isorthogonal::Bool=true) where T <: Number

    size(lat) == (3, 3) || error("Size of lat is not correct.")
    size(orbital_positions, 1) == 3 || error("orbital_positions should have three rows.")

    norbits = size(orbital_positions, 2)
    rlat = 2π * inv(lat)'
    tm = TBModel{T}(norbits, lat, rlat, Dict(), Dict(), Dict(), isorthogonal,
        missing, missing, missing, missing, missing, missing)
    tm.overlaps[R0] = I(norbits)
    for n in 1:norbits, α in 1:3
        setposition!(tm, R0, n, n, α, orbital_positions[α, n])
    end

    return tm
end


@doc raw"""
```julia
TBModel{T=ComplexF64}(lat::AbstractMatrix{Float64}, site_positions::Matrix{Float64},
    orbital_types::Vector{Vector{Int64}}; isspinful::Bool=false, isorthogonal::Bool=true,
    is_canonical_ordered::Bool=false) where T <: Number
```

Construct a TBModel.

Overlap matrix is automatically set as identity for R = [0, 0, 0] no matter the model
is orthogonal or not. Diagonal elements of the position matrix for R = [0, 0, 0] is
set according to `site_positions`.
"""
function TBModel{T}(lat::AbstractMatrix{Float64}, site_positions::Matrix{Float64},
    orbital_types::Vector{Vector{Int64}}; isspinful::Bool=false, isorthogonal::Bool=true,
    is_canonical_ordered::Bool=false) where T <: Number
    size(lat) == (3, 3) || error("Size of lat is not correct.")

    rlat = 2π * inv(lat)'
    nsites = size(site_positions, 2)
    nspins = 1 + isspinful
    site_norbits = [sum(2 * orbital_types[i] .+ 1) for i in 1:nsites] * nspins
    norbits = sum(site_norbits)
    tm = TBModel{T}(norbits, lat, rlat, Dict(), Dict(), Dict(), isorthogonal, nsites,
        site_norbits, site_positions, orbital_types, isspinful, is_canonical_ordered)
    tm.overlaps[R0] = I(norbits)
    tm.positions[R0] = [zeros(T, tm.norbits, tm.norbits) for α in 1:3]
    for i in 1:nsites, p in 1:site_norbits[i], α in 1:3
        n = _to_orbital_index(tm, (i, p))
        tm.positions[R0][α][n, n] = site_positions[α, i]
    end
    return tm
end


TBModel(norbits::Int64, lat::AbstractMatrix{Float64}; isorthogonal::Bool=true) =
    TBModel{ComplexF64}(norbits, lat, isorthogonal=isorthogonal)


TBModel(lat::AbstractMatrix{Float64}, orbital_positions::Matrix{Float64}; isorthogonal::Bool=true) =
    TBModel{ComplexF64}(lat, orbital_positions, isorthogonal=isorthogonal)


TBModel(lat::AbstractMatrix{Float64}, site_positions::Matrix{Float64},
    orbital_types::Vector{Vector{Int64}}; isspinful::Bool=false,
    isorthogonal::Bool=true, is_canonical_ordered::Bool=false) =
    TBModel{ComplexF64}(lat, site_positions,
    orbital_types; isspinful=isspinful, isorthogonal=isorthogonal,
    is_canonical_ordered=is_canonical_ordered)


function Base.convert(::Type{TBModel{T}}, tm::TBModel) where T <: Number
    return TBModel{T}(
        tm.norbits,
        tm.lat,
        tm.rlat,
        convert(Dict{SVector{3,Int64},Matrix{T}}, tm.hoppings),
        convert(Dict{SVector{3,Int64},Matrix{T}}, tm.overlaps),
        convert(Dict{SVector{3,Int64},SVector{3,Matrix{T}}}, tm.positions),
        tm.isorthogonal,
        tm.nsites,
        tm.site_norbits,
        tm.site_positions,
        tm.orbital_types,
        tm.isspinful,
        tm.is_canonical_ordered
    )
end


function set_orbital_types!(
    tm::TBModel,
    orbital_types::Vector{Vector{Int64}};
    isspinful::Bool=false, 
    is_canonical_ordered::Bool=false,
    position_tolerance::Real=1.0e-4
)
    nsites = length(orbital_types)
    site_norbits = [sum(2 * orbital_types[i] .+ 1) for i in 1:nsites] * (1 + isspinful)
    sum(site_norbits) == tm.norbits || error("orbital_types is not compatible with this tm: wrong norbits.")
    tm.nsites = nsites
    tm.site_norbits = site_norbits
    tm.isspinful = isspinful
    tm.orbital_types = orbital_types
    site_positions = zeros(3, nsites)
    tm.is_canonical_ordered = is_canonical_ordered
    for α in 1:3, i in 1:nsites
        n = _to_orbital_index(tm, (i, 1))
        site_positions[α, i] = real(tm.positions[R0][α][n, n] / tm.overlaps[R0][n, n])
        for p in 1:site_norbits[i]
            n = _to_orbital_index(tm, (i, p))
            if !isapprox(site_positions[α, i], tm.positions[R0][α][n, n] / tm.overlaps[R0][n, n], atol=position_tolerance)
                remove_extra_information!(tm)
                error("orbital_types is not compatible with this tm: wrong positions for site $(i).")
            end
        end
    end
    tm.site_positions = site_positions
    return nothing
end


function set_is_canonical_ordered!(tm::TBModel, v::Bool)
    tm.is_canonical_ordered = v
    return nothing
end


function has_full_information(tm)
    for field in [:nsites, :site_norbits, :site_positions, :orbital_types,
        :isspinful, :is_canonical_ordered]
        if ismissing(getfield(tm, field))
            return false
        end
    end
    return true
end


function remove_extra_information!(tm)
    for field in [:nsites, :site_norbits, :site_positions, :orbital_types,
        :isspinful, :is_canonical_ordered]
        setfield!(tm, field, missing)
    end
end


@doc raw"""
```julia
sethopping!(tm::TBModel{T}, R::AbstractVector{<:Integer}, n::Int64, m::Int64,
    hopping::Number) where T
```

Set ⟨0n|H|Rm⟩ to `hopping` for `tm`.
"""
function sethopping!(tm::TBModel{T}, R::AbstractVector{<:Integer}, n::Int64, m::Int64, hopping::Number) where T <: Number
    ((n in 1:tm.norbits) && (m in 1:tm.norbits)) || error("n or m not in range 1-norbits.")
    length(R) == 3 || error("R should be a 3-element vector")
    R == R0 && n == m && imag(hopping) > IMAG_TOL && error("On site energy should be real.")


    if !(R in keys(tm.hoppings))
        tm.hoppings[R] = zeros(T, tm.norbits, tm.norbits)
        tm.hoppings[-R] = zeros(T, tm.norbits, tm.norbits)
    end

    if R == R0 && n == m
        tm.hoppings[R][n, m] = real(hopping)
    else
        tm.hoppings[R][n, m] = hopping
        tm.hoppings[-R][m, n] = conj(hopping)
    end
end


function sethopping!(tm::TBModel, R::AbstractVector{<:Integer}, (i, p)::Tuple{Int64,Int64},
    (j, q)::Tuple{Int64,Int64}, hopping::Number)
    has_full_information(tm) || error("No site information is provided in the model.")
    ((i in 1:tm.nsites) && (j in 1:tm.nsites)) || error("i or j not in range 1-nsites.")
    ((p in 1:tm.site_norbits[i]) && (q in 1:tm.site_norbits[j])) || error("n or m not in range 1-site_norbits.")
    length(R) == 3 || error("R should be a 3-element vector.")

    sethopping!(tm, R, _to_orbital_index(tm, (i, p)), _to_orbital_index(tm, (j, q)), hopping)
end


function sethopping!(tm::TBModel, R::AbstractVector{<:Integer}, i::Int64, j::Int64, hopping::Matrix{<:Number})
    has_full_information(tm) || error("No site information is provided in the model.")
    ((i in 1:tm.nsites) && (j in 1:tm.nsites)) || error("i or j not in range 1-nsites.")
    size(hopping) == (tm.site_norbits[i], tm.site_norbits[j]) || error("hopping: wrong size")
    length(R) == 3 || error("R should be a 3-element vector.")

    for p in 1:size(hopping, 1), q in 1:size(hopping, 2)
        sethopping!(tm, R, (i, p), (j, q), hopping[p, q])
    end
end


function sethopping!(tm::TBModel{T}, R::AbstractVector{<:Integer}, hopping::Matrix{T}) where T <: Number
    size(hopping) == (tm.norbits, tm.norbits) || error("size of hopping matrix inconsistent with norbits.")
    length(R) == 3 || error("R should be a 3-element vector")
    R == R0 && maximum(abs.(hopping - hopping')) > HERM_TOL && error("hopping matrix $hopping at $R should be hermitian")

    if !(R in keys(tm.hoppings))
        tm.hoppings[R] = zeros(T, tm.norbits, tm.norbits)
        tm.hoppings[-R] = zeros(T, tm.norbits, tm.norbits)
    end

    if R == R0
        tm.hoppings[R] = (hopping + hopping') / 2
    else
        tm.hoppings[R] = hopping
        tm.hoppings[-R] = hopping'
    end
end


function addhopping!(tm::TBModel{T}, R::AbstractVector{<:Integer}, n::Int64, m::Int64, hopping::Number) where T <: Number
    ((n in 1:tm.norbits) && (m in 1:tm.norbits)) || error("n or m not in range 1-norbits.")
    length(R) == 3 || error("R should be a 3-element vector.")
    R == R0 && n == m && imag(hopping) > IMAG_TOL && error("On site energy should be real.")

    if !(R in keys(tm.hoppings))
        tm.hoppings[R] = zeros(T, tm.norbits, tm.norbits)
        tm.hoppings[-R] = zeros(T, tm.norbits, tm.norbits)
    end

    if R == R0 && n == m
        tm.hoppings[R][n, m] += real(hopping)
    else
        tm.hoppings[R][n, m] += hopping
        tm.hoppings[-R][m, n] += conj(hopping)
    end
end


function addhopping!(tm::TBModel, R::AbstractVector{<:Integer}, i::Int64, j::Int64, hopping::Matrix{<:Number})
    has_full_information(tm) || error("No site information is provided in the model.")
    ((i in 1:tm.nsites) && (j in 1:tm.nsites)) || error("i or j not in range 1-nsites.")
    size(hopping) == (tm.site_norbits[i], tm.site_norbits[j]) || error("hopping: wrong size")
    length(R) == 3 || error("R should be a 3-element vector.")

    if R == R0 && i == j
        for p in 1:tm.site_norbits[i], q in p:tm.site_norbits[i]
            addhopping!(tm, R, (i, p), (j, q), hopping[p, q])
        end
    else
        for p in 1:tm.site_norbits[i], q in 1:tm.site_norbits[j]
            addhopping!(tm, R, (i, p), (j, q), hopping[p, q])
        end
    end
end


@doc raw"""
```julia
addhopping!(tm::TBModel, R::AbstractVector{<:Integer}, hopping::Matrix{<:Number})
```

Add `hopping` to the Hamiltonian of `tm`.

`hopping` should be a `(tm.norbits, tm.norbits)` array containing the matrix
⟨0n|H|Rm⟩ with n and m indices. If R is [0, 0, 0], then `hopping` needs to
be Hermitian.
"""
function addhopping!(tm::TBModel, R::AbstractVector{<:Integer}, hopping::Matrix{<:Number})
    length(R) == 3 || error("R should be a 3-element vector.")
    size(hopping) == (tm.norbits, tm.norbits) || error("hopping not compatible with tm.")
    if R == R0
        ishermitian(hopping) || error("The Hamiltonian is not Hermition")
        for n in 1:tm.norbits, m in n:tm.norbits
            addhopping!(tm, R, n, m, hopping[n, m])
        end
    else
        for n in 1:tm.norbits, m in 1:tm.norbits
            addhopping!(tm, R, n, m, hopping[n, m])
        end
    end
end


function addhopping!(tm::TBModel{T}, R::AbstractVector{<:Integer}, (i, p)::Tuple{Int64,Int64},
    (j, q)::Tuple{Int64,Int64}, hopping::Number) where T <: Number
    has_full_information(tm) || error("No site information is provided in the model.")
    ((i in 1:tm.nsites) && (j in 1:tm.nsites)) || error("i or j not in range 1-nsites.")
    ((p in 1:tm.site_norbits[i]) && (q in 1:tm.site_norbits[j])) || error("n or m not in range 1-site_norbits.")
    length(R) == 3 || error("R should be a 3-element vector.")

    addhopping!(tm, R, _to_orbital_index(tm, (i, p)), _to_orbital_index(tm, (j, q)), hopping)
end


function gethopping(tm::TBModel, R::AbstractVector{<:Integer}, (i, p)::Tuple{Int64,Int64},
    (j, q)::Tuple{Int64,Int64})
    return tm.hoppings[R][_to_orbital_index(tm, (i, p)), _to_orbital_index(tm, (j, q))]
end


@doc raw"""
```julia
setoverlap!(tm::TBModel{T}, R::AbstractVector{<:Integer}, n::Int64, m::Int64,
    overlap::Number) where T
```

Set ⟨0n|Rm⟩ to `overlap` for `tm`.
"""
function setoverlap!(tm::TBModel{T}, R::AbstractVector{<:Integer}, n::Int64, m::Int64, overlap::Number) where T <: Number
    ((n in 1:tm.norbits) && (m in 1:tm.norbits)) || error("n or m not in range 1-norbits.")
    length(R) == 3 || error("R should be a 3-element vector")
    tm.isorthogonal && error("tm is orthogonal.")
    if R == R0 && n == m
        if abs(overlap) < 0.1
            error("An orbital should have substantial overlap with itself.")
        end
        if imag(overlap) > IMAG_TOL
            error("Overlap of one orbital with itself should be real.")
        end
    end

    # populate the matrix
    if !(R in keys(tm.overlaps))
        tm.overlaps[R] = zeros(T, tm.norbits, tm.norbits)
        tm.overlaps[-R] = zeros(T, tm.norbits, tm.norbits)
    end

    # set overlap
    if R == R0 && n == m
        tm.overlaps[R][n, m] = real(overlap)
    else
        tm.overlaps[R][n, m] = overlap
        tm.overlaps[-R][m, n] = conj(overlap)
    end

    # set position according to site_positions
    if R == R0 && n == m && has_full_information(tm)
        i = _to_site_index(tm, n)[1]
        setposition!(tm, R0, n, n, tm.site_positions[:, i] * tm.overlaps[R0][n, n])
    end
end


@doc raw"""
```julia
setposition!(
    tm::TBModel,
    R::AbstractVector{<:Integer},
    n::Int64,
    m::Int64,
    α::Int64,
    pos::Number;
    position_tolerance::Real=1.0e-4
)
```

Set ⟨0n|``r_α``|Rm⟩ to `pos` for `tm`.
Overlap matrices must be set before this method is called if `tm` is not orthogonal.

`position_tolerance` is used to check wehther the value is compatible with
`tm.site_positions` (if not missing).
"""
function setposition!(
    tm::TBModel{T},
    R::AbstractVector{<:Integer},
    n::Int64,
    m::Int64,
    α::Int64,
    pos::Number;
    position_tolerance::Real=1.0e-4
) where T <: Number
    ((n in 1:tm.norbits) && (m in 1:tm.norbits)) || error("n or m not in range 1-norbits.")
    length(R) == 3 || error("R should be a 3-element vector")
    (α in 1:3) || error("α not in 1-3.")
    R == R0 && n == m && imag(pos) > IMAG_TOL && error("Position of one orbital with itself should be real.")
    if R == R0 && n == m && has_full_information(tm)
        i = _to_site_index(tm, n)[1]
        if !isapprox(tm.site_positions[α, i] * tm.overlaps[R0][n, n], pos, atol=position_tolerance)
            error("pos incompatible with site_positions and overlaps.")
        end
    end

    # populate the matrix
    if !(R in keys(tm.positions))
        tm.positions[R] = [zeros(T, tm.norbits, tm.norbits) for i in 1:3]
        tm.positions[-R] = [zeros(T, tm.norbits, tm.norbits) for i in 1:3]
    end

    if R == R0 && n == m
        tm.positions[R][α][n, m] = real(pos)
    else
        tm.positions[R][α][n, m] = pos
        tmp = (-R in keys(tm.overlaps)) ? tm.overlaps[-R][m, n] : convert(T, 0.0)
        tm.positions[-R][α][m, n] = conj(pos) - (tm.lat * R)[α] * tmp
    end
end

function setposition!(
    tm::TBModel{T},
    R::AbstractVector{<:Integer},
    n::Int64,
    m::Int64, 
    pos::Vector{<:Number}
) where T <: Number
    ((n in 1:tm.norbits) && (m in 1:tm.norbits)) || error("n or m not in range 1-norbits.")
    length(R) == 3 || error("R should be a 3-element vector")
    length(pos) == 3 || error("pos should be a 3-element vector.")

    for α in 1:3
        setposition!(tm, R, n, m, α, pos[α])
    end
end


@doc raw"""
```julia
change_energy_reference(tm::TBModel, μ::Number)::TBModel
```

change zero energy reference to μ.
"""
function change_energy_reference(tm::TBModel, μ::Number)
    ntm = deepcopy(tm)
    if ntm.isorthogonal
        ntm.hoppings[[0, 0, 0]] -= μ * I
    else
        for (R, hopping) in tm.hoppings
            ntm.hoppings[R] -= μ * tm.overlaps[R]
        end
    end
    return ntm
end


function changebasis(tm::TBModel{ComplexF64}, Us::Dict{Int64,Matrix{ComplexF64}};
    set_canonical_ordered::Bool=false)
    ntm = deepcopy(tm)
    flatten_orbital_types = vcat(tm.orbital_types...)
    nspinlessorbits = tm.norbits ÷ (1 + tm.isspinful)
    U = zeros(ComplexF64, nspinlessorbits, nspinlessorbits)
    cnt = 0
    for l in flatten_orbital_types
        U[(cnt + 1):(cnt + 2l + 1), (cnt + 1):(cnt + 2l + 1)] = Us[l]
        cnt += 2l + 1
    end
    @assert cnt == nspinlessorbits
    if tm.isspinful U = kron([1 0; 0 1], U) end
    for (R, hopping) in tm.hoppings
        ntm.hoppings[R] = U * hopping * U'
    end
    for (R, overlap) in tm.overlaps
        ntm.overlaps[R] = U * overlap * U'
    end
    for (R, pos) in tm.positions
        ntm.positions[R] = [U * pos[α] * U' for α in 1:3]
    end
    if set_canonical_ordered
        set_is_canonical_ordered!(ntm, true)
    end
    return ntm
end


function _to_orbital_index(tm::TBModel, (i, n)::Tuple{Int64,Int64})
    if tm.isspinful
        tmp1 = div(sum(tm.site_norbits[1:(i - 1)]), 2)
        tmp2 = div(tm.site_norbits[i], 2)
        tmp3 = div(tm.norbits, 2)
        if n > tmp2
            return tmp3 + tmp1 + n - tmp2
        else
            return tmp1 + n
        end
    else
        return sum(tm.site_norbits[1:(i - 1)]) + n
    end
end


function _to_orbital_index(tm::TBModel, i::Int64)
    N = div(tm.norbits, (1 + tm.isspinful))
    l = div.(tm.site_norbits, (1 + tm.isspinful))
    r = sum(l[1:i - 1]) + 1:sum(l[1:i])
    if tm.isspinful
        return [r; r .+ N]
    else
        return collect(r)
    end
end


function _to_site_index(tm::TBModel, n::Int64)
    l = div.(tm.site_norbits, 1 + tm.isspinful)
    N = div(tm.norbits, 1 + tm.isspinful)
    s = 0
    if n > N
        n = n - N
        s = 1
    end
    i = 0
    ss = 0
    for cnt in 1:length(l)
        ss += l[cnt]
        if n <= ss
            i = cnt
            break
        end
    end
    return (i, n - sum(l[1:(i - 1)]) + s * l[i])
end

function get_orbital_position(tm::TBModel, m::Int64)
    position = zeros(Float64, 3)
    for α = 1:3
        position[α] = tm.positions[R0][α][m,m]
    end
    return position
end

function get_orbital_position(tm::TBModel)
    n = tm.norbits
    position = zeros(Float64, (3, n))
    for m = 1:n
        for α = 1:3
            position[α,m] = tm.positions[R0][α][m,m]
        end
    end
    return position
end

function prune!(tm::TBModel, tol::Float64=1.0e-10)
    for (R, hopping) in tm.hoppings
        if !(R == [0, 0, 0]) && all(abs.(hopping) .< tol)
            pop!(tm.hoppings, R)
        end
    end
    for (R, overlap) in tm.overlaps
        if !(R == [0, 0, 0]) && all(abs.(overlap) .< tol)
            pop!(tm.overlaps, R)
        end
    end
    for (R, pos) in tm.positions
        if !(R == [0, 0, 0]) && all([all(abs.(pos[α]) .< tol) for α in 1:3])
            pop!(tm.positions, R)
        end
    end
end
