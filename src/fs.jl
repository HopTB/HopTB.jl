using LinearAlgebra, StaticArrays

export FermiSurface

################################################################################
##  Fermi Surface Types
################################################################################


@doc raw"""
FermiSurface encodes information about one Fermi surface. It is usually constructed from Marching
Tetrahedra method.

# Fields
  - `rlat::SMatrix{3,3,Float64,9}`, the reciprocal lattice vectors stored in columns
  - `ks::Matrix{Float64}`, a list of reduced k points on the Fermi surface.
  - `faces::Matrix{Int64}`, a set of triangles connecting neighbouring k points on
    the Fermi surface.
  - `weights::Vector{Float64}`, the area associated with each k point.
  - `bandidx::Union{Missing,Int64}`, the band index.

# Construction method
```julia
FermiSurface(rlat::AbstractMatrix{<:Real}, ks::AbstractMatrix{<:Real}, faces::AbstractMatrix{<:Integer})
```
"""
struct FermiSurface
    rlat::SMatrix{3,3,Float64,9}
    ks::Matrix{Float64}
    faces::Matrix{Int64}
    weights::Vector{Float64}
    bandidx::Union{Missing,Int64}
end

function _get_triangle_area(vertices)
    # Heron's formula
    a = norm(vertices[:, 1] - vertices[:, 2])
    b = norm(vertices[:, 2] - vertices[:, 3])
    c = norm(vertices[:, 3] - vertices[:, 1])
    s = (a + b + c) / 2
    return √(s * (s - a) * (s - b) * (s - c))    
end

function FermiSurface(
    rlat::AbstractMatrix{<:Real},
    ks::AbstractMatrix{<:Real},
    faces::AbstractMatrix{<:Integer};
    bandidx::Union{Missing,Int64}=missing
)
    nks = size(ks, 2)
    nfaces = size(faces, 2)
    (minimum(faces) >= 0 && maximum(faces) <= nks) || error("Wrong face indices.")
    weights = zeros(nks)
    for iface in 1:nfaces
        face = faces[:, iface]
        area = _get_triangle_area(rlat * ks[:, face])
        for i in 1:3
            weights[face[i]] += area / 3
        end
    end
    return FermiSurface(rlat, ks, faces, weights, bandidx)
end

function Base.show(io::IO, fs::FermiSurface)
    print(io, "FermiSurface: band $(fs.bandidx)")
end
