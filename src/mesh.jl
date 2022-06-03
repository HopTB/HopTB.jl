module Meshes

export UniformMesh, find_point_in_mesh

@doc raw"""
```julia
UniformMesh(
    meshsize::Vector{Int64};
    startpoint::Vector{Float64}=zeros(length(meshsize)),
    regionsize::Vector{Float64}=ones(length(meshsize)),
    endboundary::Bool=false
)
```

A struct representing uniform mesh without storing everything in memory. If `endboundary` are true,
the mesh contains points in the end boundary (start boundary is always included). In the mesh, the
first dimension changes fastest.

Example:

```julia
julia> m = Hop.Meshes.UniformMesh([3, 4, 5]; startpoint=[0.1, 0.2, 0.3], regionsize=[1.0, 2.0, 3.0])
UniformMesh: [3, 4, 5]

julia> m[29]
3-element Array{Float64,1}:
 0.43333333333333335
 0.7
 1.5

julia> sm = m[3:5]
UniformSubMesh: 3 points in a [3, 4, 5] mesh

julia> sm[3]
3-element Array{Float64,1}:
 0.43333333333333335
 0.7
 0.3
```
"""
struct UniformMesh
    meshsize::Vector{Int64}
    startpoint::Vector{Float64}
    regionsize::Vector{Float64}
    dimension::Int64
    endboundary::Bool
    meshsize_cumprod::Vector{Int64}
end

function UniformMesh(
    meshsize::Vector{Int64};
    startpoint::Vector{Float64}=zeros(length(meshsize)),
    regionsize::Vector{Float64}=ones(length(meshsize)),
    endboundary::Bool=false
)
    length(meshsize) == length(startpoint) == length(regionsize) || 
        error("Arguments are of incompatible shapes.")
    all(meshsize .> 0) || error("meshsize should contain positive integers.")
    all(regionsize .> 0) || error("regionsize should contain positive numbers.")
    return UniformMesh(meshsize, startpoint, regionsize, length(meshsize), endboundary, cumprod(meshsize))
end

function Base.show(io::IO, mesh::UniformMesh)
    print(io, "UniformMesh: $(mesh.meshsize)")
end

function Base.length(m::UniformMesh)
    return m.meshsize_cumprod[end]
end

function Base.getindex(m::UniformMesh, idx::Int64)
    (idx >= 1 && idx <= prod(m.meshsize)) || error("Index is out of range.")
    point = zeros(Float64, m.dimension)
    remaining_idx = idx - 1
    for i in m.dimension:-1:2
        foo, remaining_idx = divrem(remaining_idx, m.meshsize_cumprod[i-1])
        point[i] = m.regionsize[i] * foo / (m.meshsize[i] - m.endboundary) + m.startpoint[i]
    end
    point[1] = m.regionsize[1] * remaining_idx / (m.meshsize[1] - m.endboundary) + m.startpoint[1]
    return point
end


struct UniformSubMesh{T<:AbstractArray{Int64}}
    mesh::UniformMesh
    indices::T
end

function Base.show(io::IO, submesh::UniformSubMesh)
    print(io, "UniformSubMesh: $(length(submesh.indices)) points in a $(submesh.mesh.meshsize) mesh")
end

function Base.getindex(mesh::UniformMesh, indices::AbstractVector{Int64})
    length(indices) > 0 || error("indices should have at least one elements.")
    return UniformSubMesh{typeof(indices)}(mesh, indices)
end

function Base.length(submesh::UniformSubMesh)
    return length(submesh.indices)
end

function Base.getindex(submesh::UniformSubMesh, idx::Int64)
    return submesh.mesh[submesh.indices[idx]]
end

function Base.iterate(mesh::Union{UniformSubMesh,UniformMesh})
    return (mesh[1], 2)
end

function Base.iterate(mesh::Union{UniformSubMesh,UniformMesh}, state::Int64)
    if state > length(mesh)
        return nothing
    else
        return (mesh[state], state + 1)
    end
end


@doc raw"""
```julia
find_point_in_mesh(point::AbstractVector{<:Real}, mesh::UniformMesh) --> Vector{Int64}
```

Given a point in a mesh, find its position in the mesh.

Example
```julia
julia> mesh = Hop.Meshes.UniformMesh([3, 4, 5]; startpoint=[0.1, 0.2, 0.3], regionsize=[1.0, 2.0, 3.0]);

julia> find_point_in_mesh(mesh[29], mesh)
3-element Array{Int64,1}:
 2
 2
 3
```
"""
function find_point_in_mesh(point::AbstractVector{<:Real}, mesh::UniformMesh)
    mesh.dimension == length(point) || error("Length of the point is not compatible with the mesh.")
    ret = zeros(Int64, mesh.dimension)
    for i in 1:mesh.dimension
        foo = (point[i] - mesh.startpoint[i]) * (mesh.meshsize[i] - mesh.endboundary) / mesh.regionsize[i]
        bar = round(Int64, foo)
        if abs(foo - bar) > 1.0e-4
            error("The point is not in the mesh.")
        end
        ret[i] = bar + 1
    end
    return ret 
end

end