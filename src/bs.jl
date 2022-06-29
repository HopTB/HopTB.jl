module BandStructure

using Distributed, LinearAlgebra, SharedArrays
using StaticArrays, Meshing
using Base.Threads: nthreads, threadid, @threads
using ..HopTB
using ..HopTB.Meshes: UniformMesh, find_point_in_mesh
using ..HopTB.Utilities: constructmeshkpts, distance_on_circle, splitkpts, splitvector
using ..HopTB.Parallel: ParallelFunction, claim!, stop!, parallel_do

export clteig, getbs, getjdos, getdos
export get_fermi_surfaces


################################################################################
##  Energy mesh
################################################################################

"""
A uniform mesh of band energies in the Brillouin zone.
"""
struct EnergyMesh
    energies::Array{Float64,4}
    band_indices::Vector{Int64}
    meshsize::Vector{Int64}
    endboundary::Bool
end

function EnergyMesh(
    energies::Array{Float64,4},
    band_indices::AbstractVector{Int64},
    meshsize::Vector{Int64};
    endboundary::Bool=false
)
    size(energies, 1) == length(band_indices) || error("energies is not compatible with band_indices.")
    collect(size(energies)[2:end]) == meshsize || error("energies is not compatible with meshsize.")
    length(unique(band_indices)) ==  length(band_indices) || error("band_indices contains repeated indices.")
    return EnergyMesh(energies, band_indices, meshsize, endboundary)
end

function Base.show(io::IO, em::EnergyMesh)
    print(io, "EnergyMesh: $(length(em.band_indices)) band(s) on a $(em.meshsize) mesh")
end

@doc raw"""
```julia
get_energy_mesh(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    band_indices::AbstractVector{Int64};
    endboundary::Bool=false,
    batchsize::Int64=1
)::EnergyMesh
```

Collect the eigenvalues for bands indicated by `band_indices` on a mesh denoted by
`meshsize`. `Es[n, i, j, k]` correspond to the energy of the band n at the k point denoted by 
(i, j, k). If `endboundary` is true, the end boundary of the Brillouin zone in included.
"""
function get_energy_mesh(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    band_indices::AbstractVector{Int64};
    endboundary::Bool=false,
    batchsize::Int64=1
)
    energies = SharedArray{Float64}(length(band_indices), meshsize...)
    mesh = UniformMesh(meshsize, endboundary=endboundary)
    parallel_do(
        k -> energies[:, find_point_in_mesh(k, mesh)...] = geteig(tm, k).values[band_indices],
        mesh;
        batchsize=batchsize
    )
    return EnergyMesh(sdata(energies), band_indices, meshsize; endboundary=endboundary)
end


"""
```julia
remove_endboundary(em::EnergyMesh)
```
Remove the endboundary from `em`.
"""
function remove_endboundary(em::EnergyMesh)
    em.endboundary || error("The energy mesh does not contain end boundaries.")
    return EnergyMesh(em.energies[:, 1:end-1, 1:end-1, 1:end-1],
        deepcopy(em.band_indices), em.meshsize .- 1; endboundary=false)
end


"""
```julia
add_endboundary(em::EnergyMesh)
```
Add the endboundary to `em`.
"""
function add_endboundary(em::EnergyMesh)
    !(em.endboundary) || error("The energy mesh already contains end boundaries.")
    N, M, L = em.meshsize
    energies = zeros(length(em.band_indices), N + 1, M + 1, L + 1)
    for k in 1:L+1, j in 1:M+1, i in 1:N+1
        energies[:, i, j, k] = em.energies[:, mod(i, 1:N), mod(j, 1:M), mod(k, 1:L)]
    end
    return EnergyMesh(energies, deepcopy(em.band_indices), em.meshsize .+ 1; endboundary=true)
end


################################################################################
##  Original Code
################################################################################

function _clteig(atm::AbstractTBModel, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)
    allegvals = zeros((atm.norbits, nkpts))
    for ik in 1:nkpts
        k = kpts[:, ik]
        allegvals[:, ik] = geteig(atm, k).values
    end
    return allegvals
end


function clteig(atm::AbstractTBModel, kpts::AbstractMatrix{Float64})
    nkpts = size(kpts, 2)

    kptslist = HopTB.Utilities.splitkpts(kpts, nworkers())
    jobs = Vector{Future}()
    for iw in 1:nworkers()
        job = @spawn _clteig(atm, kptslist[iw])
        append!(jobs, [job])
    end

    allegvals = zeros((atm.norbits, 0))
    for iw in 1:nworkers()
        allegvals = cat(allegvals, HopTB.Utilities.safe_fetch(jobs[iw]), dims=(2,))
    end
    return allegvals
end


function clteig_worker(tm::AbstractTBModel, bandinds::Vector{Int64}, kinds::Vector{CartesianIndex{3}}, Es::SharedArray{Float64})
    nks = length(kinds)
    for kind in kinds
        k = [(kind.I[i]-1)/size(Es, i+1) for i in 1:3]
        Es[:, kind] = geteig(tm, k).values[bandinds]
    end
end

"""
```julia
clteig(tm::AbstractTBModel, nkmesh::Vector{Int64}; bandinds::Vector{Int64}=collect(1:tm.norbits))
```

Collect eigenvalues of bands denoted by `bandinds` on the kmesh defined by `nkmesh`.
This function returns a 4D SharedArray:
    first index -> bands;
    second-fourth index -> k point.

For example, with `nkmesh` being `[10, 10, 10]`, `Es[3, 4, 5, 6]` is the third eigenvalue
at k point `[0.3, 0.4, 0.5]`.
"""
function clteig(tm::AbstractTBModel, nkmesh::Vector{Int64}; bandinds::Vector{Int64}=collect(1:tm.norbits))
    nks = prod(nkmesh)
    nbands = length(bandinds)
    Es = SharedArray{Float64}((nbands, nkmesh...))
    kinds_list = splitvector(CartesianIndices(view(Es, 1, :, :, :))[:], nworkers())
    jobs = Vector{Future}()
    for iw in 1:nworkers()
        job = @spawn clteig_worker(tm, bandinds, kinds_list[iw], Es)
        append!(jobs, [job])
    end
    for job in jobs fetch(job) end
    return Es
end


@doc raw"""
```julia
getbs(atm::AbstractTBModel, kpath::AbstractMatrix{Float64}, pnkpts::Int64;
    connect_end_points::Bool=false)::(Vector{Float64}, Matrix{Float64})
```

Calculate band structure along a `kpath` for `atm`. `kpath` is a matrix that marks start
and end points of each line in columns, and `pnkpts` is the number of points in each line.
if `connect_end_points` is true, end point of the previous segment is used
as the start point of the next segment.
This function returns (`kdist`, `egvals`), where `kdist` is the distance of k points in
reciprocal space and `egvals` contains band energies stored in column for each k point.
"""
function getbs(atm::AbstractTBModel, kpath::AbstractMatrix{Float64}, pnkpts::Int64;
    connect_end_points::Bool=false)
    kdist, kpts = HopTB.Utilities.constructlinekpts(kpath,
        convert(Matrix{Float64}, atm.rlat), pnkpts, connect_end_points=connect_end_points)
    nkpts = size(kpts, 2)
    egvals = zeros(atm.norbits, nkpts)
    for ik in 1:nkpts
        k = kpts[:, ik]
        egvals[:, ik] = geteig(atm, k).values
    end
    return (kdist, egvals)
end


function _check_in_ks(k, ks)
    for ik in 1:size(ks, 2)
        if norm(distance_on_circle.(ks[:, ik], k)) < 1.0e-4
            return true
        end
    end
    return false
end


function jdos_worker(ks::Matrix{Float64}, tm::AbstractTBModel,
    ωs::Vector{Float64}, μ::Float64, ϵ::Float64=0.1)
    nωs = size(ωs, 1)
    jdos = zeros(nωs)
    nks = size(ks, 2)
    for ik in 1:nks
        k = ks[:, ik]
        egvals = geteig(tm, k).values
        for n in 1:tm.norbits, m in 1:tm.norbits
            en = egvals[n]
            em = egvals[m]
            fn = (en < μ) ? 1 : 0
            fm = (em < μ) ? 1 : 0
            for iω in 1:nωs
                ω = ωs[iω]
                jdos[iω] += (fm-fn)*exp(-(en-em-ω)^2/ϵ^2)*sqrt(1/π)/ϵ
            end
        end
    end
    return jdos
end


@doc raw"""
```julia
getjdos(tm::AbstractTBModel, ωs::Vector{Float64}, μ::Float64,
    nkmesh::Vector{Int64}; ϵ::Float64=0.1)
```

Calculate joint density of states between valence and conduction band.

Joint density of states is defined as
```math
∫\frac{d\mathbf{k}}{(2π)^3}\sum_{n,m}f_{mn}δ(E_n-E_m-ω).
```
"""
function getjdos(tm::AbstractTBModel, ωs::Vector{Float64}, μ::Float64,
    nkmesh::Vector{Int64}; ϵ::Float64=0.1)
    nks = prod(nkmesh)
    nωs = size(ωs, 1)
    klist = splitkpts(constructmeshkpts(nkmesh), nworkers())
    pf = ParallelFunction(jdos_worker, tm, ωs, μ, ϵ, len=nworkers())
    for iworker in 1:nworkers()
        pf(klist[iworker])
    end
    jdos = zeros(nωs)
    for iworker in 1:nworkers()
        jdos += claim!(pf)
    end
    stop!(pf)
    bzvol = abs(tm.rlat[:, 1]×tm.rlat[:, 2]⋅tm.rlat[:, 3])
    return jdos*bzvol/nks/(2π)^3
end


@doc raw"""
```julia
clt_jdos(tm::AbstractTBModel, ωs::Vector{Float64}, μ::Float64,
    ks::Matrix{Float64}; ϵ::Float64=0.1)
```

Collect contributions of joint density of states between valence and
conduction band at each k point.

The contribution of joint density of states is defined as
```math
\sum_{n,m}f_{mn}δ(En-Em-ω).
```
"""
function clt_jdos(tm::AbstractTBModel, ωs::Vector{Float64}, μ::Float64,
    ks::AbstractMatrix{Float64}; ϵ::Float64=0.1)
    nωs = size(ωs, 1)
    nks = size(ks, 2)
    jdos = zeros(nks, nωs)
    for ik in 1:nks
        k = ks[:, ik]
        egvals = geteig(tm, k).values
        for n in 1:tm.norbits, m in 1:tm.norbits
            en = egvals[n]
            em = egvals[m]
            fn = (en < μ) ? 1 : 0
            fm = (em < μ) ? 1 : 0
            for iω in 1:nωs
                ω = ωs[iω]
                jdos[ik, iω] += (fm-fn)*exp(-(en-em-ω)^2/ϵ^2)*sqrt(1/π)/ϵ
            end
        end
    end
    return jdos
end


################################################################################
##  Density of states
################################################################################

function getdos_k!(
    dos::Vector{Float64},
    tm::AbstractTBModel,
    ωs::AbstractVector{Float64},
    k::Vector{Float64},
    ϵ::Float64=0.1,
)
    Es = geteig(tm, k).values
    factor = 1/((2π)^3*ϵ*√π)
    for (iω, ω) in enumerate(ωs), n in 1:tm.norbits
        dos[iω] += exp(-(Es[n]-ω)^2/ϵ^2)*factor
    end
end

function get_dos_ks(
    ks::Matrix{Float64},
    tm::AbstractTBModel,
    ωs::AbstractVector{Float64},
    ϵ::Float64=0.1,
)
    nks = size(ks, 2); nωs = length(ωs)
    dos = zeros(nωs)
    for i = 1:nks
        getdos_k!(dos, tm, ωs, ks[:, i], ϵ)
    end
    return dos
end

@doc raw"""
```julia
getdos(tm::AbstractTBModel, ωs::AbstractVector{Float64}, nkmesh::Vector{Int64};
    ϵ::Float64=0.1) --> dos::Vector{Float64}
```

Calculate density of states by a `nkmesh` mesh in reciprocal space.

The density of states is defined as
```math
ν(ω) = \sum_n \int \frac{\mathrm{d}\boldsymbol{k}}{(2\pi)^3} δ(E_{n\boldsymbol{k}}-ω).
```

`ϵ` is the broadening of the delta function.

`dos` is in the unit of 1/(eV*Å^3).
"""
function getdos(
    tm::AbstractTBModel,
    ωs::AbstractVector{Float64},
    nkmesh::Vector{Int64};
    ϵ::Float64=0.1
)
    nks = prod(nkmesh); nωs = size(ωs, 1)
    ks = constructmeshkpts(nkmesh)
    nsplits = 10*nworkers()
    ks_list = HopTB.Utilities.splitkpts(ks, nsplits)
    dos = zeros(nωs)

    pf = ParallelFunction(get_dos_ks, tm, ωs, ϵ, len=nsplits)
    @async for i in 1:nsplits
        pf(ks_list[i])
    end
    @sync @async for i in 1:nsplits
        dos += claim!(pf)
    end
    stop!(pf)

    bzvol = abs(det(tm.rlat))
    return dos*bzvol/nks
end


"""
```julia
getdos(
    tm::AbstractTBModel,
    em::EnergyMesh,
    ωs::AbstractVector{Float64};
    ϵ::Float64=0.1
)
```

Compute density of states from an energy mesh.
"""
function getdos(
    tm::AbstractTBModel,
    em::EnergyMesh,
    ωs::AbstractVector{Float64};
    ϵ::Float64=0.1
)
    Es = em.endboundary ? remove_endboundary(em).energies : em.energies
    dos = zeros(length(ωs))
    for (iω, ω) in enumerate(ωs), E in Es
        dos[iω] += exp(-(E-ω)^2 / ϵ^2)
    end
    bzvol = abs(det(tm.rlat))
    nks = prod(size(Es)[2:end])
    return dos * bzvol / ((2π)^3 * ϵ * √π * nks)
end


"""
```julia
get_fermi_energy(
    tm::AbstractTBModel,
    ωs::AbstractRange,
    dos::AbstractVector{<:Real},
    n::Integer
)
```

Calculate fermi energy for `tm` from the density of states `dos`. The number
of electrons in a unit cell is `n`. `ωs` is the energies corresponding to `dos`.
"""
function get_fermi_energy(
    tm::AbstractTBModel,
    ωs::AbstractRange,
    dos::AbstractVector{<:Real},
    n::Integer
)
    n > 0 || error("The number of electrons should be larger than zero.")
    length(ωs) > 2 || error("The length of ωs should be larger than zero.")
    ucvol = det(tm.lat)
    dω = ωs[2] - ωs[1]
    s = 0.0
    for i in 1:length(dos)
        s += dos[i]
        if s * ucvol * dω > n
            return ωs[i]
        end
    end
    error("The integration of dos is smaller than the number of electrons.")
end


################################################################################
##  Extract Fermi surface
################################################################################

"""
```julia
function get_fermi_surfaces(
    tm::AbstractTBModel,
    em::EnergyMesh;
    fermi_energy::Real=0.0,
    band_indices::AbstractVector{Int64}=em.band_indices
)::Vector{FermiSurface}
```

Calculate Fermi surfaces from the energy mesh `em` for bands designated with `band_indices`.

The Fermi surfaces are extracted by Marching Tetrahedra method.
"""
function get_fermi_surfaces(
    tm::AbstractTBModel,
    em::EnergyMesh;
    fermi_energy::Real=0.0,
    band_indices::AbstractVector{Int64}=em.band_indices
)
    issubset(band_indices, em.band_indices) || error("The energy mesh does not contain required bands.")
    Es = em.endboundary ? em.energies : add_endboundary(em).energies
    results = Vector{FermiSurface}()
    for bandidx in band_indices
        ks, faces = isosurface(
            Es[findfirst(isequal(bandidx), em.band_indices), :, :, :],
            MarchingTetrahedra(iso=fermi_energy); 
            origin=[0.0, 0.0, 0.0],
            widths=[1.0, 1.0, 1.0]
        )
        if !isempty(ks)
            push!(results, FermiSurface(tm.rlat, reduce(hcat, ks), reduce(hcat, faces); bandidx=bandidx))
        end
    end
    return results
end


"""
```julia
function get_fermi_surfaces(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    band_indices::AbstractVector{Int64};
    fermi_energy::Real=0.0,
    batchsize::Int64=1
)::Vector{FermiSurface}
```

Calculate Fermi surfaces (with Fermi energy `fermi_energy`) for bands (specified by `band_indices`) for `tm`.
The Fermi surfaces are extracted by Marching Tetrahedra method on a uniform mesh specified by `meshsize`.
"""
function get_fermi_surfaces(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    band_indices::AbstractVector{Int64};
    fermi_energy::Real=0.0,
    batchsize::Int64=1
)
    em = get_energy_mesh(tm, meshsize, band_indices; batchsize=batchsize, endboundary=true)
    return get_fermi_surfaces(tm, em; fermi_energy=fermi_energy)
end

"""
```julia
function get_fermi_surface(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    bandidx::Integer;
    fermi_energy::Real=0.0,
    batchsize::Int64=1
)::Union{FermiSurface,Nothing}
```

Calculate Fermi surfaces (with Fermi energy `fermi_energy`) for the band specified by `bandidx` for `tm`.
The Fermi surfaces are extracted by Marching Tetrahedra method on a uniform mesh specified by `meshsize`.
"""
function get_fermi_surface(
    tm::AbstractTBModel,
    meshsize::AbstractVector{Int64},
    bandidx::Integer;
    fermi_energy::Real=0.0,
    batchsize::Int64=1
)
    fss = get_fermi_surfaces(tm, meshsize, [bandidx]; fermi_energy=fermi_energy, batchsize=batchsize)
    if isempty(fss)
        return nothing
    else
        return fss[1]
    end
end

end
