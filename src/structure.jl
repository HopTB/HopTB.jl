module Structure

using ..HopTB
using LinearAlgebra

"""
```julia
permute_orbits!(tm::TBModel, p::Vector{Int64})
```

Permute the orbits of `tm` by a permutation vector `p`.

Example of `p`: `p=[2, 4, 3, 1]` denotes the first orbit of the permuted `tm`
is the second orbit of the original `tm`; the second orbit of the permuted `tm`
is the fourth orbit of the original `tm`...

After the orbits are permuted, any field of `tm` that is allowd to be missing
is set to missing. This is due to a permutation of the orbitals generally loses
site information.
"""
function permute_orbits!(tm::TBModel, p::Vector{Int64})
    @assert length(p) == tm.norbits
    P = zeros(Int64, tm.norbits, tm.norbits)
    for i in 1:tm.norbits
        P[i, p[i]] = 1
    end
    for (R, hopping) in tm.hoppings
        tm.hoppings[R] = P*hopping*P'
    end
    for (R, overlap) in tm.overlaps
        tm.overlaps[R] = P*overlap*P'
    end
    for (R, position) in tm.positions
        tm.positions[R] = [P*position[α]*P' for α in 1:3]
    end
    HopTB.remove_extra_information!(tm)
    return tm
end


"""
```julia
shift!(tm, d::Vector{Float64}=[0.0, 0.0, 0.0])
```

shift `tm` by distance `d`.

`d` is in reduced coordinates.
"""
function shift!(tm, d::Vector{Float64}=[0.0, 0.0, 0.0])
    dc = tm.lat*d
    for (R, position) in tm.positions
        tm.positions[R] = [position[α]+dc[α]*tm.overlaps[R] for α in 1:3]
    end
    if !ismissing(tm.site_positions)
        tm.site_positions = tm.site_positions.+dc
    end
end


"""
```julia
makeslab(tm::TBModel, dir::Int64, nlayer::Int64)
```

Constructs a (d-1)-dimensional tight-binding model out of a d-dimensional one
by repeating the unit cell a given number of times (`nlayer`) along one of the
periodic lattice vectors (`dir`).
"""
function makeslab(tm::TBModel, dir::Int64, nlayer::Int64)
    dir in 1:3 || error("dir should be an integer between 1 to 3.")
    nlayer > 0 || error("nlayer should be positive.")

    lat_slab = convert(Array, tm.lat)
    lat_slab[:, dir] = lat_slab[:, dir]*nlayer

    slab = TBModel(tm.norbits*nlayer, lat_slab, isorthogonal=tm.isorthogonal)

    for (R, hopping) in tm.hoppings
        if abs(R[dir]) < nlayer
            R′ = convert(Array, R)
            R′[dir] = 0
            for n in 1:tm.norbits, m in 1:tm.norbits
                for i in 0:nlayer-1
                    n′ = n+i*tm.norbits
                    m′ = m+i*tm.norbits+R[dir]*tm.norbits
                    if m′ in 1:slab.norbits
                        sethopping!(slab, R′, n′, m′, hopping[n, m])
                    end
                end
            end
        end
    end

    if !tm.isorthogonal
        for (R, overlap) in tm.overlaps
            if abs(R[dir]) < nlayer
                R′ = convert(Array, R)
                R′[dir] = 0
                for n in 1:tm.norbits, m in 1:tm.norbits
                    for i in 0:nlayer-1
                        n′ = n+i*tm.norbits
                        m′ = m+i*tm.norbits+R[dir]*tm.norbits
                        if m′ in 1:slab.norbits
                            setoverlap!(slab, R′, n′, m′, overlap[n, m])
                        end
                    end
                end
            end
        end
    end

    for (R, position) in tm.positions
        if abs(R[dir]) < nlayer
            overlap = get(tm.overlaps, R, zeros(tm.norbits, tm.norbits))
            R′ = convert(Array, R)
            R′[dir] = 0
            for n in 1:tm.norbits, m in 1:tm.norbits
                for i in 0:nlayer-1
                    n′ = n+i*tm.norbits
                    m′ = m+i*tm.norbits+R[dir]*tm.norbits
                    if m′ in 1:slab.norbits
                        for α in 1:3
                            setposition!(slab, R′, n′, m′, α,
                                position[α][n, m]+i*tm.lat[α, dir]*overlap[n, m])
                        end
                    end
                end
            end
        end
    end

    if !ismissing(tm.isspinful)
        if tm.isspinful
            p = vcat([collect(1:tm.norbits÷2).+tm.norbits*(i-1) for i in 1:nlayer]...)
            p = vcat(p, p.+tm.norbits÷2)
            permute_orbits!(slab, p)
        end
        slab.isspinful = tm.isspinful
    end

    if !ismissing(tm.nsites)
        slab.nsites = tm.nsites*nlayer
        slab.site_norbits = vcat([tm.site_norbits for i in 1:nlayer]...)
        slab.site_positions = hcat([tm.site_positions.+(i-1)*tm.lat[:, dir] for i in 1:nlayer]...)
    end

    if !ismissing(tm.orbital_types)
        slab.orbital_types = vcat([tm.orbital_types for i in 1:nlayer]...)
    end

    slab.is_canonical_ordered = tm.is_canonical_ordered

    return slab
end


"""
```julia
make_superlattice(tm::TBModel, sc::Matrix{Int64})
```

Construct superlattice. `sc` is the supercell lattice vector expressed
in the basis of the original lattice vector stored in column.

Returns `(Rs, orbindices)`

Every new orbital corresponds to the original orbital `orbindices[i]` at `Rs[:, i]`.

Notices:
 (a) For the original model `tm`, all components of the reduced coordinates of
  every orbital should have absolute value <= 1.
 (b) If the original model has `isspinful=true`, it will be guaranteed for the
  new superlattice, first half orbitals are spin up and second half are spin
  down and the two halves are in one to one correspondence.
"""
function make_superlattice(tm::TBModel, sc::AbstractMatrix{Int64})
    nvol = round(Int64, det(sc))
    norbits_sc = tm.norbits*nvol
    lat_sc = tm.lat*sc
    orbpos = HopTB.get_orbital_position(tm)
    Rs = zeros(Int64, (3, norbits_sc))
    orbindices = zeros(Int64, norbits_sc)
    cnt = 0
    nspins = 1
    if !ismissing(tm.isspinful)
        if tm.isspinful
            nspins = 2
        end
    end
    for i=-1:maximum(sc[1, :])+1, j=-1:maximum(sc[2, :])+1, k=-1:maximum(sc[3, :])+1
        for iorb=1:tm.norbits÷nspins
            pos = orbpos[:, iorb]+tm.lat*[i, j, k]
            (x, y, z) = inv(lat_sc)*pos # reduced coordinate in lat_sc
            if 0<=x<1.0 && 0<=y<1.0 && 0<=z<1.0
                cnt += 1
                Rs[:, cnt] = [i, j, k]
                orbindices[cnt] = iorb
            end
        end
    end
    cnt == norbits_sc÷nspins || error("Failing to find orbitals.")
    if nspins == 2
        tmp = norbits_sc÷2
        Rs[:, tmp+1:end] = Rs[:, 1:tmp]
        orbindices[tmp+1:end] = orbindices[1:tmp].+tm.norbits÷nspins
    end

    return Rs, orbindices
end


function _find_sc_ind(R, n, Rs, orbinds, sc)
    for i in 1:size(Rs, 2)
        if n == orbinds[i]
            tmp = inv(sc)*(R-Rs[:, i])
            R′ = round.(Int64, tmp)
            if norm(R′-tmp) < 1.0e-8
                return (R′, i)
            end
        end
    end
    error("Cannot find the orbit")
end


"""
```julia
make_supercell(tm, sc::AbstractMatrix{Int64})
```

Construct supercell TBModel. `sc` is the supercell lattice vector expressed
in the basis of the original lattice vector stored in column.

Notices:
 (a) For the original model `tm`, all components of the reduced coordinates of
  every orbital should have absolute value <= 1.
 (b) If the original model has `isspinful=true`, it will be guaranteed for the
  new superlattice, first half orbitals are spin up and second half are spin
  down and the two halves are in one to one correspondence.
 (c) All extra information except for `isspinful` is lost during the construction of
  the supercell.
"""
function make_supercell(tm::TBModel{T}, sc::AbstractMatrix{Int64}) where T<:Number
    Rs, orbinds = make_superlattice(tm, sc)
    nvol = round(Int64, det(sc))
    tm_sc = TBModel(tm.norbits*nvol, tm.lat*sc, isorthogonal=tm.isorthogonal)

    for n in 1:tm_sc.norbits
        n′ = orbinds[n]
        for (R′, hopping) in tm.hoppings
            for m′ in 1:tm.norbits
                R, m = _find_sc_ind(Rs[:, n]+R′, m′, Rs, orbinds, sc)
                sethopping!(tm_sc, R, n, m, hopping[n′, m′])
            end
        end
    end

    if !tm.isorthogonal
        for n in 1:tm_sc.norbits
            n′ = orbinds[n]
            for (R′, overlap) in tm.overlaps
                for m′ in 1:tm.norbits
                    R, m = _find_sc_ind(Rs[:, n]+R′, m′, Rs, orbinds, sc)
                    setoverlap!(tm_sc, R, n, m, overlap[n′, m′])
                end
            end
        end
    end


    for n in 1:tm_sc.norbits
        n′ = orbinds[n]
        for (R′, position) in tm.positions
            for m′ in 1:tm.norbits
                R, m = _find_sc_ind(Rs[:, n]+R′, m′, Rs, orbinds, sc)
                if R′ in keys(tm.overlaps)
                    tmp = tm.overlaps[R′][n′, m′]
                else
                    tmp = zero(T)
                end
                for α in 1:3
                    setposition!(tm_sc, R, n, m, α, position[α][n′, m′]+(tm.lat*Rs[:, n])[α]*tmp)
                end
            end
        end
    end

    tm_sc.isspinful=tm.isspinful

    return tm_sc
end


end
