module Group

import Base.*, Base.==, Base.show
using LinearAlgebra
using ..HopTB, ..HopTB.Memoize
export Symmetry
export getrotation, getTRS, getinversion, getmirror
export generategroup, symmetrize, get_bloch_rep


struct Symmetry
    rotation_matrix::Matrix{Float64}
    translation::Vector{Float64}
    spin_rep::Matrix{ComplexF64}
    isantilinear::Bool
    isdouble::Bool
end

"""
```julia
Symmetry(rotation_matrix, translation, spin_rep;
    isantilinear::Bool=false, isdouble::Bool=false)
```

Type representing symmetries. The symmetry is a magnetic space group symmetry, i.e.,
symmetry containing rotation, reflection, inversion, translation and time reversal.
The defining feature for magnetic space group symmetry is its represention in real
space (`rotation_matrix` and `translation`, where translation operation follows the
the action of `rotation_matrix`) and its representation in spinors. The symmetry may
belong to a double group (`isdouble=true`) or not (`isdouble=false`). If the symmetry
does not belong to a double group, it is only relevant for spinless fermions and the
`spin_rep` information is not used at all. The symmetry may
be antilinear denoted by the keyword `isantilinear`.
"""
function Symmetry(rotation_matrix, translation, spin_rep;
    isantilinear::Bool=false, isdouble::Bool=false)
    return Symmetry(rotation_matrix, translation, spin_rep, isantilinear, isdouble)
end

function *(s1::Symmetry, s2::Symmetry)
    if !(s1.isdouble == s2.isdouble)
        error("cannot multiply ordinary elements with double group elements.")
    end
    nrotation_matrix = s1.rotation_matrix*s2.rotation_matrix
    ntranslation = s1.translation+s1.rotation_matrix*s2.translation
    nspin_rep = s1.spin_rep*(s1.isantilinear ? conj(s2.spin_rep) : s2.spin_rep)
    nantilinear = xor(s1.isantilinear, s2.isantilinear)
    return Symmetry(nrotation_matrix, ntranslation, nspin_rep,
        isantilinear = nantilinear, isdouble = s1.isdouble)
end

function ==(s1::Symmetry, s2::Symmetry)
    if !(s1.isantilinear == s2.isantilinear)
        return false
    end
    if !(s1.isdouble == s2.isdouble)
        return false
    end
    if !isapprox(s1.rotation_matrix, s2.rotation_matrix)
        return false
    end
    if !isapprox(s1.translation, s2.translation)
        return false
    end
    if s1.isdouble
        if !isapprox(s1.spin_rep, s2.spin_rep)
            return false
        end
    end
    return true
end

function Base.show(io::IO, s::Symmetry)
    dcpn = decompose_to_primitive(s)
    str = ""
    if dcpn[:TRS] str = "T"*str end
    if abs(dcpn[:θ]) > 0
        rationalized_θ = rationalize(dcpn[:θ]/(2π), tol=√eps())
        str = "C($(dcpn.axis), $(rationalized_θ))"*str
    end
    if dcpn[:inversion] str = "I"*str end
    if sum(abs.(dcpn[:translation])) > 0
        str = "t($(dcpn[:translation]))"*str
    end
    if s.isdouble && dcpn[:Ē] str = "Ē"*str end
    if str == "" str = "E"*str end
    print(io, "Symmetry: $str")
end

"""
```julia
inverse(s::Symmetry)
```

Get inverse symmetry operation for `s`.
"""
function inverse(s::Symmetry)
    rotation_matrix = inv(s.rotation_matrix)
    translation = -rotation_matrix*s.translation
    spin_rep = s.isantilinear ? conj(inv(s.spin_rep)) : inv(s.spin_rep)
    return Symmetry(rotation_matrix, translation, spin_rep, isantilinear=s.isantilinear, isdouble=s.isdouble)
end


"""
```julia
getJ(l::Real)::Vector{Matrix{ComplexF64}}
```

Calculate matrix represention of angular momentum operator [Jx, Jy, Jz] for
spin J particles.

Ref: wikipedia, Spin_(physics)#Higher_spins.
"""
function getJ(l::Real)::Vector{Matrix{ComplexF64}}
    dim = Int64(2l + 1)
    Jx, Jy, Jz = [zeros(ComplexF64, dim, dim) for i in 1:3]
    for n in 1:dim, m in 1:dim
        Jx[m, n] = ((m==n+1)+(m+1==n))/2*sqrt((l+1)*(m+n-1)-m*n)
        Jy[m, n] = im*((m==n+1)-(m+1==n))/2*sqrt((l+1)*(m+n-1)-m*n)
        Jz[m, n] = (l+1-m)*(m==n)
    end
    return [Jx, Jy, Jz]
end


function cross_product_matrix(u)
    return [0 -u[3] u[2]; u[3] 0 -u[1]; -u[2] u[1] 0]
end

function _get_rotation_matrix(θ::Real, axisin::Vector{<:Real})
    axis = normalize(axisin)
    rotation_matrix = cos(θ)*I+sin(θ)*cross_product_matrix(axis)+(1-cos(θ))*(axis*(axis'))
end

"""
```julia
gettranslation(v::Vector{<:Real}; isdouble::Bool=false)
```

get translation along `v`.
"""
function gettranslation(v::Vector{<:Real}; isdouble::Bool=false)
    return Symmetry(I(3), v, I(2), isdouble=isdouble)
end

"""
```julia
getrotation(θ::Real, axis::Vector{<:Real}; isdouble::Bool=false)
```

get `θ` rotation around the `axis`.
"""
function getrotation(θ::Real, axisin::Vector{<:Real}; isdouble::Bool=false)
    axis = normalize(axisin)
    rotation_matrix = _get_rotation_matrix(θ, axis)
    s = getJ(1/2)
    spin_rep = exp(-im*θ*reduce(+, [axis[i]*s[i] for i = 1:3]))
    return Symmetry(rotation_matrix, [0.0, 0.0, 0.0], spin_rep, isdouble=isdouble)
end

"""
```julia
getinversion(; isdouble::Bool=false)
```

get inversion symmetry.
"""
function getinversion(; isdouble::Bool=false)
    return Symmetry(-I(3), [0.0, 0.0, 0.0], I(2), isdouble=isdouble)
end

"""
```julia
getmirror(axis::Vector{<:Real}; isdouble::Bool=false) -> Symmetry

get mirror symmetry along `axis`.
```
"""
function getmirror(axis::Vector{<:Real}; isdouble::Bool=false)
    return getrotation(pi, axis, isdouble=isdouble)*getinversion(isdouble=isdouble)
end

"""
```julia
getTRS(; isdouble::Bool=false)
```

get time reversal symmetry.
"""
function getTRS(; isdouble::Bool=false)
    return Symmetry(I(3), [0.0, 0.0, 0.0], [0 -1; 1 0], isantilinear=true, isdouble=isdouble)
end

"""
```julia
isimproper(s::Symmetry)
```

check whether the rotation part of `s` are improper or not.
"""
isimproper(s::Symmetry) = det(s.rotation_matrix)<0

"""
```julia
decompose_to_primitive(s::Symmetry)
```

Every magnetic space group symmetry can be decomposed into a canonical
representation, i.e., Ē*t*I*C*T, where T is time reversal, C is rotation, I
is inversion, t is translation and Ē is 2π rotation. This function decomposes
the symmetry to these atomic actions. This function always treats the symmetry
as a double group symmetry.

The returned named tuple are in the following format
```julia
(TRS::Bool, inversion::Bool, translation::Vector{Float64}, θ::Float64,
    axis::Vector{Float64}, Ē::Bool)
```
"""
function decompose_to_primitive(s::Symmetry)
    TRS = s.isantilinear
    inversion = isimproper(s)
    proper_rotation_matrix = (-1)^(isimproper(s))*s.rotation_matrix
    axis = normalize(nullspace(proper_rotation_matrix-I, atol=√eps())[:, 1])
    if axis[3] < 0 axis = -axis end
    myacos(x) = abs(abs(x)-1.0)<√eps() ? (x>0 ? 0.0 : π) : acos(x)
    θ = myacos((tr(proper_rotation_matrix)-1)/2)
    if !(proper_rotation_matrix ≈ _get_rotation_matrix(θ, axis))
        @assert proper_rotation_matrix ≈ _get_rotation_matrix(-θ, axis)
        θ = -θ
    end
    tmp = getrotation(θ, axis, isdouble=s.isdouble).spin_rep*(TRS ? [0 -1; 1 0] : I(2))
    if isapprox(s.spin_rep, tmp)
        Ē = false
    else
        @assert isapprox(s.spin_rep, -tmp)
        Ē = true
    end
    return (TRS=TRS, inversion=inversion, translation=deepcopy(s.translation), θ=θ, axis=axis, Ē=Ē)
end

"""
```julia
get_orb_rep(s::Symmetry, l::Integer)::Matrix{ComplexF64}
```

get representation of `s` in `l` type orbitals. `l` is an integer representing
the eigenvalues of angular momentum.
"""
function get_orb_rep(s::Symmetry, l::Integer)::Matrix{ComplexF64}
    dcpn = decompose_to_primitive(s)
    L = getJ(l)
    # for rotation
    rep = exp(-im*dcpn[:θ]*reduce(+, [dcpn[:axis][i]*L[i] for i=1:3]))
    # for inversion
    if dcpn[:inversion]
        rep = (-1)^l*rep
    end
    # for TRS
    if dcpn[:TRS]
        TRSrep = zeros(ComplexF64, 2l+1, 2l+1)
        for m = -l:l
            TRSrep[l+m+1, l-m+1] = (-1)^m
        end
        rep = rep*TRSrep
    end
    if s.isdouble
        return kron(s.spin_rep, rep)
    else
        return rep
    end
end


function _populate_matrix!(mat, xind, yind, v; isspinful=false)
    if isspinful
        matdim = size(mat, 1)÷2
        vdim = size(v, 1)÷2
        for m in [0, 1], n in [0, 1]
            mat[xind.+m*matdim, yind.+n*matdim] = v[(1:vdim).+m*vdim, (1:vdim).+n*vdim]
        end
    else
        mat[xind, yind] = v
    end
end

"""
```julia
get_orb_rep(s::Symmetry, ls::Vector{Int64})::Matrix{ComplexF64}
```

get representation of `s` in `ls` type orbitals. `ls` is a list of integers representing
the eigenvalues of angular momentum.
"""
@memoize function get_orb_rep(s::Symmetry, ls::Vector{<:Integer})::Matrix{ComplexF64}
    dim = sum((2ls.+1))
    rep = zeros(ComplexF64, dim*(s.isdouble+1), dim*(s.isdouble+1))
    cnt = 0
    for l in ls
        rep_l = get_orb_rep(s, l)
        _populate_matrix!(rep, (cnt+1):(cnt+2l+1), (cnt+1):(cnt+2l+1), rep_l, isspinful=s.isdouble)
        cnt += 2l+1
    end
    return rep
end

function equal_mod_translation(s1::Symmetry, s2::Symmetry)
    if s1.rotation_matrix≈s2.rotation_matrix && s1.isantilinear==s2.isantilinear && s1.isdouble==s2.isdouble
        if s1.isdouble
            if s1.spin_rep ≈ s2.spin_rep
                return true
            else
                return false
            end
        end
        return true
    end
    return false
end

function in_mod_translation(s, symmetries::Vector{Symmetry})
    for s′ in symmetries
        if equal_mod_translation(s, s′)
            return true
        end
    end
    return false
end

"""
```julia
generategroup(symmetries::Vector{Symmetry}; maxiter=100)::Vector{Symmetry}
```

generate a group with generators `symmetries`.
"""
function generategroup(symmetries::Vector{Symmetry}; maxiter=100)::Vector{Symmetry}
    group = deepcopy(symmetries)
    n = 0
    while true
        lastgroup = deepcopy(group)
        for s1 in symmetries, s2 in lastgroup
            for ns in [s1*s2, s2*s1]
                if !in_mod_translation(ns, group)
                    push!(group, ns)
                end
            end
        end
        if length(lastgroup) == length(group)
            return group
        end
        n += 1
        if n > maxiter
            error("Maximal iteration number reached. Is this a finite group?")
        end
    end
end

function _get_transformed_site(
    s::Symmetry,
    lat::AbstractMatrix{Float64},
    site_positions::Matrix{Float64},
    R::AbstractVector{<:Integer},
    i::Integer;
    position_tolerance::Float64=1.0e-2
)::Tuple{Vector{Int64},Int64}
    pos = lat*R+site_positions[:, i]
    npos = s.rotation_matrix*pos+s.translation
    nsites = size(site_positions, 2)
    ni = 0
    nR = [0, 0, 0]
    for cnt in 1:nsites
        tmp = inv(lat)*(npos-site_positions[:, cnt])
        nR = round.(tmp)
        if sum(abs.(tmp - nR)) < position_tolerance
            ni = cnt
            break
        end
    end
    if ni == 0
        error("cannot find transformed site for site $i !")
    end
    return (nR, ni)
end

function _get_transformed_site(
    s::Symmetry,
    tm::TBModel,
    R::AbstractVector{<:Integer},
    i::Integer;
    position_tolerance::Float64=1.0e-2
)::Tuple{Vector{Int64},Int64}
    return _get_transformed_site(s, tm.lat, tm.site_positions, R, i, position_tolerance=position_tolerance)
end

function transform_hamiltonian!(
    ntm::TBModel,
    s::Symmetry,
    tm::TBModel;
    position_tolerance::Float64=1.0e-2
)
    empty!(ntm.hoppings)
    for (R, hopping) in tm.hoppings, i in 1:tm.nsites, j in 1:tm.nsites
        indi = HopTB._to_orbital_index(tm, i)
        indj = HopTB._to_orbital_index(tm, j)
        site_hopping = hopping[indi, indj]
        nR1, ni = _get_transformed_site(inverse(s), tm, [0, 0, 0], i; position_tolerance=position_tolerance)
        nR2, nj = _get_transformed_site(inverse(s), tm, R, j; position_tolerance=position_tolerance)
        ui = get_orb_rep(s, tm.orbital_types[i])
        uj = get_orb_rep(s, tm.orbital_types[j])
        nindi = HopTB._to_orbital_index(tm, ni)
        nindj = HopTB._to_orbital_index(tm, nj)
        if !(nR2-nR1 in keys(ntm.hoppings))
            ntm.hoppings[nR2 - nR1] = zeros(ComplexF64, ntm.norbits, ntm.norbits)
        end
        ntm.hoppings[nR2 - nR1][nindi, nindj] =
            s.isantilinear ? conj(ui' * site_hopping * uj) : ui' * site_hopping * uj
    end
    return ntm
end

function transform_overlap!(
    ntm::TBModel,
    s::Symmetry,
    tm::TBModel;
    position_tolerance::Float64=1.0e-2
)
    empty!(ntm.overlaps)
    for (R, overlap) in tm.overlaps, i in 1:tm.nsites, j in 1:tm.nsites
        indi = HopTB._to_orbital_index(tm, i)
        indj = HopTB._to_orbital_index(tm, j)
        site_overlap = overlap[indi, indj]
        nR1, ni = _get_transformed_site(inverse(s), tm, [0, 0, 0], i; position_tolerance=position_tolerance)
        nR2, nj = _get_transformed_site(inverse(s), tm, R, j; position_tolerance=position_tolerance)
        ui = get_orb_rep(s, tm.orbital_types[i])
        uj = get_orb_rep(s, tm.orbital_types[j])
        nindi = HopTB._to_orbital_index(tm, ni)
        nindj = HopTB._to_orbital_index(tm, nj)
        if !(nR2-nR1 in keys(ntm.overlaps))
            ntm.overlaps[nR2 - nR1] = zeros(ComplexF64, ntm.norbits, ntm.norbits)
        end
        ntm.overlaps[nR2 - nR1][nindi, nindj] =
            s.isantilinear ? conj(ui' * site_overlap * uj) : ui' * site_overlap * uj
    end
    return ntm
end

"""
Transform position matrix of tm and store the new matrix into ntm. The original position
matrix of ntm will be discarded.

Notice that ntm should have a transformed overlap matrix of tm.
"""
function transform_position!(
    ntm::TBModel{T},
    s::Symmetry,
    tm::TBModel{T};
    position_tolerance::Float64=1.0e-2
) where T <: Number
    empty!(ntm.positions)
    for (R, pos) in tm.positions, i in 1:tm.nsites, j in 1:tm.nsites
        indi = HopTB._to_orbital_index(tm, i)
        indj = HopTB._to_orbital_index(tm, j)
        site_pos = [pos[α][indi, indj] for α in 1:3]
        nR1, ni = _get_transformed_site(inverse(s), tm, [0, 0, 0], i; position_tolerance=position_tolerance)
        nR2, nj = _get_transformed_site(inverse(s), tm, R, j; position_tolerance=position_tolerance)
        ui = get_orb_rep(s, tm.orbital_types[i])
        uj = get_orb_rep(s, tm.orbital_types[j])
        nindi = HopTB._to_orbital_index(tm, ni)
        nindj = HopTB._to_orbital_index(tm, nj)
        if !(nR2 - nR1 in keys(ntm.positions))
            ntm.positions[nR2 - nR1] = [zeros(ComplexF64, ntm.norbits, ntm.norbits) for α in 1:3]
        end
        rminv = inv(s.rotation_matrix)
        tmp1 = site_pos - [
            s.translation[α] * (
                R in keys(tm.overlaps) ? tm.overlaps[R][indi, indj] : zeros(T, length(indi), length(indj))
            ) for α in 1:3
        ]
        tmp2 = rminv * [s.isantilinear ? conj(ui' * tmp1[α] * uj) : ui' * tmp1[α] * uj for α in 1:3]
        for α in 1:3
            ntm.positions[nR2 - nR1][α][nindi, nindj] = tmp2[α] - (tm.lat * nR1)[α] * (
                nR1 - nR2 in keys(ntm.overlaps) ? ntm.overlaps[nR2 - nR1][nindi, nindj] : zeros(T, length(nindi), length(nindj))
            )
        end
    end
    return ntm
end

function transform(tm::TBModel, s::Symmetry; position_tolerance::Float64=1.0e-2)
    ntm = deepcopy(tm)
    transform_hamiltonian!(ntm, s, tm; position_tolerance=position_tolerance)
    transform_overlap!(ntm, s, tm; position_tolerance=position_tolerance)
    transform_position!(ntm, s, tm; position_tolerance=position_tolerance)
    return ntm
end

"""
```julia
symmetrize(tm::TBModel, g::Vector{Symmetry})
```

Symmetrize `tm` by the symmetry group `g`.
"""
function symmetrize(tm::TBModel, g::Vector{Symmetry}; position_tolerance::Float64=1.0e-2)
    @assert tm.is_canonical_ordered
    ns = length(g)
    ntm = deepcopy(tm)
    empty!(ntm.hoppings)
    empty!(ntm.overlaps)
    empty!(ntm.positions)
    for s in g
        tm_tmp = transform(tm, s; position_tolerance=position_tolerance)
        for (k, hopping) in tm_tmp.hoppings
            ntm.hoppings[k] = get(ntm.hoppings, k, zeros(ComplexF64, size(hopping))) + hopping / ns
        end
        for (k, overlap) in tm_tmp.overlaps
            ntm.overlaps[k] = get(ntm.overlaps, k, zeros(ComplexF64, size(overlap))) + overlap / ns
        end
        for (k, pos) in tm_tmp.positions
            ntm.positions[k] = get(ntm.positions, k, [zeros(ComplexF64, size(pos[1])) for α in 1:3]) + pos / ns
        end
    end
    return ntm
end

function get_k_rep(s::Symmetry)
    return (-1)^(s.isantilinear)*s.rotation_matrix
end

"""
```julia
get_bloch_rep(s::Symmetry, tm::TBModel, orbital_types::Vector{Vector{Int64}},
    k::Vector{<:Real}; isspinful::Bool=false)
```
"""
function get_bloch_rep(s::Symmetry, tm::TBModel, k::Vector{<:Real})
    rep = zeros(ComplexF64, tm.norbits, tm.norbits)
    kc = tm.rlat*k
    for i in 1:tm.nsites
        R, j = _get_transformed_site(s, tm.lat, tm.site_positions, [0, 0, 0], i)
        rep[HopTB._to_orbital_index(tm, j), HopTB._to_orbital_index(tm, i)] =
            get_orb_rep(s, tm.orbital_types[i])*exp(-im*((get_k_rep(s)*kc)⋅(tm.lat*R)))
    end
    return rep
end


# This transformation is used to change the basis from openmx default to canonical basis.
# It can be directly put into HopTB.changebasis
const Us_openmx = Dict(
    0 => Matrix{ComplexF64}(I, 1, 1),
    1 => [-1/√2 im/√2 0; 0 0 1; 1/√2 im/√2 0],
    2 => [0 1/√2 -im/√2 0 0; 0 0 0 -1/√2 im/√2; 1 0 0 0 0; 0 0 0 1/√2 im/√2; 0 1/√2 im/√2 0 0],
    3 => [0 0 0 0 0 -1/√2 im/√2; 0 0 0 1/√2 -im/√2 0 0; 0 -1/√2 im/√2 0 0 0 0; 1 0 0 0 0 0 0;
        0 1/√2 im/√2 0 0 0 0; 0 0 0 1/√2 im/√2 0 0; 0 0 0 0 0 1/√2 im/√2]
);

const Us_openmx_inv = Dict{Int64,Matrix{ComplexF64}}()

for l in keys(Us_openmx)
    Us_openmx_inv[l] = Us_openmx[l]'
end

# This transformation is used to change the basis from wannier90 default to canonical basis.
# It can be directly put into HopTB.changebasis
const Us_wann = Dict(
    0 => Matrix{ComplexF64}(I, 1, 1),
    1 => [0 -1/√2 im/√2; 1 0 0; 0 1/√2 im/√2],
    2 => [0 0 0 1/√2 -im/√2; 0 -1/√2 im/√2 0 0; 1 0 0 0 0; 0 1/√2 im/√2 0 0; 0 0 0 1/√2 im/√2],
    3 => [0 0 0 0 0 -1/√2 im/√2; 0 0 0 1/√2 -im/√2 0 0; 0 -1/√2 im/√2 0 0 0 0; 1 0 0 0 0 0 0;
        0 1/√2 im/√2 0 0 0 0; 0 0 0 1/√2 im/√2 0 0; 0 0 0 0 0 1/√2 im/√2]
)

end
