module Interface

using StaticArrays, LinearAlgebra
using ..HopTB

export createmodelaims, createmodelopenmx, createmodelwannier

# FHI-aims
"""
create TBModel from FHI-aims interface.
"""
function createmodelaims(filepath::String)
    f = open(filepath)
    # number of basis
    @assert occursin("n_basis", readline(f)) # start
    norbits = parse(Int64, readline(f))
    @assert occursin("end", readline(f)) # end
    @assert occursin("n_ham", readline(f)) # start
    nhams = parse(Int64, readline(f))
    @assert occursin("end", readline(f)) # end
    @assert occursin("n_cell", readline(f)) # start
    ncells = parse(Int64, readline(f))
    @assert occursin("end", readline(f)) # end
    # lattice vector
    @assert occursin("lattice_vector", readline(f)) # start
    lat = Matrix{Float64}(I, 3, 3)
    for i in 1:3
        lat[:, i] = map(x->parse(Float64, x), split(readline(f)))
    end
    @assert occursin("end", readline(f)) # end
    # hamiltonian
    @assert occursin("hamiltonian", readline(f)) # start
    hamiltonian = zeros(nhams)
    i = 1
    while true
        @assert !eof(f)
        ln = split(readline(f))
        if occursin("end", ln[1]) break end
        hamiltonian[i:i + length(ln) - 1] = map(x->parse(Float64, x), ln)
        i += length(ln)
    end
    # overlaps
    @assert occursin("overlap", readline(f)) # start
    overlaps = zeros(nhams)
    i = 1
    while true
        @assert !eof(f)
        ln = split(readline(f))
        if occursin("end", ln[1]) break end
        overlaps[i:i + length(ln) - 1] = map(x->parse(Float64, x), ln)
        i += length(ln)
    end
    # index hamiltonian
    @assert occursin("index_hamiltonian", readline(f)) # start
    indexhamiltonian = zeros(Int64, ncells * norbits, 4)
    i = 1
    while true
        @assert !eof(f)
        ln = split(readline(f))
        if occursin("end", ln[1]) break end
        indexhamiltonian[i, :] = map(x->parse(Int64, x), ln)
        i += 1
    end
    # cell index
    @assert occursin("cell_index", readline(f)) # start
    cellindex = zeros(Int64, ncells, 3)
    i = 1
    while true
        @assert !eof(f)
        ln = split(readline(f))
        if occursin("end", ln[1]) break end
        if i <= ncells
            cellindex[i, :] = map(x->parse(Int64, x), ln)
        end
        i += 1
    end
    # column index hamiltonian
    @assert occursin("column_index_hamiltonian", readline(f)) # start
    columnindexhamiltonian = zeros(Int64, nhams)
    i = 1
    while true
        @assert !eof(f)
        ln = split(readline(f))
        if occursin("end", ln[1]) break end
        columnindexhamiltonian[i:i + length(ln) - 1] = map(x->parse(Int64, x), ln)
        i += length(ln)
    end
    # positions
    positions = zeros(nhams, 3)
    for dir in 1:3
        positionsdir = zeros(nhams)
        @assert occursin("position", readline(f)) # start
        readline(f) # skip direction
        i = 1
        while true
            @assert !eof(f)
            ln = split(readline(f))
            if occursin("end", ln[1]) break end
            positionsdir[i:i + length(ln) - 1] = map(x->parse(Float64, x), ln)
            i += length(ln)
        end
        positions[:, dir] = positionsdir
    end
    if !eof(f)
        withsoc = true
        soc_matrix = zeros(nhams, 3)
        for dir in 1:3
            socdir = zeros(nhams)
            @assert occursin("soc_matrix", readline(f)) # start
            readline(f) # skip direction
            i = 1
            while true
                @assert !eof(f)
                ln = split(readline(f))
                if occursin("end", ln[1]) break end
                socdir[i:i + length(ln) - 1] = map(x->parse(Float64, x), ln)
                i += length(ln)
            end
            soc_matrix[:, dir] = socdir
        end
    else
        withsoc = false
    end
    close(f)

    if withsoc
        σx = [0 1; 1 0]
        σy = [0 -im; im 0]
        σz = [1 0; 0 -1]
        σ0 = [1 0; 0 1]
        nm = TBModel{ComplexF64}(2*norbits, lat, isorthogonal=false)
        # convention here is first half up (spin=0); second half down (spin=1).
        for i in 1:size(indexhamiltonian, 1)
            for j in indexhamiltonian[i, 3]:indexhamiltonian[i, 4]
                for nspin in 0:1
                    for mspin in 0:1
                        sethopping!(nm,
                            cellindex[indexhamiltonian[i, 1], :],
                            columnindexhamiltonian[j] + norbits * nspin,
                            indexhamiltonian[i, 2] + norbits * mspin,
                            σ0[nspin + 1, mspin + 1] * hamiltonian[j] -
                            (σx[nspin + 1, mspin + 1] * soc_matrix[j, 1] +
                            σy[nspin + 1, mspin + 1] * soc_matrix[j, 2] +
                            σz[nspin + 1, mspin + 1] * soc_matrix[j, 3]) * im)
                        setoverlap!(nm,
                            cellindex[indexhamiltonian[i, 1], :],
                            columnindexhamiltonian[j] + norbits * nspin,
                            indexhamiltonian[i, 2] + norbits * mspin,
                            σ0[nspin + 1, mspin + 1] * overlaps[j])
                    end
                end
            end
        end
        for i in 1:size(indexhamiltonian, 1)
            for j in indexhamiltonian[i, 3]:indexhamiltonian[i, 4]
                for nspin in 0:1
                    for mspin in 0:1
                        for dir in 1:3
                            setposition!(nm,
                                cellindex[indexhamiltonian[i, 1], :],
                                columnindexhamiltonian[j] + norbits * nspin,
                                indexhamiltonian[i, 2] + norbits * mspin,
                                dir,
                                σ0[nspin + 1, mspin + 1] * positions[j, dir])
                        end
                    end
                end
            end
        end
        return nm
    else
        nm = TBModel{Float64}(norbits, lat, isorthogonal=false)
        for i in 1:size(indexhamiltonian, 1)
            for j in indexhamiltonian[i, 3]:indexhamiltonian[i, 4]
                sethopping!(nm,
                    cellindex[indexhamiltonian[i, 1], :],
                    columnindexhamiltonian[j],
                    indexhamiltonian[i, 2],
                    hamiltonian[j])
                setoverlap!(nm,
                    cellindex[indexhamiltonian[i, 1], :],
                    columnindexhamiltonian[j],
                    indexhamiltonian[i, 2],
                    overlaps[j])
            end
        end
        for i in 1:size(indexhamiltonian, 1)
            for j in indexhamiltonian[i, 3]:indexhamiltonian[i, 4]
                for dir in 1:3
                    setposition!(nm,
                        cellindex[indexhamiltonian[i, 1], :],
                        columnindexhamiltonian[j],
                        indexhamiltonian[i, 2],
                        dir,
                        positions[j, dir])
                end
            end
        end
        return nm
    end
end


function _parseopenmx(filepath::String)
    # define some helper functions for mixed structure of OpenMX binary data file.
    function multiread(::Type{T}, f, size)::Vector{T} where T
        ret = Vector{T}(undef, size)
        read!(f, ret);ret
    end
    multiread(f, size) = multiread(Int32, f, size)

    function read_mixed_matrix(::Type{T}, f, dims::Vector{<:Integer}) where T
        ret::Vector{Vector{T}} = []
        for i = dims; t = Vector{T}(undef, i);read!(f, t);push!(ret, t); end; ret
    end

    function read_matrix_in_mixed_matrix(::Type{T}, f, spins, atomnum, FNAN, natn, Total_NumOrbs) where T
        ret = Vector{Vector{Vector{Matrix{T}}}}(undef, spins)
        for spin = 1:spins;t_spin = Vector{Vector{Matrix{T}}}(undef, atomnum)
            for ai = 1:atomnum;t_ai = Vector{Matrix{T}}(undef, FNAN[ai])
                for aj_inner = 1:FNAN[ai]
                    t = Matrix{T}(undef, Total_NumOrbs[natn[ai][aj_inner]], Total_NumOrbs[ai])
                    read!(f, t);t_ai[aj_inner] = t
                end;t_spin[ai] = t_ai
            end;ret[spin] = t_spin
        end;return ret
    end
    read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs) = read_matrix_in_mixed_matrix(Float64, f, spins, atomnum, FNAN, natn, Total_NumOrbs)

    read_3d_vecs(::Type{T}, f, num) where T = reshape(multiread(T, f, 4 * num), 4, Int(num))[2:4,:]
    read_3d_vecs(f, num) = read_3d_vecs(Float64, f, num)
    # End of helper functions

    bound_multiread(T, size) = multiread(T, f, size)
    bound_multiread(size) = multiread(f, size)
    bound_read_mixed_matrix() = read_mixed_matrix(Int32, f, FNAN)
    bound_read_matrix_in_mixed_matrix(spins) = read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs)
    bound_read_3d_vecs(num) = read_3d_vecs(f, num)
    bound_read_3d_vecs(::Type{T}, num) where T = read_3d_vecs(T, f, num)
    # End of bound helper functions

    f = open(filepath)
    atomnum, SpinP_switch, Catomnum, Latomnum, Ratomnum, TCpyCell, order_max = bound_multiread(7)
    @assert (SpinP_switch >> 2) == 3
    SpinP_switch &= 0x03

    atv, atv_ijk = bound_read_3d_vecs.([Float64,Int32], TCpyCell + 1)

    Total_NumOrbs, FNAN = bound_multiread.([atomnum,atomnum])
    FNAN .+= 1
    natn = bound_read_mixed_matrix()
    ncn = ((x)->x .+ 1).(bound_read_mixed_matrix()) # These is to fix that atv and atv_ijk starts from 0 in original C code.

    tv, rtv, Gxyz = bound_read_3d_vecs.([3,3,atomnum])

    Hk = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)
    iHk = SpinP_switch == 3 ? bound_read_matrix_in_mixed_matrix(3) : nothing
    OLP = bound_read_matrix_in_mixed_matrix(1)[1]
    OLP_r = []
    for dir in 1:3, order in 1:order_max
        t = bound_read_matrix_in_mixed_matrix(1)[1]
        if order == 1 push!(OLP_r, t) end
    end
    OLP_p = bound_read_matrix_in_mixed_matrix(3)
    DM = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)
    iDM = bound_read_matrix_in_mixed_matrix(2)
    solver = bound_multiread(1)[1]
    chem_p, E_temp = bound_multiread(Float64, 2)
    dipole_moment_core, dipole_moment_background = bound_multiread.(Float64, [3,3])
    Valence_Electrons, Total_SpinS = bound_multiread(Float64, 2)
    dummy_blocks = bound_multiread(1)[1]
    for i in 1:dummy_blocks
        bound_multiread(UInt8, 256)
    end

    # we suppose that the original output file(.out) was appended to the end of the scfout file.
    function strip1(s::Vector{UInt8})
        startpos = 0
        for i = 1:length(s) + 1
            if i > length(s) || s[i] & 0x80 != 0 || !isspace(Char(s[i] & 0x7f))
                startpos = i
                break
            end
        end
        return s[startpos:end]
    end
    function startswith1(s::Vector{UInt8}, prefix::Vector{UInt8})
        return length(s) >= length(prefix) && s[1:length(prefix)] == prefix
    end
    target_line = Vector{UInt8}("Fractional coordinates of the final structure")
    while !startswith1(strip1(Vector{UInt8}(readline(f))), target_line)
        if eof(f)
            error("Atom positions not found. Please check if the .out file was appended to the end of .scfout file!")
        end
    end
    for i = 1:2;@assert readline(f) == "***********************************************************";end
    @assert readline(f) == ""
    atom_frac_pos = zeros(3, atomnum)
    for i = 1:atomnum
        m = match(r"^\s*\d+\s+\w+\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)", readline(f))
        atom_frac_pos[:,i] = ((x)->parse(Float64, x)).(m.captures)
    end
    atom_pos = tv * atom_frac_pos
    close(f)

    # use the atom_pos to fix
    # TODO: Persuade wangc to accept the following code, which seems hopeless and meaningless.
    """
    for axis = 1:3
        ((x2, y2, z)->((x, y)->x .+= z * y).(x2, y2)).(OLP_r[axis], OLP, atom_pos[axis,:])
    end
    """
    for axis in 1:3,i in 1:atomnum, j in 1:FNAN[i]
        OLP_r[axis][i][j] .+= atom_pos[axis,i] * OLP[i][j]
    end

    # fix type mismatch
    atv_ijk = Matrix{Int16}(atv_ijk)

    return atomnum, SpinP_switch, atv, atv_ijk, Total_NumOrbs, FNAN, natn, ncn, tv, Hk, iHk, OLP, OLP_r
end

function _parseopenmx38(filepath::String)
    # define some helper functions for mixed structure of OpenMX binary data file.
    function multiread(::Type{T}, f, size)::Vector{T} where T
        ret = Vector{T}(undef, size)
        read!(f, ret);ret
    end
    multiread(f, size) = multiread(Int32, f, size)

    function read_mixed_matrix(::Type{T}, f, dims::Vector{<:Integer}) where T
        ret::Vector{Vector{T}} = []
        for i = dims; t = Vector{T}(undef, i);read!(f, t);push!(ret, t); end; ret
    end

    function read_matrix_in_mixed_matrix(::Type{T}, f, spins, atomnum, FNAN, natn, Total_NumOrbs) where T
        ret = Vector{Vector{Vector{Matrix{T}}}}(undef, spins)
        for spin = 1:spins;t_spin = Vector{Vector{Matrix{T}}}(undef, atomnum)
            for ai = 1:atomnum;t_ai = Vector{Matrix{T}}(undef, FNAN[ai])
                for aj_inner = 1:FNAN[ai]
                    t = Matrix{T}(undef, Total_NumOrbs[natn[ai][aj_inner]], Total_NumOrbs[ai])
                    read!(f, t);t_ai[aj_inner] = t
                end;t_spin[ai] = t_ai
            end;ret[spin] = t_spin
        end;return ret
    end
    read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs) = read_matrix_in_mixed_matrix(Float64, f, spins, atomnum, FNAN, natn, Total_NumOrbs)

    read_3d_vecs(::Type{T}, f, num) where T = reshape(multiread(T, f, 4 * num), 4, Int(num))[2:4,:]
    read_3d_vecs(f, num) = read_3d_vecs(Float64, f, num)
    # End of helper functions

    bound_multiread(T, size) = multiread(T, f, size)
    bound_multiread(size) = multiread(f, size)
    bound_read_mixed_matrix() = read_mixed_matrix(Int32, f, FNAN)
    bound_read_matrix_in_mixed_matrix(spins) = read_matrix_in_mixed_matrix(f, spins, atomnum, FNAN, natn, Total_NumOrbs)
    bound_read_3d_vecs(num) = read_3d_vecs(f, num)
    bound_read_3d_vecs(::Type{T}, num) where T = read_3d_vecs(T, f, num)
    # End of bound helper functions

    f = open(filepath)
    atomnum, SpinP_switch, Catomnum, Latomnum, Ratomnum, TCpyCell = bound_multiread(6)

    atv, atv_ijk = bound_read_3d_vecs.([Float64,Int32], TCpyCell + 1)

    Total_NumOrbs, FNAN = bound_multiread.([atomnum,atomnum])
    FNAN .+= 1
    natn = bound_read_mixed_matrix()
    ncn = ((x)->x .+ 1).(bound_read_mixed_matrix()) # These is to fix that atv and atv_ijk starts from 0 in original C code.

    tv, rtv, Gxyz = bound_read_3d_vecs.([3,3,atomnum])

    Hk = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)
    iHk = SpinP_switch == 3 ? bound_read_matrix_in_mixed_matrix(3) : nothing
    OLP, OLP_rx, OLP_ry, OLP_rz = bound_read_matrix_in_mixed_matrix(4)
    OLP_r = [OLP_rx,OLP_ry,OLP_rz]
    DM = bound_read_matrix_in_mixed_matrix(SpinP_switch + 1)

    # we suppose that the original output file(.out) was appended to the end of the scfout file.
    function strip1(s::Vector{UInt8})
        startpos = 0
        for i = 1:length(s) + 1
            if i > length(s) || s[i] & 0x80 != 0 || !isspace(Char(s[i] & 0x7f))
                startpos = i
                break
            end
        end
        return s[startpos:end]
    end
    function startswith1(s::Vector{UInt8}, prefix::Vector{UInt8})
        return length(s) >= length(prefix) && s[1:length(prefix)] == prefix
    end
    target_line = Vector{UInt8}("Fractional coordinates of the final structure")
    while !startswith1(strip1(Vector{UInt8}(readline(f))), target_line)
        if eof(f)
            error("Atom positions not found. Please check if the .out file was appended to the end of .scfout file!")
        end
    end
    for i = 1:2;@assert readline(f) == "***********************************************************";end
    @assert readline(f) == ""
    atom_frac_pos = zeros(3, atomnum)
    for i = 1:atomnum
        m = match(r"^\s*\d+\s+\w+\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)\s+([0-9+-.Ee]+)", readline(f))
        atom_frac_pos[:,i] = ((x)->parse(Float64, x)).(m.captures)
    end
    atom_pos = tv * atom_frac_pos
    close(f)

    # use the atom_pos to fix
    # TODO: Persuade wangc to accept the following code, which seems hopeless and meaningless.
    """
    for axis = 1:3
        ((x2, y2, z)->((x, y)->x .+= z * y).(x2, y2)).(OLP_r[axis], OLP, atom_pos[axis,:])
    end
    """
    for axis in 1:3,i in 1:atomnum, j in 1:FNAN[i]
        OLP_r[axis][i][j] .+= atom_pos[axis,i] * OLP[i][j]
    end

    # fix type mismatch
    atv_ijk = Matrix{Int64}(atv_ijk)

    return atomnum, SpinP_switch, atv, atv_ijk, Total_NumOrbs, FNAN, natn, ncn, tv, Hk, iHk, OLP, OLP_r
end

function _createmodelopenmx_inner(filepath::String, parserfunc::Function)
    function calcassistvars(Total_NumOrbs)
        # generate accumulated-indices
        numorb_base = Vector{Int32}(undef, length(Total_NumOrbs))
        numorb_base[1] = 0
        for i = 2:length(Total_NumOrbs)
            numorb_base[i] = numorb_base[i - 1] + Total_NumOrbs[i - 1]
        end
        return numorb_base
    end

    atomnum, SpinP_switch, atv, atv_ijk, Total_NumOrbs, FNAN, natn, ncn, tv, Hk, iHk, OLP, OLP_r = parserfunc(filepath)
    numorb_base = calcassistvars(Total_NumOrbs)
    Total_NumOrbs_sum = sum(Total_NumOrbs)
    ((x)->x .*= 0.529177249).([atv, tv]) # Bohr to Ang
    atv = nothing # atv is never used actually
    for t in [Hk,iHk]
        if !isnothing(t)
            ((x)->((y)->((z)->z .*= 27.211399).(y)).(x)).(t) # Hartree to eV
        end
    end
    ((x)->((y)->((z)->z .*= 0.529177249).(y)).(x)).(OLP_r)


    # sethopping_Hatree(R,i,j,E)=sethopping!(nm,R,i,j,E*27.211399)
    # setposition_Bohr(R,i,j,alpha,r)=setposition!(nm,R,i,j,alpha,r*0.529177249)

    if SpinP_switch == 0
        nm = TBModel{Float64}(Total_NumOrbs_sum, tv, isorthogonal = false)
        for i in 1:atomnum,j in 1:FNAN[i],ii in 1:Total_NumOrbs[i],jj in 1:Total_NumOrbs[natn[i][j]]
            sethopping!(nm, atv_ijk[:,ncn[i][j]], numorb_base[i] + ii, numorb_base[natn[i][j]] + jj, Hk[1][i][j][jj,ii])
            setoverlap!(nm, atv_ijk[:,ncn[i][j]], numorb_base[i] + ii, numorb_base[natn[i][j]] + jj, OLP[i][j][jj,ii])
        end
        for i in 1:atomnum,j in 1:FNAN[i],ii in 1:Total_NumOrbs[i],jj in 1:Total_NumOrbs[natn[i][j]]
            for alpha = 1:3
                setposition!(nm, atv_ijk[:, ncn[i][j]], numorb_base[i] + ii, numorb_base[natn[i][j]] + jj, alpha, OLP_r[alpha][i][j][jj,ii])
            end
        end
    elseif SpinP_switch == 1
        error("Collinear spin is not supported currently")
    elseif SpinP_switch == 3
        for i in 1:length(Hk[4]),j in 1:length(Hk[4][i])
            Hk[4][i][j] += iHk[3][i][j]
            iHk[3][i][j] = -Hk[4][i][j]
        end
        nm = TBModel{ComplexF64}(Total_NumOrbs_sum * 2, tv, isorthogonal = false)
        for spini in 0:1,spinj in (parserfunc === _parseopenmx ? spini : 0):1
            Hk_real, Hk_imag = spini == 0 ? spinj == 0 ? (Hk[1], iHk[1]) : (Hk[3], Hk[4]) : spinj == 0 ? (Hk[3], iHk[3]) : (Hk[2], iHk[2])
            for i in 1:atomnum,j in 1:FNAN[i],ii in 1:Total_NumOrbs[i],jj in 1:Total_NumOrbs[natn[i][j]]
                sethopping!(nm, atv_ijk[:,ncn[i][j]],
                            spini * Total_NumOrbs_sum + numorb_base[i] + ii,
                            spinj * Total_NumOrbs_sum + numorb_base[natn[i][j]] + jj,
                            Hk_real[i][j][jj,ii] + im * Hk_imag[i][j][jj,ii])
                if spini == spinj
                    setoverlap!(nm, atv_ijk[:,ncn[i][j]],
                            spini * Total_NumOrbs_sum + numorb_base[i] + ii,
                            spinj * Total_NumOrbs_sum + numorb_base[natn[i][j]] + jj,
                            OLP[i][j][jj,ii])
                end
            end
        end
        for spini in 0:1,spinj in spini,i in 1:atomnum,j in 1:FNAN[i],ii in 1:Total_NumOrbs[i],jj in 1:Total_NumOrbs[natn[i][j]]
            for alpha = 1:3
                setposition!(nm, atv_ijk[:, ncn[i][j]],
                            spini * Total_NumOrbs_sum + numorb_base[i] + ii,
                            spinj * Total_NumOrbs_sum + numorb_base[natn[i][j]] + jj,
                            alpha, OLP_r[alpha][i][j][jj,ii])
            end
        end
    else
        error("SpinP_switch is $SpinP_switch, rather than valid values 0, 1 or 3")
    end

    return nm
end
function createmodelopenmx(filepath::String)
    return _createmodelopenmx_inner(filepath, _parseopenmx)
end
function createmodelopenmx38(filepath::String)
    return _createmodelopenmx_inner(filepath, _parseopenmx38)
end


# Wannier90
"""
create TBModel from Wannier90 interface.
"""
function createmodelwannier(filepath::String)
    f = open(filepath)
    readline(f) # This line is comment
    lat = zeros(3, 3)
    for i in 1:3
        lat[:, i] = map(s->parse(Float64, s), split(readline(f)))
    end
    norbits = parse(Int64, readline(f))
    nrpts = parse(Int64, readline(f))
    rndegen = zeros(0)
    while true
        line = readline(f)
        if line == "" break end
        rndegen = [rndegen; map(s->parse(Int64, s), split(line))]
    end
    @assert length(rndegen) == nrpts

    om = TBModel{ComplexF64}(norbits, lat, isorthogonal=true)

    for irpt in 1:nrpts
        R = map(s->parse(Int64, s), split(readline(f)))
        for m in 1:norbits
            for n in 1:norbits
                line = readline(f)
                tmp = map(s->parse(Float64, s), split(line)[end - 1:end])
                sethopping!(om, R, n, m, (tmp[1] + im * tmp[2]) / rndegen[irpt])
            end
        end
        @assert readline(f) == ""
    end

    for irpt in 1:nrpts
        R = map(s->parse(Int64, s), split(readline(f)))
        for m in 1:norbits
            for n in 1:norbits
                line = readline(f)
                tmp = map(s->parse(Float64, s), split(line)[end - 5:end])
                setposition!(om, R, n, m, 1, (tmp[1] + im * tmp[2]) / rndegen[irpt])
                setposition!(om, R, n, m, 2, (tmp[3] + im * tmp[4]) / rndegen[irpt])
                setposition!(om, R, n, m, 3, (tmp[5] + im * tmp[6]) / rndegen[irpt])
            end
        end
        @assert readline(f) == ""
    end
    close(f)

    return om
end


"""
Create TBModel from Wannier90 interface.

`tbfile` should be `seedname_tb.dat` and `wsvecfile` should be `seedname_wsvec.dat`.

This interface accounts for the distance between orbitals, the effect of which is the
same as `use_ws_distance = true` in the Wannier90 input file.
"""
function createmodelwannier(tbfile::String, wsvecfile::String)
    # read wsvec file
    wsvecs = Dict{Vector{Int64},Vector{Vector{Int64}}}()
    wsndegen = Dict{Vector{Int64},Int64}()
    open(wsvecfile) do f
        readline(f) # this line is comment
        while !eof(f)
            foo = readline(f)
            foo == "" && break
            key = map(x -> parse(Int64, x), split(foo))
            ndegen = parse(Int64, readline(f))
            vecs = [map(x -> parse(Int64, x), split(readline(f))) for _ in 1:ndegen]
            wsvecs[key] = vecs
            wsndegen[key] = ndegen
        end
    end

    lat = zeros(3, 3)
    norbits = 0
    hoppings = Dict{Vector{Int64},Matrix{ComplexF64}}()
    positions = Dict{Vector{Int64},Vector{Matrix{ComplexF64}}}()

    open(tbfile) do f
        readline(f) # this line is comment
        for i in 1:3
            lat[:, i] = map(s->parse(Float64, s), split(readline(f)))
        end
        norbits = parse(Int64, readline(f))
        nrpts = parse(Int64, readline(f))
        rndegen = zeros(0)
        while true
            line = readline(f)
            line == "" && break
            rndegen = [rndegen; map(s->parse(Int64, s), split(line))]
        end
        @assert length(rndegen) == nrpts

        for irpt in 1:nrpts
            R = map(s -> parse(Int64, s), split(readline(f)))
            for m in 1:norbits, n in 1:norbits
                tmp = map(s -> parse(Float64, s), split(readline(f))[(end - 1):end])
                wskey = [R..., n, m]
                for R′ in wsvecs[wskey]
                    hopping = get!(hoppings, R + R′, zeros(ComplexF64, norbits, norbits))
                    hopping[n, m] += (tmp[1] + im * tmp[2]) / rndegen[irpt] / wsndegen[wskey]
                end
            end
            @assert readline(f) == ""
        end

        for irpt in 1:nrpts
            R = map(s -> parse(Int64, s), split(readline(f)))
            for m in 1:norbits, n in 1:norbits
                tmp = map(s -> parse(Float64, s), split(readline(f))[(end - 5):end])
                wskey = [R..., n, m]
                for R′ in wsvecs[wskey]
                    position = get!(positions, R + R′, [zeros(ComplexF64, norbits, norbits) for _ in 1:3])
                    position[1][n, m] += (tmp[1] + im * tmp[2]) / rndegen[irpt] / wsndegen[wskey]
                    position[2][n, m] += (tmp[3] + im * tmp[4]) / rndegen[irpt] / wsndegen[wskey]
                    position[3][n, m] += (tmp[5] + im * tmp[6]) / rndegen[irpt] / wsndegen[wskey]
                end
            end
            @assert readline(f) == ""
        end
    end

    tm = TBModel{ComplexF64}(norbits, lat, isorthogonal=true)

    for (R, hopping) in hoppings
        for m in 1:norbits, n in 1:norbits
            sethopping!(tm, R, n, m, hopping[n, m])
        end
    end

    for (R, position) in positions
        for m in 1:norbits, n in 1:norbits
            for α in 1:3
                setposition!(tm, R, n, m, α, position[α][n, m])
            end
        end
    end

    return tm
end

end
