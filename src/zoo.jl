module Zoo

using ..Hop

"""
A two band model for boron nitride.

Δ is onsite energy and t is hopping between nearest neighbours.
"""
function getBN(;Δ::Float64=0.5, t::Float64=1.0)
    lat = [1 1/2 0; 0 √3/2 0; 0 0 10.0]
    orbital_positions = lat * ([1/3 1/3 0; 2/3 2/3 0]')
    tm = TBModel(lat, orbital_positions, isorthogonal=true)
    # onsite energy
    addhopping!(tm, [0, 0, 0], 1, 1, -Δ)
    addhopping!(tm, [0, 0, 0], 2, 2, Δ)
    # nearest neighbour hopping
    for R in [[0, 0, 0], [-1, 0, 0], [0, -1, 0]]
        addhopping!(tm, R, 1, 2, t)
    end
    return tm
end


"""
A minimal 3D TB model: cubic lattice, one orbital per unit cell and
nearest neighbour hopping `t`.
"""
function getcube(;t::Real=1.0)
    lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    orbital_positions = reshape([0.0, 0.0, 0.0], (3, 1))
    tm = TBModel(lat, orbital_positions, isorthogonal=true)
    # nearest neighbour hopping
    for R in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        addhopping!(tm, R, 1, 1, t)
    end
    return tm
end

end