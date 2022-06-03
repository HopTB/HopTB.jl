using HopTB, Test

function getBN()
    lat = [1 1/2 0; 0 √3/2 0; 0 0 10.0]
    site_positions = lat*([1/3 1/3 0; 2/3 2/3 0]')
    tm = TBModel(lat, site_positions, [[0], [0]], isorthogonal=true)
    addhopping!(tm, [0, 0, 0], (1, 1), (1, 1), -0.5)
    addhopping!(tm, [0, 0, 0], (2, 1), (2, 1), 0.5)
    addhopping!(tm, [0, 0, 0], (1, 1), (2, 1), 1.0)
    addhopping!(tm, [-1, 0, 0], (1, 1), (2, 1), 1.0)
    addhopping!(tm, [0, -1, 0], (1, 1), (2, 1), 1.0)
    return tm
end


"""
Kane-Mele QSH model
"""
function get_Kane_Mele(;include_rashba=true)
    lat = [1 1/2 0; 0 √3/2 0; 0 0 10.0]
    site_positions = lat*([1/3 1/3 0; 2/3 2/3 0]')
    tm = TBModel(lat, site_positions, [[0], [0]], isorthogonal=true, isspinful=true)
    σ0 = [1 0; 0 1]
    σ1 = [0 1; 1 0]
    σ2 = [0 -im; im 0]
    σ3 = [1 0; 0 -1]

    # onsite energy
    addhopping!(tm, [0, 0, 0], 1, 1, σ0)
    addhopping!(tm, [0, 0, 0], 2, 2, -σ0)

    t = 1.0
    soc = 0.6*0.5*t
    rashba = 0.25*t
    # spin independent hopping
    addhopping!(tm, [0, 0, 0], 1, 2, t*σ0)
    addhopping!(tm, [0, -1, 0], 1, 2, t*σ0)
    addhopping!(tm, [-1, 0, 0], 1, 2, t*σ0)
    # soc
    addhopping!(tm, [0, 1, 0], 1, 1, -im*soc*σ3)
    addhopping!(tm, [1, 0, 0], 1, 1, im*soc*σ3)
    addhopping!(tm, [1, -1, 0], 1, 1, -im*soc*σ3)
    addhopping!(tm, [0, 1, 0], 2, 2, im*soc*σ3)
    addhopping!(tm, [1, 0, 0], 2, 2, -im*soc*σ3)
    addhopping!(tm, [1, -1, 0], 2, 2, im*soc*σ3)
    # Rashba
    if include_rashba
        addhopping!(tm, [0, 0, 0], 1, 2, im*rashba*(0.5*σ1-√3/2*σ2))
        addhopping!(tm, [0, -1, 0], 1, 2, -im*rashba*σ1)
        addhopping!(tm, [-1, 0, 0], 1, 2, im*rashba*(0.5*σ1+√3/2*σ2))
    end

    return tm
end


function getHaldane()
    lat = [1 1/2 0; 0 √3/2 0; 0 0 10.0]
    site_positions = lat*([1/3 1/3 0; 2/3 2/3 0]')
    haldane = HopTB.TBModel(lat, site_positions, [[0], [0]])

    Δ = 0.0
    t1 = -1.0
    t2 = 0.15im
    t2c = conj(t2)

    addhopping!(haldane, [0, 0, 0], (1, 1), (1, 1), -Δ)
    addhopping!(haldane, [0, 0, 0], (2, 1), (2, 1), Δ)

    addhopping!(haldane, [0, 0, 0], (1, 1), (2, 1), t1)
    addhopping!(haldane, [1, 0, 0], (2, 1), (1, 1), t1)
    addhopping!(haldane, [0, 1, 0], (2, 1), (1, 1), t1)

    addhopping!(haldane, [1, 0, 0], (1, 1), (1, 1), t2)
    addhopping!(haldane, [1, -1, 0], (2, 1), (2, 1), t2)
    addhopping!(haldane, [0, 1, 0], (2, 1), (2, 1), t2)
    addhopping!(haldane, [1, 0, 0], (2, 1), (2, 1), t2c)
    addhopping!(haldane, [1, -1, 0], (1, 1), (1, 1), t2c)
    addhopping!(haldane, [0, 1, 0], (1, 1), (1, 1), t2c)

    return haldane
end

function getGraphene()
    # graphene model
    lat = [1 1/2 0; 0 √3/2 0; 0 0 1]
    site_positions = lat*([1/3 1/3 0; 2/3 2/3 0]');
    tm = TBModel(lat,site_positions,[[0], [0]])
    @test tm.rlat ≈ [2π 0 0; -2π/√3 4π/√3 0; 0 0 2π]
    #sethopping!(tm, [0, 0, 0], 1, 1, 0.0)
    #sethopping!(tm, [0, 0, 0], 2, 2, 0.0)
    sethopping!(tm, [0, 0, 0], 2, 1, -3.0)
    sethopping!(tm, [1, 0, 0], 2, 1, -3.0)
    sethopping!(tm, [0, 1, 0], 2, 1, -3.0)
    return tm
end

"""
```julia
function getBHZ(;m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
```
Construct a 4-band tight-binding BHZ model. This model has same band structure
near Γ point as the 4-band k⋅p BHZ model, see function [`BHZ`](@ref) below.
"""
function getBHZ(;m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
    lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    tm = TBModel(4, lat)
    for i=1:4
        setposition!(tm, [0, 0, 0], i, i, [0.0,0.0,0.0])
    end

    sethopping!(tm, [0, 0, 0], 1, 1, m-4b)
    sethopping!(tm, [0, 0, 0], 2, 2, 4b-m)
    sethopping!(tm, [0, 0, 0], 3, 3, m-4b)
    sethopping!(tm, [0, 0, 0], 4, 4, 4b-m)

    sethopping!(tm, [1, 0, 0], 1, 1, b)
    sethopping!(tm, [1, 0, 0], 2, 2, -b)
    sethopping!(tm, [1, 0, 0], 3, 3, b)
    sethopping!(tm, [1, 0, 0], 4, 4, -b)
    sethopping!(tm, [1, 0, 0], 1, 2, -im/2*a)
    sethopping!(tm, [1, 0, 0], 2, 1, -im/2*a)
    sethopping!(tm, [1, 0, 0], 3, 4, im/2*a)
    sethopping!(tm, [1, 0, 0], 4, 3, im/2*a)

    sethopping!(tm, [0, 1, 0], 1, 1, b)
    sethopping!(tm, [0, 1, 0], 2, 2, -b)
    sethopping!(tm, [0, 1, 0], 3, 3, b)
    sethopping!(tm, [0, 1, 0], 4, 4, -b)
    sethopping!(tm, [0, 1, 0], 1, 2, -a/2)
    sethopping!(tm, [0, 1, 0], 2, 1, a/2)
    sethopping!(tm, [0, 1, 0], 3, 4, -a/2)
    sethopping!(tm, [0, 1, 0], 4, 3, a/2)

    return tm
end

"""
```julia
function getHalfBHZ(;m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
```
Construct a 2-band tight-binding BHZ model. This model has same band structure
near Γ point as the 2-band k⋅p BHZ model, see function [`HalfBHZ`](@ref) below.
"""
function getHalfBHZ(;m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
    lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    tm = TBModel(2, lat)
    for i=1:2
        setposition!(tm, [0, 0, 0], i, i, [0.0,0.0,0.0])
    end

    sethopping!(tm, [0, 0, 0], 1, 1, m-4b)
    sethopping!(tm, [0, 0, 0], 2, 2, 4b-m)

    sethopping!(tm, [1, 0, 0], 1, 1, b)
    sethopping!(tm, [1, 0, 0], 2, 2, -b)
    sethopping!(tm, [1, 0, 0], 1, 2, -im/2*a)
    sethopping!(tm, [1, 0, 0], 2, 1, -im/2*a)

    sethopping!(tm, [0, 1, 0], 1, 1, b)
    sethopping!(tm, [0, 1, 0], 2, 2, -b)
    sethopping!(tm, [0, 1, 0], 1, 2, -a/2)
    sethopping!(tm, [0, 1, 0], 2, 1, a/2)

    return tm
end

"""
```julia
function getTI3D(;e0=-0.0068, e1=-1.3, e2=-19.6,
    a1=2.2,a2=4.1,m0=0.28,m1=-10,m2=-56.6,)
```
Construct a 4-band tight-binding 3DTI model. Its k⋅p conterpart see function
[`TI3D`](@ref) below.
The default values are from Bi2Se3 paper.
"""
function getTI3D(;e0=-0.0068, e1=-1.3, e2=-19.6, a1=2.2,a2=4.1,m0=0.28,m1=-10,m2=-56.6,)
    lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    tm = TBModel(4, lat)
    for i=1:4
        setposition!(tm, [0, 0, 0], i, i, [0.0,0.0,0.0])
    end

    sethopping!(tm, [0, 0, 0], 1, 1, e0+2*e1+4*e2+m0+2m1+4m2)
    sethopping!(tm, [0, 0, 0], 2, 2, e0+2*e1+4*e2-m0-2m1-4m2)
    sethopping!(tm, [0, 0, 0], 3, 3, e0+2*e1+4*e2+m0+2m1+4m2)
    sethopping!(tm, [0, 0, 0], 4, 4, e0+2*e1+4*e2-m0-2m1-4m2)

    sethopping!(tm, [1, 0, 0], 1, 4, -im*a2/2)
    sethopping!(tm, [1, 0, 0], 2, 3, -im*a2/2)
    sethopping!(tm, [1, 0, 0], 3, 2, -im*a2/2)
    sethopping!(tm, [1, 0, 0], 4, 1, -im*a2/2)
    sethopping!(tm, [1, 0, 0], 1, 1, -m2-e2)
    sethopping!(tm, [1, 0, 0], 2, 2, m2-e2)
    sethopping!(tm, [1, 0, 0], 3, 3, -m2-e2)
    sethopping!(tm, [1, 0, 0], 4, 4, m2-e2)

    sethopping!(tm, [0, 1, 0], 1, 4, -a2/2)
    sethopping!(tm, [0, 1, 0], 2, 3, -a2/2)
    sethopping!(tm, [0, 1, 0], 3, 2, a2/2)
    sethopping!(tm, [0, 1, 0], 4, 1, a2/2)
    sethopping!(tm, [0, 1, 0], 1, 1, -m2-e2)
    sethopping!(tm, [0, 1, 0], 2, 2, m2-e2)
    sethopping!(tm, [0, 1, 0], 3, 3, -m2-e2)
    sethopping!(tm, [0, 1, 0], 4, 4, m2-e2)

    sethopping!(tm, [0, 0, 1], 1, 2, -im*a1/2)
    sethopping!(tm, [0, 0, 1], 2, 1, -im*a1/2)
    sethopping!(tm, [0, 0, 1], 3, 4, im*a1/2)
    sethopping!(tm, [0, 0, 1], 4, 3, im*a1/2)
    sethopping!(tm, [0, 0, 1], 1, 1, -m1-e1)
    sethopping!(tm, [0, 0, 1], 2, 2, m1-e1)
    sethopping!(tm, [0, 0, 1], 3, 3, -m1-e1)
    sethopping!(tm, [0, 0, 1], 4, 4, m1-e1)
    return tm
end

"""
```julia
function getTI2D(;vf=2.36, m0=-0.029, m1=12.9)
```
Construct a 4-band tight-binding 2DTI model. Its k⋅p conterpart see function
[`TI2D`](@ref) below.
The default values are from wang jing's paper.
To keep with the tradition of wannier, the basis is arranged as
u↑, l↑, u↓, l↓.
"""
function getTI2D(;vf=2.36, m0=-0.029, m1=12.9)
    lat = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    tm = TBModel(4, lat)
    for i=1:4
        setposition!(tm, [0, 0, 0], i, i, [0.0,0.0,0.0])
    end

    sethopping!(tm, [0, 0, 0], 1, 2, m0+4m1)
    #sethopping!(tm, [0, 0, 0], 2, 1, m0+4m1)
    sethopping!(tm, [0, 0, 0], 3, 4, m0+4m1)
    #sethopping!(tm, [0, 0, 0], 4, 3, m0+4m1)

    sethopping!(tm, [1, 0, 0], 1, 3, vf/2)
    sethopping!(tm, [1, 0, 0], 2, 4, -vf/2)
    sethopping!(tm, [1, 0, 0], 3, 1, -vf/2)
    sethopping!(tm, [1, 0, 0], 4, 2, vf/2)
    sethopping!(tm, [1, 0, 0], 1, 2, -m1)
    sethopping!(tm, [1, 0, 0], 2, 1, -m1)
    sethopping!(tm, [1, 0, 0], 3, 4, -m1)
    sethopping!(tm, [1, 0, 0], 4, 3, -m1)

    sethopping!(tm, [0, 1, 0], 1, 3, -im*vf/2)
    sethopping!(tm, [0, 1, 0], 2, 4, im*vf/2)
    sethopping!(tm, [0, 1, 0], 3, 1, -im*vf/2)
    sethopping!(tm, [0, 1, 0], 4, 2, im*vf/2)
    sethopping!(tm, [0, 1, 0], 1, 2, -m1)
    sethopping!(tm, [0, 1, 0], 2, 1, -m1)
    sethopping!(tm, [0, 1, 0], 3, 4, -m1)
    sethopping!(tm, [0, 1, 0], 4, 3, -m1)

    return tm
end

σ0 = [1 0; 0 1]
σ1 = [0 1; 1 0]
σ2 = [0 -im; im 0]
σ3 = [1 0; 0 -1]
"""
Construct k⋅p model for a three dimensional topological insulators. The default
values are from Bi2Se3 paper.
"""
function TI3D(kpt; e0=-0.0068, e1=-1.3, e2=-19.6, a1=2.2,a2=4.1,m0=0.28,m1=-10,m2=-56.6,)
    (kx,ky,kz) = kpt
    ham = zeros(ComplexF64, (4, 4))
    ham += (e0+e1*kz^2+e2*(kx^2+ky^2))*kron(σ0, σ0)
    ham += (m0+m1*kz^2+m2*(kx^2+ky^2))*kron(σ0, σ3)
    ham += a1*kz*kron(σ3, σ1)
    ham += a2*kx*kron(σ1, σ1)
    ham += a2*ky*kron(σ2, σ1)
    return ham
end
"""
Construct k⋅p model for Half of BHZ model. The default
values are topological trivial with m*b<0.
"""
function HalfBHZ(kpt; m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
    (kx,ky) = kpt
    ham = zeros(ComplexF64, (2, 2))
    ham += a*kx*σ1
    ham += a*ky*σ2
    ham += (m-b*(kx^2+ky^2))*σ3
    return ham
end
"""
Construct k⋅p model for 4-band BHZ model. The default
values are topological trivial with m*b<0.
"""
function BHZ(kpt, m::Float64=1.0, b::Float64=-0.2, a::Float64=0.2)
    (kx,ky) = kpt
    ham = zeros(ComplexF64, (4, 4))
    ham += (m-b*(kx^2+ky^2))*kron(σ0,σ3)
    ham += a*kx*kron(σ3,σ1)
    ham += a*ky*kron(σ0,σ2)
    return ham
end

"""
Construct k⋅p model for two-dimensional TI model. The default
values are from wang jing's paper.
To keep with the tradition of wannier, the basis is arranged as
u↑, l↑, u↓, l↓.
"""
function  TI2D(kpt; vf=2.36, m0=-0.029, m1=12.9)
    (kx,ky) = kpt
    ham = zeros(ComplexF64, (4, 4))
    ham += vf*ky*kron(σ1, σ3)
    ham += -vf*kx*kron(σ2, σ3)
    ham += (m0+m1*(kx^2+ky^2))*kron(σ0, σ1)
    return ham
end
