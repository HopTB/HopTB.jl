using Hop, Test
using LinearAlgebra
include("zoo.jl")

let
    # stability test
    tm = getBN()
    # This observable is simply the charge density on the A sublattice
    function getO(k)
        return [1.0 0.0; 0.0 0.0]
    end
    # This perturbation simply add onsite energy to the A sublattice
    function getV(k)
        return [0.2 0.0; 0.0 0.0]
    end
    # Fermi distribution for chemical potental=0
    function getn(ω)
        return 1/(exp(400*ω)+1)
    end
    # Charge density on A sublattice should decrease since the energy is higher
    @test Hop.Response.getresponse(k->getH(tm, k), getO, getV, getn, [10, 10, 1]) ≈ -0.029 atol=1.0e-3
end


let
    function getO(k)
        return kron([1.0 0.0; 0.0 -1.0], I(2))
    end
    # This perturbation is also s_z.
    function getV(k)
        return kron([1.0 0.0; 0.0 -1.0], I(2))
    end
    function getn(ω)
        return 1/(exp(400*ω)+1)
    end
    tm = get_Kane_Mele()
    ftm = Hop.Floquet.addlight(tm, [1.0,-im, 0], 2.0, 1)
    ks = Hop.Utilities.constructmeshkpts([5, 5, 1])
    @test Hop.Response.get_Floquet_response(k->getH(ftm, k), getO, getV, getn, 0, [5, 5, 1])[1] ≈
        sum(Hop.Response.clt_Floquet_response(k->getH(ftm, k), getO, getV, getn, 0, ks)[2])/25
end


let
    # This observable is s_z.
    # The basis for TI2D is arranged by uup, dup, udn, ldn
    function getO(k)
        return kron([1.0 0.0; 0.0 -1.0], I(2))
    end
    # This perturbation is also s_z.
    function getV(k)
        return kron([1.0 0.0; 0.0 -1.0], I(2))
    end
    # Fermi distribution for chemical potental=0
    function getn(ω)
        return 1/(exp(400*ω)+1)
    end

    # test on Floquet Tight-binding TI2D
    tm=getTI2D()
    w=2.0;
    forder=1;
    avec=0.0*[1.0,-im,0];
    ftm = Hop.Floquet.addlight(tm, avec, w, forder)
    function getH′(k)
        b′=inv(ftm.tm.rlat)
        return getH(ftm,Vector(b′*k))
    end
    # tests without light, check if consistent with without light
    @test Hop.Response.get_Floquet_response(getH′, getO, getV, getn, 0, [0.1,0.0,0.0],Γ=0.0001)[1] ≈ -6.62313 atol=1e-4
    @test Hop.Response.get_Floquet_response(getH′, getO, getV, getn, 1, [0.1,0.0,0.0],Γ=0.0001)[1] ≈ 0.0 atol=1e-4
    @test Hop.Response.get_Floquet_response(getH′, getO, getV, getn, 0, [10, 10, 1],Γ=0.001;
        k1=[-0.1,-0.1,0.0],k2=[0.1,0.1,0.0])[1]*0.04 ≈ -0.45955 atol=1e-4

    avec=0.1*[1.0,-im,0];
    ftm = Hop.Floquet.addlight(tm, avec, w, forder)
    function getH′′(k)
        b′=inv(ftm.tm.rlat)
        return getH(ftm,Vector(b′*k))
    end
    @test Hop.Response.get_Floquet_response(getH′′, getO, getV, getn, 0, [0.1,0.0,0.0],Γ=0.0001)[1] ≈ -2.72846 atol=1e-4
    @test Hop.Response.get_Floquet_response(getH′′, getO, getV, getn, 1, [0.1,0.0,0.0],Γ=0.0001)[1] ≈ -0.17655 atol=1e-4
    @test Hop.Response.get_Floquet_response(getH′′, getO, getV, getn, 0, [10, 10, 1],Γ=0.001;
        k1=[-0.1,-0.1,0.0],k2=[0.1,0.1,0.0])[1]*0.04 ≈ -0.13986 atol=1e-4
end
