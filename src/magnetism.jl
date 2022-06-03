module Magnetism

using LinearAlgebra, Distributed
using ..HopTB

@doc raw"""
```julia
get_orbital_moment(tm::AbstractTBModel, α::Int64, k::Vector{Float64})
    --> M::Matrix{ComplexF64}
```

This function returns orbital magnetic moment in the unit of Bohr magneton

```math
M_{nm}=i m ϵ^{αβγ} \sum_{\bar{l}} \frac{v_{nl}^β v_{lm}^γ}{ϵ_n-ϵ_l},
```

where either n = m or band n and band m are completely degenerate.

`M[n, m]` would be set to zero if the above condition is not fulfilled.

!!! warning
    This function has not been tested.
"""
function get_orbital_moment(tm::AbstractTBModel, α::Int64, k::Vector{Float64})
    if α == 1
        result = _get_orbital_moment_core(tm, 2, 3, k)-_get_orbital_moment_core(tm, 3, 2, k)
    elseif α == 2
        result = _get_orbital_moment_core(tm, 3, 1, k)-_get_orbital_moment_core(tm, 1, 3, k)
    elseif α == 3
        result = _get_orbital_moment_core(tm, 1, 2, k)-_get_orbital_moment_core(tm, 2, 1, k)
    else
        error("Invalid value of α.")
    end
    return 0.1312342118*im*result # m*e/ħ^2*10^-20
end

function _get_orbital_moment_core(tm::AbstractTBModel, β::Int64, γ::Int64, k::Vector{Float64})
    result = zeros(ComplexF64, tm.norbits, tm.norbits)
    vβ = getvelocity(tm, β, k); vγ = getvelocity(tm, γ, k)
    Es = geteig(tm, k).values
    for n in 1:tm.norbits, m in 1:tm.norbits
        if abs(Es[n]-Es[m]) <= HopTB.DEGEN_THRESH[1]
            for l in 1:tm.norbits
                if abs(Es[n]-Es[l]) > HopTB.DEGEN_THRESH[1]
                    result[n, m] += vβ[n, l]*vγ[l, m]/(Es[n]-Es[l])
                end
            end
        end
    end
    return result
end


"""
```julia
get_field_modified_Es(tm::AbstractTBModel, α::Int64, B::Float64,
    k::Vector{Float64}; double_degenerate=false)
    --> mEs::Vector{Float64}
```

Calculate modified energy due to orbital and spin Zeeman correction.
If double_degenerate is true, the bands are assumed to be doubly degenerate.
Otherwise the bands are assumed to be nondegenerate.

!!! warning
    This function has not been tested.
"""
function get_field_modified_Es(tm::AbstractTBModel, α::Int64, B::Float64,
    k::Vector{Float64}; double_degenerate=false)
    spin_moment = -5.8e-5*getspin(tm, α, k) # negative charge of electron and Bohr magneton
    orbital_moment = 5.8e-5*get_orbital_moment(tm, α, k) # Bohr magneton
    mEs = deepcopy(geteig(tm, k).values)
    if double_degenerate
        for n in 1:tm.norbits÷2
            moment = spin_moment[2n-1:2n, 2n-1:2n]+orbital_moment[2n-1:2n, 2n-1:2n]
            mEs[2n-1:2n] += -real(eigvals(moment)*B)
        end
    else
        for n in 1:tm.norbits
            mEs[n] += -real((spin_moment[n, n]+orbital_moment[n, n])*B)
        end
    end
    return mEs
end

end
