module Memoize

using MacroTools

export @memoize

# settings

const MEMOIZE_ENABLED = [true]

@doc raw"""
```julia
enable_memoization()
```

Enables memoization.
"""
function enable_memoization()
    MEMOIZE_ENABLED[1] = true
end

@doc raw"""
```julia
disable_memoization()
```

Disables memoization.
"""
function disable_memoization()
    MEMOIZE_ENABLED[1] = false
end

# memoization

@doc raw"""
Gather all the input variables from a MacroTools-splitted function.
All the variables are escaped.

# Example

```julia
julia> f = :(
           function foo(x, y::Int64=2; z::Int64=3)::Int64
               return x + y + z
           end
       );

julia> gather_inputs(MacroTools.splitdef(f))
:(($(Expr(:escape, :x)), $(Expr(:escape, :y)), $(Expr(:escape, :z))))
```
"""
function gather_inputs(f_dict)
    key = :(())
    for arg in f_dict[:args]
        push!(key.args, esc(MacroTools.splitarg(arg)[1]))
    end
    for kwarg in f_dict[:kwargs]
        push!(key.args, esc(MacroTools.splitarg(kwarg)[1]))
    end
    return key
end


@doc raw"""
Construct an expression from a MacroTools-splitted function
definition to call the function.

# Example

```julia
julia> f = :(
           function foo(x, y::Int64=2; z::Int64=3)::Int64
               return x + y + z
           end
       );

julia> call_function(MacroTools.splitdef(f))
:(foo(x, y, z = z))
```
"""
function call_function(f_dict)
    call_expr = Expr(:call, f_dict[:name])
    for arg in f_dict[:args]
        argname = MacroTools.splitarg(arg)[1]
        push!(call_expr.args, argname)
    end
    for kwarg in f_dict[:kwargs]
        kwargname = MacroTools.splitarg(kwarg)[1]
        push!(call_expr.args, Expr(:kw, kwargname, kwargname))
    end
    return call_expr
end


mutable struct Cache
    leadkey::Any
    data::Dict{Any,Any}
end


function _memoize(leadkey, f::Expr)
    f_unmemoized_dict = MacroTools.splitdef(f)
    f_memoized_dict = deepcopy(f_unmemoized_dict)
    f_unmemoized_dict[:name] = Symbol("##_", f_unmemoized_dict[:name], "_unmemoized")
    cache = esc(Symbol("##_", f_memoized_dict[:name], "_cache"))
    call_expr = esc(call_function(f_unmemoized_dict))
    key = gather_inputs(f_memoized_dict)

    # escape most parts of f_memoized
    for item in [:args, :kwargs, :whereparams]
        if haskey(f_memoized_dict, item)
            f_memoized_dict[item] = map(esc, f_memoized_dict[item])
        end
    end
    for item in [:name, :rtype]
        if haskey(f_memoized_dict, item)
            f_memoized_dict[item] = esc(f_memoized_dict[item])
        end
    end

    f_memoized_dict[:body] = quote
        if !MEMOIZE_ENABLED[1]
            return $call_expr
        else
            if !isequal($(esc(leadkey)), $(cache).leadkey)
                empty!($(cache).data)
                $(cache).leadkey = $(esc(leadkey))
            end
            return Base.get!($(cache).data, $key) do
                $call_expr
            end
        end
    end
    return quote
        $(cache) = Cache(nothing, Dict())
        $(esc(MacroTools.combinedef(f_unmemoized_dict)))
        Base.@__doc__ $(MacroTools.combinedef(f_memoized_dict))
    end
end

@doc raw"""
Make functions remember their results.

# Usage:

## Example 1

```julia
@memoize x function f(x::Vector{Int64}, y::Int64; z::Float64=3.0)::Float64
    return sum(x) + y + z
end
```

The function `f` will remember its result for the same input. However,
if `f` is called with a different `x` from previous function calls, the cache will be emptied.

## Example 2

```julia
@memoize function f(x::Vector{Int64}, y::Int64; z::Float64=3.0)::Float64
    return sum(x) + y + z
end
```
is equivalent to
```julia
@memoize nothing function f(x::Vector{Int64}, y::Int64; z::Float64=3.0)::Float64
    return sum(x) + y + z
end
```
and `f` will remember every result forever.

# Notice

 - This macro uses a dictionary as the cache. All inputs of the function (including keyword arguments)
   are used as dictionary keys. Modifying dictionary keys in julia is currently undefined behavior.
   Therefore, function arguments should never be modified, either in or out of the function.

 - The return type of the function should always be explicitly declared. Otherwise the compiler will
   have no information of the return type since the result may be retrieved from the cache.

 - This macro does not support multithreading.
"""
macro memoize(leadkey, f)
    return _memoize(leadkey, f)
end

macro memoize(f)
    return _memoize(nothing, f)
end

end
