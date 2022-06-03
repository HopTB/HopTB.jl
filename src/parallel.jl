module Parallel

using Distributed

export ParallelFunction, claim!, stop!, parallel_do

"""
`ParallelFunction` is introduced to dynamically calculate function values in
a parallel way.

Suppose a function `f` has four arguments `x`, `a`, `b` and `c`. The purpose is to evaluate
the function `f` at several different `x`s with fixed argument `a`, `b` and `c`.
A `ParallelFunction` can be constructed as
```julia
pf = ParallelFunction(f, a, b, c; len=128, processes=workers())
```
To evaluate the function `f` at `x1`, `x2` and `x3`, the following statement can
be used
```julia
pf(x1)
pf(x2)
pf(x3)
```
The function will be evaluated in parallel on any possible processes. To obtain
the result:
```julia
for i in 1:3
    println(claim!(pf))
end
```
There is no guaranteed order of the results.
"""
mutable struct ParallelFunction
    jobs::RemoteChannel
    results::RemoteChannel
    nworkers::Int64
    nreceived::Int64
    ntaken::Int64
    isstopped::Bool
end


function do_work(f, jobs, results, args...)
    while true
        x = take!(jobs)
        if isnothing(x) return nothing end
        result = try
            f(x, args...)
        catch err
            RemoteException(CapturedException(err, stacktrace(catch_backtrace())))
        end
        put!(results, result)
    end
end


function start_work(f, jobs, results, args...)
    @async do_work(f, jobs, results, args...)
    return nothing
end


function ParallelFunction(f::Function, args...; len=5*nworkers(), processes=workers())
    jobs = RemoteChannel(()->Channel{Any}(len))
    results = RemoteChannel(()->Channel{Any}(len))
    for p in processes
        remotecall_fetch(start_work, p, f, jobs, results, args...)
    end
    return ParallelFunction(jobs, results, length(processes), 0, 0, false)
end


function (pf::ParallelFunction)(x)
    if isnothing(x)
        error("nothing is the stop instruction to the worker!")
    end
    put!(pf.jobs, x)
    pf.nreceived += 1
    return nothing
end


function claim!(pf::ParallelFunction)
    if pf.isstopped error("The functions is stopped.") end
    result = take!(pf.results)
    if result isa Exception
        throw(result)
    end
    pf.ntaken += 1
    return result
end


function stop!(pf::ParallelFunction)
    for i in 1:pf.nworkers
        put!(pf.jobs, nothing)
    end
    close(pf.jobs)
    close(pf.results)
    pf.isstopped = true
    return nothing
end


function parallel_do(f::Function, collection; batchsize::Int64=1)
    batchsize > 0 || error("batchsize should be larger than 0.")
    if batchsize == 1
        parallel_do_single(f, collection)
    else
        parallel_do_batch(f, collection, batchsize)
    end
end

function parallel_do_single(f::Function, collection)
    pf = ParallelFunction(f)
    @async for item in collection
        pf(item)
    end
    @sync @async for _ in collection
        claim!(pf)
    end
    stop!(pf)
    return nothing
end

function parallel_do_batch(f::Function, collection, batchsize)
    function f_worker(xs)
        for x in xs
            f(x)
        end
        return nothing
    end
    nsplits = length(collection) ÷ batchsize
    pf = ParallelFunction(f_worker)
    @async begin
        for i in 1:nsplits
            pf(collection[(batchsize * (i - 1) + 1):(batchsize * i)])
        end
        if nsplits * batchsize < length(collection)
            pf(collection[(batchsize * nsplits + 1):length(collection)])
        end
    end
    @sync @async begin
        for i in 1:nsplits
            claim!(pf)
        end
        if nsplits * batchsize < length(collection)
            claim!(pf)
        end
    end
    stop!(pf)
    return nothing
end


function parallel_sum(f::Function, collection, initial_value; batchsize::Int64=1)
    if batchsize == 1
        return parallel_sum_single(f, collection, initial_value)
    else
        return parallel_sum_batch(f, collection, initial_value, batchsize)
    end
end

function parallel_sum_single(f::Function, collection, initial_value)
    pf = ParallelFunction(f)
    @async for item in collection
        pf(item)
    end
    result = initial_value
    @sync @async for _ in collection
        result += claim!(pf)
    end
    stop!(pf)
    return result
end

function parallel_sum_batch(f::Function, collection, initial_value, batchsize)
    function f_worker(xs)
        worker_result = initial_value
        for x in xs
            worker_result += f(x)
        end
        return worker_result
    end
    nsplits = length(collection) ÷ batchsize
    pf = ParallelFunction(f_worker)
    @async begin
        for i in 1:nsplits
            pf(collection[(batchsize * (i - 1) + 1):(batchsize * i)])
        end
        if nsplits * batchsize < length(collection)
            pf(collection[(batchsize * nsplits + 1):length(collection)])
        end
    end
    result = initial_value
    @sync @async begin
        for i in 1:nsplits
            result += claim!(pf)
        end
        if nsplits * batchsize < length(collection)
            result += claim!(pf)
        end
    end
    stop!(pf)
    return result
end

end
