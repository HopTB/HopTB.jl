cd(dirname(dirname(@__FILE__)))
@info "CWD is " * pwd()

@info "Loading Package ..."
using HopTB, Test
using Downloads

@info "Loading Data ..."

models = Dict([("Si.dat", :aims),
    ("WS2.scfout", :openmx38),
    ("WS2-soc.scfout", :openmx38),
    ("WS2.w90", :wannier),
    ("WS2.dat", :aims),
    ("SnH.openmx", :openmx38),
    ("WS2.openmx", :openmx38),
    ("Fe.openmx", :openmx38),
    ("GaAs.openmx", :openmx38),
    ("GaAs.openmx39", :openmx)])
tasks = Dict{String,Vector{NamedTuple{(:name, :func),Tuple{String,Function}}}}()
for k in keys(models)
    tasks[k] = Vector()
end
macro register(datafile::String, callback)
    k = datafile
    cb = eval(callback)
    if !haskey(tasks, k)
        throw("Data file $datafile not found!")
    end
    return quote
        push!(tasks[$k], $(esc((name = basename(string(__source__.file)), func = cb))))
    end
end

if !ispath(pwd() * "/test/data")
    mkdir(pwd() * "/test/data")
end
for filename in keys(models)
    filepath = pwd() * "/test/data/$filename"
    if !ispath(filepath)
        Downloads.download("https://raw.githubusercontent.com/HopTB/HopTestData/main/$filename", pwd() * "/test/data/$filename")
    end
end

@info "Collecting Tests ..."
include("memoize.jl")
include("model.jl")
include("bs.jl")
include("basic.jl")
include("optics.jl")
include("hall.jl")
include("group.jl")
include("wannier.jl")
include("topology.jl")
include("utilities.jl")
include("response.jl")
include("structure.jl")
include("magnetism.jl")
include("shared.jl")
include("mesh.jl")

@info "Running Tests ..."
for (datafile, taskvec) in tasks
    @info "  Loading $datafile ..."
    method = models[datafile]
    nm = getfield(HopTB.Interface, Symbol(:createmodel, method))("test/data/" * datafile)
    for (sourcefilename, testfunc) in taskvec
        @info "    Running test in $sourcefilename ..."
        testfunc(nm)
    end
    GC.gc()
end
