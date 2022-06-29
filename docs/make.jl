using Documenter, HopTB

makedocs(
    sitename="HopTB.jl",
    pages = [
        "Introduction" => "index.md",
        "Ab initio calculation" => "abinitio_tutorial.md",
        "Model calculation" => "model_tutorial.md",
        "Features and API references" => [
            "HopTB" => "hop.md",
            "Optics" => "optics.md",
            "Hall effects" => "hall.md",
            "Magnetism" => "magnetism.md",
            "Band Structure" => "bs.md",
            "Topology" => "topology.md",
            "Symmetrization of TB model" => "group.md"
        ]
    ]
)

deploydocs(
    repo = "github.com/HopTB/HopTB.jl.git",
    devbranch = "main"
)
