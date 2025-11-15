using Documenter
using GeometryToolkit  # Was: VectorUtils

makedocs(
    sitename = "GeometryToolkit.jl",
    modules = [GeometryToolkit],  # Was: VectorUtils
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md"
    ],
)