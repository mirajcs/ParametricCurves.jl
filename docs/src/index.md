# GeometryToolkit.jl

A Julia package for geometric computations.

## Overview

GeometryToolkit.jl provides a convenient and efficient suite of functions for geometric computations in Julia, including:

- Vector operations: Norm, normalization, dot products, and cross products
- Geometric computations: Angles between vectors, projections, and distance calculations
- Parametric geometry: Line and plane equations in 3D space
- Curve analysis: Arc length, tangent vectors, curvature, normal vectors, binormal vectors, torsion and Frenet-Serret frame. 
- Symbolic and numeric support: All functions work with both numeric values and symbolic expressions via SymPy

Whether you're doing computational geometry, differential geometry, or general vector calculus, GeometryToolkit.jl provides the essential tools you need.

## Installation
```julia
using Pkg
Pkg.add("GeometryToolkit")
```

For symbolic computation support, also install SymPy.jl
```julia
Pkg.add("SymPy")
```

## Quick Start
```julia
using GeometryToolkit
using SymPy

# Example usage
julia> @syms x y z l m n 
julia> Cross([x,y,z],[l,m,n])

[
y*n - z*m,
z*l - x*n,
x*m - y*l
]

```

## Documentation Contents
```@contents
Pages = [
    "api.md"
]
Depth = 2
```

## Index
```@index
```

