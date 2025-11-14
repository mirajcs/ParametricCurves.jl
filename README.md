# VectorUtils.jl

[![Build Status](https://github.com/mirajcs/VectorUtils/workflows/CI/badge.svg)](https://github.com/mirajcs/VectorUtils/actions)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mirajcs.github.io/VectorUtils/)

A Julia package for computing symbolic and numeric vector properties of parametric curves, including tangent vectors, Frenet-Serret frames, curvature, and torsion.

## Features

- ✨ **Symbolic & Numeric Computation**: Analytical expressions and numerical evaluations
- 📐 **Frenet-Serret Frame**: Complete TNB (tangent-normal-binormal) frame calculation
- 📊 **Geometric Invariants**: Curvature κ(t) and torsion τ(t)
- 🎯 **2D & 3D Support**: Handle planar and space curves
- 🚀 **High Performance**: Optimized algorithms for both modes

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mirajcs/VectorUtils")
```

## Quick Example

```julia
using VectorUtils

# Define a circular helix
r(t) = [cos(t), sin(t), t]

# Compute at t = π/4
t₀ = π/4

# Tangent vector
T = tangent(r, t₀)

# Full Frenet frame
frame = frenet_frame(r, t₀)
println("Tangent:  ", frame.tangent)
println("Normal:   ", frame.normal)
println("Binormal: ", frame.binormal)

# Curvature and torsion
κ = curvature(r, t₀)
τ = torsion(r, t₀)
println("Curvature: κ = ", κ)
println("Torsion:   τ = ", τ)
```

## Symbolic Computation

```julia
using Symbolics

@variables t
r_sym = [cos(t), sin(t), t]

# Symbolic tangent
T_sym = tangent_symbolic(r_sym, t)

# Symbolic curvature
κ_sym = curvature_symbolic(r_sym, t)
```

## Examples

### Circle
```julia
# Parametric circle in the xy-plane
circle(t) = [cos(t), sin(t), 0]

κ = curvature(circle, 0.0)  # κ = 1 (constant)
τ = torsion(circle, 0.0)    # τ = 0 (planar curve)
```

### Helix
```julia
# Circular helix
helix(t) = [a*cos(t), a*sin(t), b*t]

κ = curvature(helix, π)     # κ = a/(a² + b²)
τ = torsion(helix, π)       # τ = b/(a² + b²)
```

### Viviani's Curve
```julia
# Intersection of sphere and cylinder
viviani(t) = [cos(t)^2, cos(t)*sin(t), sin(t)]

frame = frenet_frame(viviani, π/6)
```

## Documentation

For detailed documentation, examples, and API reference, visit:

📚 **[Documentation](https://mirajcs.github.io/VectorUtils/)**

## Mathematical Background

VectorUtils implements the Frenet-Serret formulas for space curves:

- **Tangent**: T(t) = r'(t) / ||r'(t)||
- **Curvature**: κ(t) = ||r'(t) × r''(t)|| / ||r'(t)||³
- **Normal**: N(t) = T'(t) / ||T'(t)||
- **Binormal**: B(t) = T(t) × N(t)
- **Torsion**: τ(t) = (r'(t) × r''(t)) · r'''(t) / ||r'(t) × r''(t)||²

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see LICENSE file for details

## Related Packages

- [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) - Symbolic computation
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) - ODE solving
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) - Visualization

## Citation

If you use VectorUtils.jl in your research, please cite:

```bibtex
@software{vectorutils_jl,
  author = {mirajcs},
  title = {VectorUtils.jl: Vector Calculus Tools for Parametric Curves},
  year = {2024},
  url = {https://github.com/mirajcs/VectorUtils}
}
```