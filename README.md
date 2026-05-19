# ParametricCurves.jl

[![Build Status](https://github.com/mirajcs/VectorUtils/workflows/CI/badge.svg)](https://github.com/mirajcs/VectorUtils/actions)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mirajcs.github.io/VectorUtils/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20288830.svg)](https://doi.org/10.5281/zenodo.20288830)

A Julia package for computing symbolic and numeric vector properties of parametric curves, including vector operations, Frenet–Serret frames, curvature, and torsion.

## Features

- **Symbolic & Numeric Computation**: Analytical expressions and numerical evaluations
- **Frenet-Serret Frame**: Complete TNB (tangent-normal-binormal) frame calculation
- **Geometric Invariants**: Curvature κ(t) and torsion τ(t)
- **2D & 3D Support**: Handle planar and space curves
- **High Performance**: Optimized algorithms for both modes

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mirajcs/ParametricCurves.jl")
```

## Quick Example

```julia
using ParametricCurves

# Define a circular helix
r = [cos(t), sin(t), t]

# Compute at t = π/4
t₀ = π/4


T = Tangent(r, t, t₀) 
N = Norm(r, t, t₀) 

# Curvature 
κ = Curvature(r, t, t₀)


```

## Symbolic Computation

```julia
using SymPy

@syms t
r_sym = [cos(t), sin(t), t]

# Symbolic tangent
T_sym = Tangent(r_sym, t)

# Symbolic curvature
κ_sym = Curvature(r_sym, t)
```
<!--
## Examples

### Circle
```julia
# Parametric circle in the xy-plane
circle = [cos(t), sin(t), 0]

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
-->

## Documentation

For detailed documentation, examples, and API reference, visit:

📚 **[Documentation](https://mirajcs.github.io/ParametricCurves.jl/)**



## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

MIT License - see LICENSE file for details

<!--

## Related Packages

-->

## Citation

If you use ParametricCurves.jl in your research, please cite:

```bibtex
@software{parametric_curves_jl,
  author = {Miraj Samarakkody},
  title = {ParametricCurves.jl: A Toolkit for Parametric Curves},
  year = {2026},
  publisher = {Zenodo},
  doi = {https://doi.org/10.5281/zenodo.20288830},
  url = {https://github.com/mirajcs/ParametricCurves.jl}
}
```

## Acknowledgement

This work was supported by the HBCU UP Implementation project, Award No. 2510537.
