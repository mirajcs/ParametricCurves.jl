# VectorUtils.jl

*A Julia package for computing vector calculus properties of parametric curves*

## Overview

VectorUtils.jl provides comprehensive tools for analyzing parametric vector curves, including:

- **Tangent vectors** (symbolic and numeric)
- **Frenet-Serret frames** (tangent, normal, binormal vectors)
- **Curvature** computation
- **Torsion** computation
- Support for both symbolic (analytical) and numeric computations

## Features

-  **Dual Mode**: Work with both symbolic expressions and numeric evaluations
-  **Complete Frenet Frame**: Compute T, N, B vectors for any parametric curve
-  **Geometric Properties**: Calculate curvature κ(t) and torsion τ(t)
-  **Flexible Input**: Handle 2D and 3D parametric curves
-  **Optimized**: Efficient algorithms for both analytical and numerical methods


## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/mirajcs/VectorUtils")
```

Or in the Pkg REPL mode (press `]`):

```
add https://github.com/mirajcs/VectorUtils
```

## Mathematical Background

For a parametric curve **r**(t) = [x(t), y(t), z(t)], VectorUtils computes:

**Tangent Vector**: **T**(t) = **r**'(t) / ||**r**'(t)||

**Curvature**: κ(t) = ||**r**'(t) × **r**''(t)|| / ||**r**'(t)||³

**Normal Vector**: **N**(t) = **T**'(t) / ||**T**'(t)||

**Binormal Vector**: **B**(t) = **T**(t) × **N**(t)

**Torsion**: τ(t) = (**r**'(t) × **r**''(t)) · **r**'''(t) / ||**r**'(t) × **r**''(t)||²

## Contents

```@contents
Pages = ["getting_started.md", "api/core.md", "examples/basic.md", "theory.md"]
Depth = 2
```

## Index

```@index
```