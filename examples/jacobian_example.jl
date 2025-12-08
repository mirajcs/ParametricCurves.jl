using ParametricCurves
using SymPy

println("=== Jacobian Matrix Examples ===\n")

# Example 1: Basic 3D transformation
println("Example 1: Jacobian of a 3D vector function")
@syms x y z
funcs = [x^2 + y, x*y*z, z^2]
J = Jacobian(funcs, [x, y, z])
println("Function: f(x,y,z) = [x² + y, xyz, z²]")
println("Jacobian matrix:")
for row in J
    println("  ", row)
end
println()

# Example 2: Parametric curve derivative
println("Example 2: Jacobian of a parametric curve (equivalent to derivative)")
@syms t
curve = [cos(t), sin(t), t^2]
J_curve = Jacobian(curve, [t])
println("Curve: r(t) = [cos(t), sin(t), t²]")
println("Derivative dr/dt:")
for row in J_curve
    println("  ", row)
end
println()

# Example 3: Polar to Cartesian transformation
println("Example 3: Jacobian of polar to Cartesian coordinates")
@syms r θ
polar_to_cartesian = [r*cos(θ), r*sin(θ)]
J_polar = Jacobian(polar_to_cartesian, [r, θ])
println("Transformation: (r,θ) → (r·cos(θ), r·sin(θ))")
println("Jacobian matrix:")
for row in J_polar
    println("  ", row)
end
println()

# Example 4: Spherical to Cartesian transformation
println("Example 4: Jacobian of spherical to Cartesian coordinates")
@syms ρ φ
spherical = [ρ*sin(φ)*cos(θ), ρ*sin(φ)*sin(θ), ρ*cos(φ)]
J_spherical = Jacobian(spherical, [ρ, θ, φ])
println("Transformation: (ρ,θ,φ) → spherical coordinates")
println("Jacobian matrix:")
for (i, row) in enumerate(J_spherical)
    println("  Row $i: ", row)
end
println()

# Example 5: Linear transformation
println("Example 5: Jacobian of a linear transformation (constant)")
funcs_linear = [2*x + y, x - z, y + 3*z]
J_linear = Jacobian(funcs_linear, [x, y, z])
println("Linear transformation: [2x+y, x-z, y+3z]")
println("Jacobian matrix (constant):")
for row in J_linear
    println("  ", row)
end
