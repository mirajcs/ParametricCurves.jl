module VectorUtils

greet() = print("Hello Geometrician")

export Norm, Normalize, Dot, Cross, Angle, Projection, ParametricLine, PlaneEquation, ArcLength, ArcLengthParametrization, Tangent, Curvature, Normal, Binormal, Torsion, FrenetSerret

using LinearAlgebra, SymPy, QuadGK

"""
    Norm(v::AbstractVector)

Compute the Euclidean norm (magnitude) of vector `v`.

Supports both numeric and symbolic vectors.

# Examples
```julia
julia> Norm([3.0, 4.0])
5.0

julia> @syms x y
julia> Norm([x, y])
sqrt(x^2 + y^2)
```
"""
Norm(v::AbstractVector{<:Number}) = norm(v)

function Norm(v::AbstractVector{<:Sym})
    return simplify(sqrt(sum(x->x^2, v)))
end

"""
    Normalize(v::AbstractVector)

Return the unit vector in the direction of `v`.

Supports both numeric and symbolic vectors.

# Examples
```julia
julia> Normalize([3.0, 4.0])
2-element Vector{Float64}:
 0.6
 0.8
```
"""
function Normalize(v::AbstractVector{<:Real})
    n = Norm(v)
    if n == 0
        error("Zero divisor")
    end 
    normalize = [v[i] / n for i in 1:length(v)]
    return normalize
end 

function Normalize(v::AbstractVector{<:Sym})
    n = Norm(v)
    normalize = [v[i] / n for i in 1:length(v)]
    return simplify(normalize)
end

"""
    Dot(a::AbstractVector, b::AbstractVector)

Compute the dot product of two vectors `a` and `b`.

Supports both numeric and symbolic vectors.

# Examples
```julia
julia> Dot([1, 2, 3], [4, 5, 6])
32
```
"""
Dot(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) = dot(a, b)

function Dot(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    @assert length(a) == length(b)
    return sum(a[i]*b[i] for i in 1:length(a))
end

"""
    Cross(a::AbstractVector, b::AbstractVector)

Compute the cross product of two 3D vectors `a` and `b`.

Supports both numeric and symbolic vectors.

# Examples
```julia
julia> Cross([1, 0, 0], [0, 1, 0])
3-element Vector{Int64}:
 0
 0
 1
```
"""
function Cross(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    if length(a) != 3 || length(b) != 3
        error("Cross product is only defined for 3D Vectors")
    end 
    return cross(a, b)
end

function Cross(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    @assert length(a) == 3 && length(b) == 3
    return [
        a[2]*b[3] - a[3]*b[2],
        a[3]*b[1] - a[1]*b[3],
        a[1]*b[2] - a[2]*b[1]
    ]
end

"""
    Angle(a::AbstractVector, b::AbstractVector)

Compute the angle between two vectors `a` and `b` in radians.

Supports both numeric and symbolic vectors.

# Examples
```julia
julia> Angle([1, 0], [0, 1])
1.5707963267948966  # π/2
```
"""
function Angle(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return acos(clamp(Dot(a, b) / (Norm(a)*Norm(b)), -1.0, 1.0))
end

function Angle(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    return acos(Dot(a, b) / (Norm(a)*Norm(b)))
end

"""
    Projection(a::AbstractVector, b::AbstractVector)

Project vector `a` onto vector `b`.

Returns the vector projection of `a` onto `b`.

# Examples
```julia
julia> Projection([3, 4], [1, 0])
2-element Vector{Float64}:
 3.0
 0.0
```
"""
function Projection(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    if Norm(b) == 0
        error("Cannot project onto the zero vector")
    end
    return Dot(a, b)/Dot(b, b)*b
end

function Projection(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    return Dot(a, b)/Dot(b, b)*b
end

"""
    ParametricLine(x::AbstractVector, v::AbstractVector)

Generate the parametric equation for a line through point `x` parallel to direction vector `v` in 3D.

Returns a parametric expression in terms of parameter `t`.

# Examples
```julia
julia> ParametricLine([1, 2, 3], [1, 0, 0])
3-element Vector:
 t + 1
 2
 3
```
"""
function ParametricLine(x::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    @syms t
    if length(x) != 3 || length(v) != 3
        error("Enter the 3 dimensional point coordinate and the vector")
    end
    return x .+ t .* v 
end

function ParametricLine(x::AbstractVector{<:Sym}, v::AbstractVector{<:Sym})
    @syms t
    if length(x) != 3 || length(v) != 3
        error("Enter the 3 dimensional point coordinate and the vector")
    end
    return x .+ t .* v 
end

"""
    PlaneEquation(p::AbstractVector, n::AbstractVector)

Generate the equation for a plane given a point `p` on the plane and normal vector `n`.

Returns the plane equation in the form: n₁x + n₂y + n₃z - d = 0

# Examples
```julia
julia> PlaneEquation([0, 0, 1], [0, 0, 1])
z - 1
```
"""
function PlaneEquation(p::AbstractVector{<:Real}, n::AbstractVector{<:Real})
    @syms x y z
    if length(p) != 3 || length(n) != 3
        error("Enter the 3 dimensional point on the plane coordinate and the normal vector")
    end
    return n[1]*x + n[2]*y + n[3]*z - (n[1]*p[1] + n[2]*p[2] + n[3]*p[3])
end

function PlaneEquation(p::AbstractVector{<:Sym}, n::AbstractVector{<:Sym})
    @syms x y z
    if length(p) != 3 || length(n) != 3
        error("Enter the 3 dimensional point on the plane coordinate and the normal vector")
    end
    return n[1]*x + n[2]*y + n[3]*z - (n[1]*p[1] + n[2]*p[2] + n[3]*p[3])
end

"""
    ArcLength(curve::Vector, t::Sym, a, b; symbolic=true)

Compute the arc length of a parametrized curve from parameter value `a` to `b`.

# Arguments
- `curve::Vector`: Parametric curve components [x(t), y(t), z(t)]
- `t::Sym`: Parameter symbol
- `a`: Starting parameter value
- `b`: Ending parameter value
- `symbolic::Bool`: If true, attempts symbolic integration; otherwise uses numerical integration

# Returns
Arc length value (symbolic or numeric)
"""
function ArcLength(curve::Vector, t::Sym, a, b; symbolic=true)
    
    # Compute derivatives
    derivatives = [diff(component, t) for component in curve]
    
    # Compute integrand: sqrt(∑(dx_i/dt)²)
    integrand = sqrt(sum(d^2 for d in derivatives))
    integrand_simplified = simplify(integrand)
    
    println("Arc length integral: ∫ from $a to $b of ", integrand_simplified, " d$t")
    
    if symbolic
        try
            # Attempt symbolic integration
            result = integrate(integrand_simplified, (t, a, b))
            
            if !(result isa SymPy.Integral)  # If integration was successful
                println("Symbolic result: ", result)
                return result
            else
                println("No closed-form solution found, falling back to numerical")
            end
        catch e
            println("Symbolic integration failed: ", e)
        end
    end
    
    # Numerical integration fallback
    println("Using numerical integration")
    integrand_func = lambdify(integrand_simplified, [t])
    numeric_result, error = quadgk(integrand_func, a, b, rtol=1e-8)
    println("Numerical result: ", numeric_result)
    return numeric_result
end

"""
    ArcLengthParametrization(curve::Vector, s::Sym; symbolic=true)

Compute the arc length parametrization of a curve.

# Arguments
- `curve::Vector`: Parametric curve components
- `s::Sym`: Parameter symbol for the curve
- `symbolic::Bool`: If true, attempts symbolic integration

# Returns
Arc length function s(t)
"""
function ArcLengthParametrization(curve::Vector, s::Sym; symbolic=true)
    @syms t 
    
    # Compute derivatives
    derivatives = [diff(component, s) for component in curve]
    
    # Compute integrand: sqrt(∑(dx_i/dt)²)
    integrand = sqrt(sum(d^2 for d in derivatives))
    integrand_simplified = simplify(integrand)
    
    if symbolic
        result = integrate(integrand_simplified, (s, 0, t))
        
        # Check if result contains "Integral" (unevaluated integral)
        if occursin("Integral", string(result))
            error("No closed-form solution found for arc length integral")
        end
        
        println("Symbolic result: ", result)
        return result
    end
end

"""
    Tangent(curve::Vector, t::Sym, t_val=nothing)

Compute the unit tangent vector T(t) = r'(t)/|r'(t)| for a parametric curve.

# Arguments
- `curve::Vector`: Parametric curve components
- `t::Sym`: Parameter symbol
- `t_val`: Optional numeric value to evaluate at

# Returns
Unit tangent vector
"""
function Tangent(curve::Vector, t::Sym, t_val=nothing)

    # Compute derivative of each component
    derivative = [diff(c, t) for c in curve]

    # Compute magnitude (norm) of derivative
    integrand = sqrt(sum(d^2 for d in derivative))

    # Compute unit tangent vector
    tangent = [d / integrand for d in derivative]

    # If t_val is provided, evaluate numerically
    if t_val !== nothing
        tangent = [simplify(subs(d, t => t_val)) for d in tangent]
    end

    return tangent
end

"""
    Curvature(curve::Vector, t::Sym, t_val=nothing)

Compute the curvature κ = |r'(t) × r''(t)|/|r'(t)|³ of a 3D parametric curve.

# Arguments
- `curve::Vector`: 3D parametric curve components
- `t::Sym`: Parameter symbol
- `t_val`: Optional numeric value to evaluate at

# Returns
Curvature value
"""
function Curvature(curve::Vector, t::Sym, t_val=nothing)
    @assert length(curve) == 3

    #compute the first derivative
    first_derivative = [diff(c,t) for c in curve]
    #compute the second derivative
    second_derivative = [diff(c1, t) for c1 in first_derivative]
    cross_product = Cross(first_derivative, second_derivative)
    curvature = Norm(cross_product)/(Norm(first_derivative))^3

    if t_val !== nothing
        curvature = subs(curvature, t => t_val)
    end

    return curvature
end

"""
    Normal(curve::Vector, t::Sym, t_val=nothing)

Compute the principal normal vector N(t) = T'(t)/|T'(t)| for a 3D parametric curve.

# Arguments
- `curve::Vector`: 3D parametric curve components
- `t::Sym`: Parameter symbol
- `t_val`: Optional numeric value to evaluate at

# Returns
Principal normal vector
"""
function Normal(curve::Vector, t::Sym, t_val=nothing) 

    @assert length(curve) == 3

    T_vec = Tangent(curve,t)

    #compute the derivative of tangent 
    derivative = [diff(Ti,t) for Ti in T_vec]
    #Normalize
    N = [n/ Norm(derivative) for n in derivative]

    #substitute numeric value if provided
    if t_val !== nothing
        N = [subs(Ni, t => t_val) for Ni in N]
    end
    return N 
end

"""
    Binormal(curve::Vector, t::Sym, t_val=nothing)

Compute the binormal vector B(t) = T(t) × N(t) for a 3D parametric curve.

# Arguments
- `curve::Vector`: 3D parametric curve components
- `t::Sym`: Parameter symbol
- `t_val`: Optional numeric value to evaluate at

# Returns
Binormal vector
"""
function Binormal(curve::Vector, t::Sym, t_val=nothing)
    @assert length(curve) == 3
    B_vec = simplify(Cross(Tangent(curve,t), Normal(curve, t)))

    if t_val !== nothing
        B_vec = [simplify(subs(Bi, t => t_val)) for Bi in B_vec]
    end
    return B_vec 
end

"""
    Torsion(curve::Vector, t::Sym, t_val=nothing)

Compute the torsion τ = (r'(t) × r''(t)) · r'''(t) / |r'(t) × r''(t)|² of a 3D parametric curve.

# Arguments
- `curve::Vector`: 3D parametric curve components
- `t::Sym`: Parameter symbol
- `t_val`: Optional numeric value to evaluate at

# Returns
Torsion value
"""
function Torsion(curve::Vector, t::Sym, t_val=nothing)
    @assert length(curve) == 3
    
    #compute the first derivative
    first_derivative = [diff(c,t) for c in curve]
    #compute the second derivative
    second_derivative = [diff(c1, t) for c1 in first_derivative]
    #compute the third derivative
    third_derivative = [diff(c2, t) for c2 in second_derivative]
    #cross product of first and second derivative
    cross_product = Cross(first_derivative, second_derivative)
    if Norm(cross_product) == 0
        error("Torsion is undefined at this point as the cross product of first and second derivative is zero.")
    else
        torsion = Dot(cross_product, third_derivative)/(Norm(cross_product))^2

        if t_val !== nothing
            torsion = subs(torsion, t => t_val)
        end
    end

    return torsion
end

"""
    FrenetSerret(curve::Vector, t::Sym, t_val=nothing)

Compute the Frenet-Serret frame (Tangent, Normal, Binormal).

#Arguments 
- 'curve::Vector': 3D Parametric curve components
- 't::Sym': Parameter symbol
- 't_val': Optional numeric value to evaluate at

#Returns 
Tuple of (Tangent, Normal, Binormal) vectors. 
"""

function FrenetSerret(curve::Vector, t::Sym, t_val=nothing)
    T = Tangent(curve, t, t_val)
    N = Normal(curve, t, t_val)
    B = Binormal(curve, t, t_val)

    return (T, N, B)
end

end # module VectorUtils