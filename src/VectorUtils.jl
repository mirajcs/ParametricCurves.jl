module VectorUtils

greet() = print("Hello Geometrician")

export Norm, Normalize, Dot, Cross, Angle, Projection, ParametricLine, PlaneEquation, ArcLength, ArcLengthParametrization, Tangent, Curvature, Normal, Binormal


using LinearAlgebra, SymPy, QuadGK

# Compute the Euclidean norm of vector 'v'

# numeric version
Norm(v::AbstractVector{<:Number}) = norm(v)

# symbolic version
function Norm(v::AbstractVector{<:Sym})
    return simplify(sqrt(sum(x->x^2, v)))
end

# Return the unit vector in the direction of `v`

# Numerical 
function Normalize(v::AbstractVector{<:Real})
    n = Norm(v)
    if n == 0
        error("Zero divisor")
    end 
    normalize = [v[i] / n for i in 1:length(v)]
    return normalize
end 

# symbolic
function Normalize(v::AbstractVector{<:Sym})
    n = Norm(v)
    normalize = [v[i] / n for i in 1:length(v)]
    return simplify(normalize)
end

# compute the dot product of two vectors 

Dot(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) = dot(a, b)

function Dot(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    @assert length(a) == length(b)
    return sum(a[i]*b[i] for i in 1:length(a))
end

# compute the cross product 

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

# compute the angle between two vectors, the unit is radian 
function Angle(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return acos(clamp(Dot(a, b) / (Norm(a)*Norm(b)), -1.0, 1.0))
end

function Angle(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    return acos(Dot(a, b) / (Norm(a)*Norm(b)))
end

# Project vector a onto b 
function Projection(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    if Norm(b) == 0
        error("Cannot project onto the zero vector")
    end
    return Dot(a, b)/Dot(b, b)*b
end

function Projection(a::AbstractVector{<:Sym}, b::AbstractVector{<:Sym})
    return Dot(a, b)/Dot(b, b)*b
end

# parametric equation for the line through the point x and parallel to the direction of the vector v in 3D
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
If we have the point of a plane and normal vector we have the following procedure. 
"""

# Equation for a plane 
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
Arc length for parametrized curve
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
Arc Length prametrization 
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
The tangent vector for the curve is given by r'(t)/|r'(t)|
"""

function Tangent(curve::Vector, t::Sym,  t_val=nothing)

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
Curvature of a curve is given by k  = |r'(t) x r''(t)|/|r'(t)|^3
"""

function Curvature(curve::Vector, t::Sym ,t_val = nothing)
    @assert length(curve) == 3

    #compute the first deirvative
    first_derivative = [diff(c,t) for c in curve]
    #compute the second derivative
    second_derivative = [diff(c1, t) for c1 in first_derivative]
    cross_product = Cross(first_derivative, second_derivative)
    curvature = Norm(cross_product)/(Norm(first_derivative))^3

    if t_val !==nothing
        curvature = subs(curvature, t => t_val)
    end

    return curvature

end

"""
Normal is given by N(t) = T'(t)/|T'(t)|
"""


function Normal(curve::Vector, t::Sym, t_val = nothing) 

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
Binormal vector is given by B(t) = T(t) x N(T)
"""

function Binormal(curve::Vector, t::Sym, t_val = nothing)
    @assert length(curve) == 3
    B_vec = simplify(Cross(Tangent(curve,t) , Normal(curve, t)))

    if t_val !== nothing
        B_vec = [simplify(subs(Bi, t => t_val)) for Bi in B_vec]
    end
    return B_vec 
end




end # module VectorUtils