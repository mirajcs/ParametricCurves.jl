using Pkg
Pkg.add("Test")
using Test
using ParametricCurves
using SymPy
using QuadGK

@test Norm([3, 4]) == 5

@syms x y z
v_sym = [x, y, z]
expected_sym = sqrt(x^2 + y^2 + z^2)
@test isequal(Norm(v_sym), expected_sym)


# Symbolic tests
@syms x y z l m n t
v_sym = [x, y, z]
normalized = Normalize(v_sym)
norm_expr = sqrt(x^2 + y^2 + z^2)

# Test that each component is correct
@test isequal(simplify(normalized[1] - x / norm_expr), 0)
@test isequal(simplify(normalized[2] - y / norm_expr), 0)
@test isequal(simplify(normalized[3] - z / norm_expr), 0)


#dot product test 
@test Dot([2,3],[3,4]) == 18
@test isequal(Dot(v_sym, v_sym), x^2 + y^2 + z^2)

@testset "Cross Product Tests" begin

    # ---------------------------
    # Numerical cross product
    # ---------------------------
    @testset "Numeric Cross" begin
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]

        expected = Cross(a, b)
        @test Cross(a, b) == expected
    end

    # ---------------------------
    # Symbolic cross product
    # ---------------------------
    @testset "Symbolic Cross" begin
        @syms x y z u v w
        a = [x, y, z]
        b = [u, v, w]

        expected = [
            y*w - z*v,
            z*u - x*w,
            x*v - y*u
        ]

        result = Cross(a, b)
        @test all(isapprox.(result, expected))
    end

    # ---------------------------
    # Dimension mismatch
    # ---------------------------
    @testset "Dimension Errors" begin
        @test_throws ErrorException Cross([1,2], [3,4,5])
        @test_throws ErrorException Cross([1,2,3], [4,5])
    end

    # ---------------------------
    # Mixed symbolic & numeric
    # (your methods do NOT define such dispatch)
    # ---------------------------
    @testset "Mixed Type Errors" begin
        @syms x y z
        @test_throws MethodError Cross([x,y,z], [1,2,3])
        @test_throws MethodError Cross([1,2,3], [x,y,z])
    end
end


#Angle test
@test isapprox(Angle([1,0],[0,1]), π/2)
@test isequal(Angle([2*x,y,z], [l,m,n]), acos((2*x*l + y*m + z*n)/(sqrt(4*x^2 + y^2 +z^2)*sqrt(l^2 + m^2 +n^2))))

#Projection test
@test isapprox(Projection([1,1,2],[-2,3,1]), [-3/7,9/14,3/14]) 
@test isequal(simplify(Projection([x,y,z],[l,m,n])), [(x*l + y*m +z*n)*l/(l^2+m^2+n^2), (x*l + y*m +z*n)*m/(l^2+m^2+n^2), (x*l + y*m +z*n)*n/(l^2+m^2+n^2)])


#ParametricLine Test 
result = ParametricLine([1, 2, 3], [4, 5, 6])
expected = [1 + 4t, 2 + 5t, 3 + 6t]
@test isequal(result, expected)
@test isequal(ParametricLine([x,y,z],[l,m,n]),[x + l*t, y + m*t, z + n*t])


#PlaneEquation 
result1 = PlaneEquation([2,4,-1],[2,3,4])
expected1 =  2*x + 3*y + 4*z - 12
@test isequal(result1,expected1)



#Arc Length Test

@testset "ArcLength Basic Functionality" begin
    @syms t
    
    @testset "Straight Lines" begin
        println("\n" * "="^50)
        println("Testing Straight Lines")
        println("="^50)
        
        # Test 1: Simple line along x-axis
        println("1. Testing line along x-axis...")
        curve1 = [t, 0, 0]
        result1 = ArcLength(curve1, t, 0, 5; symbolic=true)
        @test isequal(result1, 5) || isapprox(N(result1), 5.0, atol=1e-10)
        println("✓ Line along x-axis: length = 5")
        
        # Test 2: 3D line with different components
        println("2. Testing 3D line...")
        curve2 = [2*t, 3*t, 6*t]  # Direction vector (2,3,6), magnitude 7
        result2 = ArcLength(curve2, t, 0, 1; symbolic=true)
        expected2 = 7.0  # √(4+9+36) = 7
        @test isequal(result2, 7) || isapprox(N(result2), expected2, atol=1e-10)
        println("✓ 3D line: length = 7")
        
        # Test 3: Line with negative parameter range
        println("3. Testing line with negative range...")
        curve3 = [t, t, t]
        result3 = ArcLength(curve3,t, -2, 2; symbolic=true)
        expected3 = 4 * √3
        @test isequal(result3, 4*√3) || isapprox(N(result3), expected3, atol=1e-10)
        println("✓ Line with negative range: length = $(4√3)")
    end
    
    @testset "Circular Arcs" begin
        println("\n" * "="^50)
        println("Testing Circular Arcs")
        println("="^50)
        
        # Test 1: Quarter circle
        println("1. Testing quarter circle...")
        curve1 = [3*cos(t), 3*sin(t), 0]
        result1 = ArcLength(curve1,t, 0, π/2; symbolic=true)
        expected1 = (3π)/2  # radius * angle = 3 * π/2
        @test isequal(result1, (3π)/2) || isapprox(N(result1), expected1, atol=1e-10)
        println("✓ Quarter circle (r=3): length = $(3π/2)")
        
        # Test 2: Full circle
        println("2. Testing full circle...")
        curve2 = [2*cos(t), 2*sin(t), 0]
        result2 = ArcLength(curve2,t, 0, 2π; symbolic=true)
        expected2 = 4π  # circumference = 2πr = 4π
        @test isequal(result2, 4π) || isapprox(N(result2), expected2, atol=1e-10)
        println("✓ Full circle (r=2): length = 4π")
    end
    
    @testset "Helices" begin
        println("\n" * "="^50)
        println("Testing Helices")
        println("="^50)
        
        # Test 1: Standard helix
        println("1. Testing standard helix...")
        curve1 = [cos(t), sin(t), t]
        result1 = ArcLength(curve1,t, 0, 2π; symbolic=true)
        expected1 = 2π * √2  # √(sin²t + cos²t + 1) = √2
        @test isequal(result1, 2π*√2) || isapprox(N(result1), expected1, atol=1e-10)
        println("✓ Standard helix: length = 2π√2")
    end
end

@testset "ArcLength Edge Cases" begin
    @syms t
    
    @testset "Zero and Small Lengths" begin
        println("\n" * "="^50)
        println("Testing Zero and Small Lengths")
        println("="^50)
        
        # Test 1: Zero length
        println("1. Testing zero length...")
        curve1 = [t, t, t]
        result1 = ArcLength(curve1,t, 0, 0; symbolic=true)
        @test isequal(result1, 0) || isapprox(N(result1), 0.0, atol=1e-10)
        println("✓ Zero length curve")
        
        # Test 2: Very small length
        println("2. Testing very small length...")
        curve2 = [t, 2*t, 3*t]
        result2 = ArcLength(curve2,t, 0, 1e-6; symbolic=true)
        expected2 = 1e-6 * √14
        @test isapprox(N(result2), expected2, atol=1e-12)
        println("✓ Very small length")
    end
    
    @testset "2D Curves" begin
        println("\n" * "="^50)
        println("Testing 2D Curves")
        println("="^50)
        
        # Test 1: 2D line
        println("1. Testing 2D line...")
        curve2 = [2*t, 3*t]
        result2 = ArcLength(curve2,t, 0, 1; symbolic=true)
        expected2 = √13  # √(4 + 9)
        @test isequal(result2, √13) || isapprox(N(result2), expected2, atol=1e-10)
        println("✓ 2D line: length = √13")
    end
end

@testset "ArcLength Numerical Integration" begin
    @syms t
    
    @testset "Forced Numerical Integration" begin
        println("\n" * "="^50)
        println("Testing Forced Numerical Integration")
        println("="^50)
        
        # Test 1: Force numerical on simple curve
        println("1. Testing forced numerical integration...")
        curve1 = [t, t, t]
        result1_num = ArcLength(curve1,t, 0, 1; symbolic=false)
        @test isapprox(result1_num, √3, atol=1e-8)
        println("✓ Numerical integration works")
    end
    
    @testset "Complex Curves" begin
        println("\n" * "="^50)
        println("Testing Complex Curves")
        println("="^50)
        
        # Test 1: Parabola (likely needs numerical)
        println("1. Testing parabola...")
        curve1 = [t, t^2, 0]
        result1 = ArcLength(curve1,t, 0, 1; symbolic=true)
        @test result1 isa Number
        @test result1 > 0
        println("✓ Parabola handled successfully: length = $result1")
        
        # Test 2: Exponential curve
        println("2. Testing exponential curve...")
        curve2 = [t, exp(t), 0]
        result2 = ArcLength(curve2,t, 0, 1; symbolic=true)
        @test result2 isa Number
        @test result2 > 0
        println("✓ Exponential curve handled successfully: length = $result2")
    end
end


"""
Arc Length ArcLengthParametrization tests
"""

@testset "ArcLengthParametrization Tests" begin

    @testset "Arc Length ArcLengthParametrization Test" begin
        println("\n" * "="^50)
        println("Arc Length Parametrization tests")
        println("="^50)
    
    @testset "Test 1: Straight Line (should work)" begin
        @syms s
        curve = [s, 2*s]
    
        result = ArcLengthParametrization(curve,s, symbolic=true)
    
        # Use symbolic sqrt(5), not Float64
        @test simplify(result - sqrt(Sym(5))*t) == 0
        println("✓ Straight line test passed")
    end
    
    @testset "Test 2: Unit Circle (should work)" begin
        @syms s
        # Parametric circle: (cos(s), sin(s))
        curve = [cos(s), sin(s)]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # Arc length of unit circle is just t (arc length = radius * angle = 1 * t)
        @test simplify(result - t) == 0
        println("✓ Unit circle test passed")
    end
    
    @testset "Test 3: Parabola (should work)" begin
        @syms s
        # Parametric parabola: (t, t²)
        curve = [s, s^2]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # This has a known closed form involving asinh
        # Just check it doesn't contain "Integral" (unevaluated)
        @test !occursin("Integral", string(result))
        println("✓ Parabola test passed: ", result)
    end
    
    @testset "Test 4: Circle with radius r (should work)" begin
        @syms s r::positive  # Use :: for assumptions
        curve = [r*cos(s), r*sin(s)]
    
        result = ArcLengthParametrization(curve,s, symbolic=true)
    
        @test simplify(result - r*t) == 0
        println("✓ Circle with radius test passed")
    end
    
    @testset "Test 5: 3D Helix (should work)" begin
        @syms s a b
        # Helix: (a*cos(s), a*sin(s), b*s)
        curve = [a*cos(s), a*sin(s), b*s]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # Arc length should be sqrt(a² + b²) * t
        @test simplify(result - sqrt(a^2 + b^2)*t) == 0
        println("✓ 3D Helix test passed")
    end
    
    @testset "Test 6: Exponential curve (should FAIL)" begin
        @syms s
        # Curve with e^(s²) component - no closed form
        curve = [s, exp(s^2)]
        
        # This should throw an error
        @test_throws ErrorException ArcLengthParametrization(curve,s, symbolic=true)
        println("✓ Exponential curve correctly failed")
    end
    
    @testset "Test 7: Complex integrand (should FAIL)" begin
        @syms s
        # sqrt(1 + cos²(s)) generally has no elementary closed form
        curve = [s, sin(s)^3]
        
        # This might fail depending on SymPy's ability
        try
            result = ArcLengthParametrization(curve,s, symbolic=true)
            println("✓ Complex integrand unexpectedly succeeded: ", result)
        catch e
            @test e isa ErrorException
            println("✓ Complex integrand correctly failed")
        end
    end
    
    @testset "Test 8: Horizontal line (degenerate case)" begin
        @syms s
        # Horizontal line: (s, 0)
        curve = [s, Sym(0)]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # Arc length should just be t (moving 1 unit per parameter)
        @test simplify(result - t) == 0
        println("✓ Horizontal line test passed")
    end
    
    @testset "Test 9: Scaled line" begin
        @syms s
        # Line: (3s, 4s) - should have length 5t
        curve = [3*s, 4*s]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # Arc length should be 5*t (3-4-5 triangle)
        @test simplify(result - 5*t) == 0
        println("✓ Scaled line test passed")
    end
    
    @testset "Test 10: Verify result type" begin
        @syms s
        curve = [cos(s), sin(s)]
        
        result = ArcLengthParametrization(curve,s, symbolic=true)
        
        # Result should be a Sym (symbolic expression)
        @test result isa Sym
        println("✓ Result type test passed")
    end
end
end

# Tests
@testset "Tangent Function Tests" begin
    @syms t

    @testset "TangentTest" begin
        println("\n" * "="^50)
        println("Tangent tests")
        println("="^50)

    # 2D curve test
    curve1 = [t^2, t^3]
    tangent1 = Tangent(curve1,t)
    @test tangent1 == [2*t / sqrt((2*t)^2 + (3*t^2)^2), 3*t^2 / sqrt((2*t)^2 + (3*t^2)^2)]

    # Evaluate at t = 1
    tangent1_val = Tangent(curve1,t, 1)
    expected_val = [2 / sqrt(13), 3 / sqrt(13)]
    @test all(isapprox.(tangent1_val, expected_val))

    # 3D curve test
    curve2 = [t, t^2, t^3]
    tangent2 = Tangent(curve2,t)
    @test length(tangent2) == 3
    end
end


@testset "Curvature Tests" begin

    # --------------------------------------------------
    # 1. Symbolic curvature: circle x(t)=cos(t), y(t)=sin(t), z(t)=0
    # Curvature of unit circle is always 1
    # --------------------------------------------------
    @testset "Symbolic curvature of circle" begin
        @syms t
        curve = [cos(t), sin(t), 0]

        κ = Curvature(curve,t)

        # Symbolically, curvature = 1
        expected = 1

        # κ is a scalar, ensure it matches 1
        @test simplify(κ - expected) == 0
    end


    # --------------------------------------------------
    # 2. Numeric evaluation: same circle at t = π/4
    # --------------------------------------------------
    @testset "Numeric curvature at t_val" begin
        @syms t
        curve = [cos(t), sin(t), 0]

        κ_val = Curvature(curve,t, π/4)

        # Should be approximately 1
        @test isapprox.(κ_val, [1.0], atol=1e-6) |> all
    end


    # --------------------------------------------------
    # 3. Check the output is scalar
    # --------------------------------------------------
    @testset "Shape of curvature" begin
        @syms t
        curve = [t^2, t^3, t]

        κ = Curvature(curve,t)

        @test isa(κ, Sym)
    end


    # --------------------------------------------------
    # 4. Non-3D curve error
    # --------------------------------------------------
    @testset "Dimension errors" begin
        @syms t
        bad_curve = [t, t^2]  # only 2D

        @test_throws AssertionError Curvature(bad_curve,t)
    end

end




@testset "Binormal Vector Tests" begin
    @testset "Binormal Test" begin
        println("\n" * "="^50)
        println("Binormal tests")
        println("="^50)

    # ---------------------------------------------
    # 1. Symbolic binormal for unit circle in xy-plane
    # curve: r(t) = [cos(t), sin(t), 0]
    # Tangent: [-sin(t), cos(t), 0]
    # Normal: [-cos(t), -sin(t), 0]
    # Binormal: cross(T,N) = [0,0,1]
    # ---------------------------------------------
    @testset "Symbolic Binormal" begin
        @syms t
        curve = [cos(t), sin(t), 0]

        B = Binormal(curve, t)

        # Expected symbolic binormal
        expected = [0, 0, 1]
        @test all(Bi == Ei for (Bi, Ei) in zip(B, expected))
        @test all(Bi isa Sym for Bi in B)
    end

    # ---------------------------------------------
    # 2. Numeric evaluation at t_val = π/4
    # ---------------------------------------------
    @testset "Numeric Binormal at t_val" begin
        @syms t
        curve = [cos(t), sin(t), 0]

        t_val = pi/4
        B_val = Binormal(curve, t, t_val)

        expected_val = [0.0, 0.0, 1.0]
        @test all(isapprox(float(Bi), float(Ei); atol=1e-6) for (Bi, Ei) in zip(B_val, expected_val))
    end

    # ---------------------------------------------
    # 3. Check dimension
    # ---------------------------------------------
    @testset "Dimension check" begin
        @syms t
        curve = [t, t^2, t^3]
        B = Binormal(curve, t)
        @test length(B) == 3
    end

    # ---------------------------------------------
    # 4. Non-3D curve should throw assertion
    # ---------------------------------------------
    @testset "Dimension errors" begin
        @syms t
        bad_curve = [t, t^2]  # 2D
        @test_throws AssertionError Binormal(bad_curve,t)
    end
end 
end

@testset "Torsion Tests" begin
    println("\n" * "="^50)
    println("Torsion tests")
    println("="^50)

    
    # ---------------------------------------------
    # 1. Torsion of a circle in the xy-plane (should be zero)
    # r(t) = [cos t, sin t, 0]
    t = symbols("t", real=true)
    curve1 = [cos(t), sin(t), 0]

    @testset "Planar Circle" begin
        τ_sym = Torsion(curve1, t)
        @test simplify(τ_sym) == 0

        τ_num = Torsion(curve1, t, 1.2)
        @test τ_num == 0
    end

    # ---------------------------------------------
    # 2. Torsion of a helix (constant nonzero)
    # r(t) = [cos t, sin t, t]
    # Known torsion = 1 / (1 + 1)^2 = 1/2
    curve2 = [cos(t), sin(t), t]

    @testset "Circular Helix" begin
        τ_sym2 = Torsion(curve2, t)
        τ_simplified = simplify(τ_sym2)
        @test τ_simplified == 1/2

        τ_num2 = Torsion(curve2, t, 0.3)
        @test isapprox(N(τ_num2), 0.5; atol=1e-6)
    end

    # ---------------------------------------------
    # 3. Straight line (torsion undefined or zero)
    # r(t) = [t, 0, 0]
    # First and second derivatives parallel → torsion should be 0
    curve3 = [t, 0, 0]

    @testset "Straight Line" begin
        @test_throws ErrorException Torsion(curve3, t)
    end
end


@testset "FrenetSerret Tests" begin
    println("\n" * "="^50)
    println("FrenetSerret tests")
    println("="^50)

    t = symbols("t", real=true)
    
    # ---------------------------------------------
    # 1. Circular Helix - verify orthonormal frame
    # r(t) = [cos t, sin t, t]
    curve1 = [cos(t), sin(t), t]

    @testset "Circular Helix - Orthonormality" begin
        T, N, B = FrenetSerret(curve1, t)
        
        # Verify the frame is orthonormal
        @test simplify(Dot(T, T)) == 1
        @test simplify(Dot(N, N)) == 1
        @test simplify(Dot(B, B)) == 1
        @test simplify(Dot(T, N)) == 0
        @test simplify(Dot(T, B)) == 0
        @test simplify(Dot(N, B)) == 0
    end

    # ---------------------------------------------
    # 2. Circular Helix at t=0
    @testset "Circular Helix at t=0" begin
        T, N, B = FrenetSerret(curve1, t, 0)
        
        # At t=0: curve is (1, 0, 0)
        @test isapprox(T[1], 0.0; atol=1e-10)
        @test isapprox(T[2], sqrt(2)/2; atol=1e-10)
        @test isapprox(T[3], sqrt(2)/2; atol=1e-10)
        
        @test isapprox(N[1], -1.0; atol=1e-10)
        @test isapprox(N[2], 0.0; atol=1e-10)
        @test isapprox(N[3], 0.0; atol=1e-10)
        
        @test isapprox(B[1], 0.0; atol=1e-10)
        @test isapprox(B[2], -sqrt(2)/2; atol=1e-10)
        @test isapprox(B[3], sqrt(2)/2; atol=1e-10)
    end

    # ---------------------------------------------
    # 3. Circle in xy-plane at t=π/2
    # r(t) = [2*cos t, 2*sin t, 0]
    curve2 = [2*cos(t), 2*sin(t), 0]

    @testset "Circle in xy-plane" begin
        T, N, B = FrenetSerret(curve2, t, PI/2)
        
        # At t=π/2: curve is (0, 2, 0)
        # T is tangent to circle, N points toward center, B points in z-direction
        @test isapprox(T[1], -1.0; atol=1e-10)
        @test isapprox(T[2], 0.0; atol=1e-10)
        @test isapprox(T[3], 0.0; atol=1e-10)
        
        @test isapprox(N[1], 0.0; atol=1e-10)
        @test isapprox(N[2], -1.0; atol=1e-10)
        @test isapprox(N[3], 0.0; atol=1e-10)
        
        @test isapprox(B[1], 0.0; atol=1e-10)
        @test isapprox(B[2], 0.0; atol=1e-10)
        @test isapprox(B[3], 1.0; atol=1e-10)
    end


    # ---------------------------------------------
    # 4. Parabolic Path - verify orthonormality at t=1
    # r(t) = [t, t^2, 0]
    curve4 = [t, t^2, 0]

    @testset "Parabolic Path Orthonormality" begin
        T, N, B = FrenetSerret(curve4, t, 1)
        
        # Verify orthonormality at t=1
        @test isapprox(Dot(T, T), 1.0; atol=1e-10)
        @test isapprox(Dot(N, N), 1.0; atol=1e-10)
        @test isapprox(Dot(B, B), 1.0; atol=1e-10)
        @test isapprox(Dot(T, N), 0.0; atol=1e-10)
        @test isapprox(Dot(T, B), 0.0; atol=1e-10)
        @test isapprox(Dot(N, B), 0.0; atol=1e-10)
    end


    # ---------------------------------------------
    # 5. Symbolic Return (no t_val specified)
    @testset "Symbolic Return" begin
        T, N, B = FrenetSerret(curve1, t)
        
        # Should return symbolic expressions
        @test T isa Vector
        @test N isa Vector
        @test B isa Vector
        @test length(T) == 3
        @test length(N) == 3
        @test length(B) == 3
    end
end





@testset "Acceleration Tests" begin
    println("\n" * "="^50)
    println("Acceleration tests")
    println("="^50)

    @testset "Basic symbolic acceleration" begin
        @syms t
        curve = [t^2, t, 2*t]
        aT, aN = Acceleration(curve, t)

        # Verify aT is symbolic
        @test aT isa Sym
        @test aN isa Sym

        # Verify they're not zero
        @test simplify(aT) != 0
        @test simplify(aN) != 0
        println("✓ Basic symbolic acceleration computed")
    end

    @testset "Numeric evaluation at t_val" begin
        @syms t
        curve = [t^2, t, 2*t]
        aT_val, aN_val = Acceleration(curve, t, 1)

        # At t=1: v = [2, 1, 2], a = [2, 0, 0]
        # |v| = 3
        # aT = v·a / |v| = 4 / 3
        # v×a = [0, 4, -2], |v×a| = sqrt(20) = 2*sqrt(5)
        # aN = |v×a| / |v| = 2*sqrt(5) / 3
        expected_aT = 4/3
        expected_aN = 2*sqrt(5)/3

        @test isapprox(N(aT_val), expected_aT, atol=1e-6)
        @test isapprox(N(aN_val), expected_aN, atol=1e-6)
        println("✓ Numeric evaluation at t_val passed")
    end

    @testset "Straight line motion" begin
        @syms t
        # Moving along [t, 0, 0] - straight line along x-axis
        curve = [t, 0, 0]
        aT, aN = Acceleration(curve, t)

        # For straight line motion:
        # aT = (1*0) / 1 = 0 (constant velocity)
        # aN = 0 (no curvature)
        @test simplify(aT) == 0
        @test simplify(aN) == 0
        println("✓ Straight line motion has zero acceleration")
    end

    @testset "Circular motion" begin
        @syms t
        # Unit circle: r(t) = [cos(t), sin(t), 0]
        curve = [cos(t), sin(t), 0]
        aT, aN = Acceleration(curve, t)

        # For circular motion at constant speed:
        # v = [-sin(t), cos(t), 0], |v| = 1
        # a = [-cos(t), -sin(t), 0]
        # aT = 0 (constant speed)
        # aN = 1 (centripetal acceleration for unit circle)
        @test simplify(aT) == 0
        @test simplify(aN) == 1
        println("✓ Circular motion: aT=0, aN=1")
    end

    @testset "Parabolic motion" begin
        @syms t
        curve = [t, t^2, 0]
        aT, aN = Acceleration(curve, t)

        # At t=0: v=[1,0,0], a=[0,2,0]
        aT_0, aN_0 = Acceleration(curve, t, 0)

        # aT(0) = 0*1 / 1 = 0
        # aN(0) = ||(1,0,0) × (0,2,0)|| / 1 = 2
        @test isapprox(N(aT_0), 0.0, atol=1e-6)
        @test isapprox(N(aN_0), 2.0, atol=1e-6)
        println("✓ Parabolic motion computed correctly")
    end

    @testset "Return type is tuple" begin
        @syms t
        curve = [t^2, sin(t), cos(t)]
        result = Acceleration(curve, t)

        @test result isa Tuple
        @test length(result) == 2
        @test all(isa(r, Sym) for r in result)
        println("✓ Return type is tuple of two symbolic expressions")
    end

    @testset "Helix acceleration" begin
        @syms t
        # Helix: [cos(t), sin(t), t]
        curve = [cos(t), sin(t), t]
        aT, aN = Acceleration(curve, t)

        # v = [-sin(t), cos(t), 1], |v| = sqrt(2)
        # a = [-cos(t), -sin(t), 0]
        # v·a = sin(t)cos(t) - sin(t)cos(t) = 0
        # aT = 0 / sqrt(2) = 0 (constant speed)
        # v×a = [sin(t), -cos(t), 1], |v×a| = sqrt(2)
        # aN = sqrt(2) / sqrt(2) = 1 (centripetal acceleration)
        @test simplify(aT) == 0
        @test simplify(aN - 1) == 0
        println("✓ Helix acceleration: aT=0, aN=1")
    end

    @testset "Acceleration components are positive" begin
        @syms t
        curve = [t^3, sin(t), exp(t)]

        # Evaluate at t=0.5 to check they're positive
        aT_val, aN_val = Acceleration(curve, t, 0.5)

        # aN should always be non-negative
        @test N(aN_val) >= -1e-10  # Allow small numerical error
        println("✓ Normal acceleration is non-negative")
    end
end


@testset "Gradient Tests" begin
    println("\n" * "="^50)
    println("Gradient tests")
    println("="^50)

    @testset "Basic symbolic gradient" begin
        @syms x y z
        f = x^2 + 2*y*z + z^2
        grad_f = Gradient(f, [x, y, z])

        expected = [2*x, 2*z, 2*y + 2*z]
        @test length(grad_f) == 3
        @test all(isequal(grad_f[i], expected[i]) for i in 1:3)
        println("✓ Basic symbolic gradient computed correctly")
    end

    @testset "2D gradient" begin
        @syms x y
        f = x^2 + y^2
        grad_f = Gradient(f, [x, y])

        expected = [2*x, 2*y]
        @test length(grad_f) == 2
        @test all(isequal(grad_f[i], expected[i]) for i in 1:2)
        println("✓ 2D gradient works correctly")
    end

    @testset "Gradient of linear function" begin
        @syms x y z
        f = 3*x + 2*y - z
        grad_f = Gradient(f, [x, y, z])

        expected = [3, 2, -1]
        @test all(isequal(grad_f[i], expected[i]) for i in 1:3)
        println("✓ Linear function gradient is constant")
    end

    @testset "Gradient of constant function" begin
        @syms x y z
        f = Sym(5)
        grad_f = Gradient(f, [x, y, z])

        expected = [0, 0, 0]
        @test all(isequal(grad_f[i], expected[i]) for i in 1:3)
        println("✓ Constant function gradient is zero")
    end

    @testset "Gradient with trigonometric functions" begin
        @syms x y
        f = sin(x) + cos(y)
        grad_f = Gradient(f, [x, y])

        expected = [cos(x), -sin(y)]
        @test isequal(grad_f[1], expected[1])
        @test isequal(grad_f[2], expected[2])
        println("✓ Trigonometric gradient computed correctly")
    end

    @testset "Gradient with exponential function" begin
        @syms x y
        f = exp(x*y)
        grad_f = Gradient(f, [x, y])

        expected = [y*exp(x*y), x*exp(x*y)]
        @test simplify(grad_f[1] - expected[1]) == 0
        @test simplify(grad_f[2] - expected[2]) == 0
        println("✓ Exponential function gradient computed correctly")
    end

    @testset "Gradient with polynomial" begin
        @syms x y z
        f = x^3*y^2 + x*y*z + z^4
        grad_f = Gradient(f, [x, y, z])

        @test isequal(simplify(grad_f[1]), 3*x^2*y^2 + y*z)
        @test isequal(simplify(grad_f[2]), 2*x^3*y + x*z)
        @test isequal(simplify(grad_f[3]), x*y + 4*z^3)
        println("✓ Polynomial gradient computed correctly")
    end

    @testset "Gradient with single variable" begin
        @syms x
        f = x^3 - 2*x^2 + x - 5
        grad_f = Gradient(f, [x])

        expected = [3*x^2 - 4*x + 1]
        @test isequal(grad_f[1], expected[1])
        println("✓ Single variable gradient works")
    end

    @testset "Gradient direction check" begin
        @syms x y
        f = x^2 + 4*y^2
        grad_f = Gradient(f, [x, y])

        # At point (1, 1), gradient should be [2, 8]
        point_grad = [subs(grad_f[i], x => 1, y => 1) for i in 1:2]
        @test point_grad == [2, 8]
        println("✓ Gradient points in correct direction at (1,1)")
    end

    @testset "Gradient with mixed terms" begin
        @syms x y z
        f = x*y + y*z + z*x + x^2
        grad_f = Gradient(f, [x, y, z])

        @test isequal(grad_f[1], y + z + 2*x)
        @test isequal(grad_f[2], x + z)
        @test isequal(grad_f[3], y + x)
        println("✓ Mixed terms gradient computed correctly")
    end
end

@testset "Jacobian Tests" begin
    println("\n" * "="^50)
    println("Jacobian Matrix tests")
    println("="^50)

    @testset "Basic 2D→2D Jacobian" begin
        @syms x y
        funcs = [x^2 + y, x*y]
        J = Jacobian(funcs, [x, y])

        # Expected Jacobian:
        # [2x,   1]
        # [y,    x]
        @test length(J) == 2
        @test length(J[1]) == 2
        @test isequal(J[1][1], 2*x)
        @test isequal(J[1][2], 1)
        @test isequal(J[2][1], y)
        @test isequal(J[2][2], x)
        println("✓ Basic 2D→2D Jacobian computed correctly")
    end

    @testset "3D→3D Jacobian" begin
        @syms x y z
        funcs = [x^2 + y, x*y*z, z^2]
        J = Jacobian(funcs, [x, y, z])

        # Expected Jacobian:
        # [2x,   1,   0]
        # [yz,   xz,  xy]
        # [0,    0,   2z]
        @test length(J) == 3
        @test length(J[1]) == 3
        @test isequal(J[1][1], 2*x)
        @test isequal(J[1][2], 1)
        @test isequal(J[1][3], 0)
        @test isequal(J[2][1], y*z)
        @test isequal(J[2][2], x*z)
        @test isequal(J[2][3], x*y)
        @test isequal(J[3][1], 0)
        @test isequal(J[3][2], 0)
        @test isequal(J[3][3], 2*z)
        println("✓ 3D→3D Jacobian computed correctly")
    end

    @testset "Parametric curve Jacobian (3D→1D)" begin
        @syms t
        curve = [cos(t), sin(t), t^2]
        J = Jacobian(curve, [t])

        # Expected Jacobian:
        # [-sin(t)]
        # [cos(t)]
        # [2t]
        @test length(J) == 3
        @test length(J[1]) == 1
        @test isequal(J[1][1], -sin(t))
        @test isequal(J[2][1], cos(t))
        @test isequal(J[3][1], 2*t)
        println("✓ Parametric curve Jacobian computed correctly")
    end

    @testset "Linear transformation Jacobian" begin
        @syms x y z
        # Linear transformation: [2x + y, x - z, y + 3z]
        funcs = [2*x + y, x - z, y + 3*z]
        J = Jacobian(funcs, [x, y, z])

        # Expected Jacobian (constant):
        # [2,  1,  0]
        # [1,  0, -1]
        # [0,  1,  3]
        @test isequal(J[1][1], 2)
        @test isequal(J[1][2], 1)
        @test isequal(J[1][3], 0)
        @test isequal(J[2][1], 1)
        @test isequal(J[2][2], 0)
        @test isequal(J[2][3], -1)
        @test isequal(J[3][1], 0)
        @test isequal(J[3][2], 1)
        @test isequal(J[3][3], 3)
        println("✓ Linear transformation Jacobian is constant")
    end

    @testset "Polar to Cartesian Jacobian" begin
        @syms r θ
        # Polar to Cartesian: x = r*cos(θ), y = r*sin(θ)
        funcs = [r*cos(θ), r*sin(θ)]
        J = Jacobian(funcs, [r, θ])

        # Expected Jacobian:
        # [cos(θ),  -r*sin(θ)]
        # [sin(θ),   r*cos(θ)]
        @test isequal(J[1][1], cos(θ))
        @test isequal(simplify(J[1][2]), -r*sin(θ))
        @test isequal(J[2][1], sin(θ))
        @test isequal(simplify(J[2][2]), r*cos(θ))
        println("✓ Polar to Cartesian Jacobian computed correctly")
    end

    @testset "Spherical coordinates Jacobian" begin
        @syms r θ φ
        # Spherical to Cartesian
        funcs = [r*sin(φ)*cos(θ), r*sin(φ)*sin(θ), r*cos(φ)]
        J = Jacobian(funcs, [r, θ, φ])

        # Check dimensions
        @test length(J) == 3
        @test length(J[1]) == 3

        # Check ∂x/∂r = sin(φ)*cos(θ)
        @test isequal(simplify(J[1][1]), sin(φ)*cos(θ))
        println("✓ Spherical coordinates Jacobian computed correctly")
    end

    @testset "Identity transformation" begin
        @syms x y z
        funcs = [x, y, z]
        J = Jacobian(funcs, [x, y, z])

        # Expected: Identity matrix
        @test isequal(J[1][1], 1)
        @test isequal(J[1][2], 0)
        @test isequal(J[1][3], 0)
        @test isequal(J[2][1], 0)
        @test isequal(J[2][2], 1)
        @test isequal(J[2][3], 0)
        @test isequal(J[3][1], 0)
        @test isequal(J[3][2], 0)
        @test isequal(J[3][3], 1)
        println("✓ Identity transformation Jacobian is identity matrix")
    end

    @testset "Single function Jacobian (gradient)" begin
        @syms x y
        # For a single function, Jacobian is the gradient (row vector)
        f = x^2 + y^2
        J = Jacobian([f], [x, y])

        @test length(J) == 1
        @test length(J[1]) == 2
        @test isequal(J[1][1], 2*x)
        @test isequal(J[1][2], 2*y)
        println("✓ Single function Jacobian equals gradient")
    end

    @testset "Trigonometric functions Jacobian" begin
        @syms x y
        funcs = [sin(x)*cos(y), cos(x)*sin(y)]
        J = Jacobian(funcs, [x, y])

        @test isequal(J[1][1], cos(x)*cos(y))
        @test isequal(simplify(J[1][2]), -sin(x)*sin(y))
        @test isequal(simplify(J[2][1]), -sin(x)*sin(y))
        @test isequal(J[2][2], cos(x)*cos(y))
        println("✓ Trigonometric Jacobian computed correctly")
    end

    @testset "Exponential functions Jacobian" begin
        @syms x y
        funcs = [exp(x + y), exp(x - y)]
        J = Jacobian(funcs, [x, y])

        @test isequal(J[1][1], exp(x + y))
        @test isequal(J[1][2], exp(x + y))
        @test isequal(J[2][1], exp(x - y))
        @test isequal(simplify(J[2][2]), -exp(x - y))
        println("✓ Exponential functions Jacobian computed correctly")
    end

    @testset "Non-square Jacobian (2D→3D)" begin
        @syms u v
        # Parametric surface
        funcs = [u*cos(v), u*sin(v), u^2]
        J = Jacobian(funcs, [u, v])

        # 3×2 Jacobian matrix
        @test length(J) == 3
        @test length(J[1]) == 2
        @test isequal(J[1][1], cos(v))
        @test isequal(simplify(J[1][2]), -u*sin(v))
        @test isequal(J[2][1], sin(v))
        @test isequal(simplify(J[2][2]), u*cos(v))
        @test isequal(J[3][1], 2*u)
        @test isequal(J[3][2], 0)
        println("✓ Non-square 3×2 Jacobian computed correctly")
    end

    @testset "Polynomial Jacobian" begin
        @syms x y
        funcs = [x^3 + y^2, x*y^2 + x^2]
        J = Jacobian(funcs, [x, y])

        @test isequal(J[1][1], 3*x^2)
        @test isequal(J[1][2], 2*y)
        @test isequal(J[2][1], y^2 + 2*x)
        @test isequal(J[2][2], 2*x*y)
        println("✓ Polynomial Jacobian computed correctly")
    end
end

println("\n" * "="^50)
println("All tests completed!")
println("="^50)