using Test
using Random
using LinearAlgebra
using StaticArrays
using GBGKnudsenGas

@testset "GBGKnudsenGas.jl" begin
  @testset "Single collision with a specular sphere" begin
    s = Sphere(SA[0.0, 0.0, 0.0], 1.0)

    # head-on collision, the example from the documentation
    p = collide_with_body(Particle(SA[1.0, 0.0, 0.0], SA[-10.0, 0.0, 0.0], 0.0), s)
    @test p.v ≈ SA[-1.0, 0.0, 0.0]
    @test p.x ≈ SA[-1.0, 0.0, 0.0]
    @test p.t ≈ 9.0

    # the particle passes by the sphere
    @test isnothing(collide_with_body(Particle(SA[1.0, 0.0, 0.0], SA[-10.0, 2.0, 0.0], 0.0), s))

    # the particle moves away from the sphere
    @test isnothing(collide_with_body(Particle(SA[-1.0, 0.0, 0.0], SA[-10.0, 0.0, 0.0], 0.0), s))
  end

  @testset "Epstein drag on a sphere" begin
    # Epstein (1924), Baines et al. (1965):
    # F / ρ = (1 + π/8 α σ_wall / σ_v) 8/3 √(2π) σ_v r² ⟨v₀⟩
    r, σ_v, σ_wall, r_dom = 0.5, 1.5, 3.0, 2.0

    for α in (0.0, 0.5, 1.0)
      Random.seed!(42)

      result = force_and_torque_uniform_flow(
        body           = Sphere(SA[0.0, 0.0, 0.0], r; α, wall_std_v = σ_wall),
        torque_point   = SA[0.0, 0.0, 0.0],
        mean_v         = SA[0.0, 0.0, 0.0],
        std_v          = σ_v,
        domain_radius  = r_dom,
        max_collisions = 1,
        sample_size    = 10^5,
      )

      Epstein_drag = (1 + π/8 * α * σ_wall / σ_v) * 8/3 * √(2π) * σ_v * r^2

      # χ² test with 9 degrees of freedom: the Jacobian must agree with the theory given its covariance;
      # the threshold 45 is exceeded by chance with probability of about 10⁻⁶
      δ = vec(result.dF_v - Epstein_drag * I)
      Σ = reshape(result.ΣdF_v, 9, 9)
      @test δ ⋅ (Σ \ δ) < 45

      # a centered sphere is hit with probability (r / r_dom)²
      @test abs(result.P_hit - (r / r_dom)^2) < 5 * √result.ΣP_hit
    end
  end
end
