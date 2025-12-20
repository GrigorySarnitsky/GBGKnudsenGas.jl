"""
    GBGKnudsenGas

Simulate Knudsen flow of particles past a stationary body, where the particle mean free path is much larger than the body size. In this regime, particles interact with the body but not with each other. The body is assumed to be massive, so particle collisions do not move it.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). The package includes tools to compute the next collision of a single particle with a body; you can use these to build multi-particle, multi-collision simulations. A simple computational domain for modelling a gas in thermal equilibrium flowing past a body is also provided.

All particle–body collisions are assumed to be specular. The package uses a Cartesian coordinate system. Note that `GBGKnudsenGas.jl` is purely kinematic: it works with distances, times, and velocities, but not masses or momenta.

"""
module GBGKnudsenGas

export Particle, Body, Sphere, CylinderSideWall, IsotropicHelicoid, CompoundBody, collide_with_body,
 sample_particle_state_uniform_flow, find_colliding_particle_uniform_flow, force_and_torque_uniform_flow


using StaticArrays
using LinearAlgebra
using Base.Threads


"""
    Particle(v::SVector{3, T}, x::SVector{3, T}, t::T) where {T <: Real}

State of a single point-particle with velocity vector `v` and position vector `x` at time `t`.
"""
struct Particle{T <: Real}
  v::SVector{3, T}
  x::SVector{3, T}
  t::T
end


"""
    Body

Abstract supertype for all bodies that can interact with a `Particle` via `collide_with_body(::Particle, ::Body)`.
"""
abstract type Body end


"""
    collide_with_body(particle::Particle, body::Body) -> Union{Nothing, Particle}

Compute the next collision of `particle` with `body`, assuming straight-line motion between collisions.

Return a new `Particle` located at the collision point, with updated velocity (after the collision) and time, or `nothing` if there is no collision.
"""
collide_with_body(::Particle, ::Body) = error("No collide_with_body defined for this body")


"""
    Sphere <: Body
    Sphere(x::SVector{3, T_x}, r::T_r) where {T_x <: Real, T_r <: Real}

Sphere of radius `r` centered at point `x`.
"""
struct Sphere{T_x <: Real, T_r <: Real} <: Body
  x::SVector{3, T_x}
  r::T_r
end


function collide_with_body(particle::Particle, sphere::Sphere)
  Δx = particle.x - sphere.x
  Δx², v² = Δx ⋅ Δx, particle.v ⋅ particle.v
  b = particle.v ⋅ Δx
  d = b^2 - v² * (Δx² - sphere.r^2)
  
  if b > 0 || d < 0 || Δx² < sphere.r^2 || iszero(v²) return nothing end
  
  time_to_collision = (-b - √d) / v²
  x_at_collision  = particle.x + particle.v * time_to_collision
  Δx_at_collision = x_at_collision - sphere.x
  v_after_collision = particle.v - 2 * Δx_at_collision * (particle.v ⋅ Δx_at_collision) / sphere.r^2
  
  return Particle(v_after_collision, x_at_collision, particle.t + time_to_collision)
end


"""
    CylinderSideWall <: Body
    CylinderSideWall(x_1::SVector{3, T_x}, x_2::SVector{3, T_x}, r::T_r) where {T_x <: Real, T_r <: Real}

Finite circular cylinder defined by its axis endpoints `x_1`, `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. End caps can be added via [`CompoundBody`](@ref), like we do in [`IsotropicHelicoid`](@ref).
"""
struct CylinderSideWall{T_x <: Real, T_r <: Real} <: Body
  x_1::SVector{3, T_x}
  x_2::SVector{3, T_x}
  r  ::T_r
end


function collide_with_body(particle::Particle, cylinder::CylinderSideWall)
  Δx = particle.x - cylinder.x_1
  x_2_1 = cylinder.x_2 - cylinder.x_1
  x_2_1² = x_2_1 ⋅ x_2_1
  Δx_⊥ = Δx - (Δx ⋅ x_2_1) / x_2_1² * x_2_1
  v_⊥  = particle.v - (particle.v ⋅ x_2_1) / x_2_1² * x_2_1
  Δx_⊥², v_⊥² = Δx_⊥ ⋅ Δx_⊥, v_⊥ ⋅ v_⊥
  b = v_⊥ ⋅ Δx_⊥
  d = b^2 - v_⊥² * (Δx_⊥² - cylinder.r^2)
  
  if b > 0 || d < 0 || Δx_⊥² < cylinder.r^2 || iszero(v_⊥²) return nothing end
  
  time_to_collision = (-b - √d) / v_⊥²
  x_at_collision  = particle.x + particle.v * time_to_collision
  Δx_at_collision = x_at_collision - cylinder.x_1
  
  if !(0 < Δx_at_collision ⋅ x_2_1 < x_2_1²) return nothing end

  Δx_⊥_at_collision = Δx_at_collision - (Δx_at_collision ⋅ x_2_1) / x_2_1² * x_2_1
  v_after_collision = particle.v - 2 * Δx_⊥_at_collision * (particle.v ⋅ Δx_⊥_at_collision) / cylinder.r^2
  
  return Particle(v_after_collision, x_at_collision, particle.t + time_to_collision)
end


"""
    CompoundBody <: Body
    CompoundBody{T <: Tuple{Vararg{<:Body}}}

Compound body whose components are stored in the tuple `bodies`. A collision with a `CompoundBody` is defined as the earliest collision with any of its components.

Convenience constructors:

    CompoundBody(bodies::Vararg{<:Body})
    CompoundBody(bodies::AbstractVector{<:Body})
"""
struct CompoundBody{T <: Tuple{Vararg{<:Body}}} <: Body
  bodies::T
end


CompoundBody(bodies::Vararg{<:Body}) = CompoundBody(bodies)
CompoundBody(bodies::AbstractVector{<:Body}) = CompoundBody(Tuple(bodies))


function earlier(a, b)
  a === nothing && return b
  b === nothing && return a
  return a.t <= b.t ? a : b
end


earliest_collision(::Particle, ::Tuple{}) = nothing
function earliest_collision(particle::Particle, bodies::Tuple)
  first_collision = collide_with_body(particle, first(bodies))
  rest_collision  = earliest_collision(particle, Base.tail(bodies))
  return earlier(first_collision, rest_collision)
end


function collide_with_body(particle::Particle, compound_body::CompoundBody)
  return earliest_collision(particle, compound_body.bodies)
end


# """
#     CompoundBodyDynamic <: Body
#     CompoundBodyDynamic(bodies::Array{Body, 1})

# A compound body whose components are stored in the array `bodies`. A collision with a `CompoundBodyDynamic` is defined as the earliest collision with any of its components.

# This representation might be preferable if the number of component bodies is very large, however the default choice should be tuple-based `CompoundBody`. Performance differences between `CompoundBodyDynamic` and `CompoundBody` should be evaluated on a case-by-case basis.
# """
# struct CompoundBodyDynamic <: Body
#   bodies::Array{Body, 1}
# end

# function collide_with_body(particle::Particle, compound_body::CompoundBodyDynamic)
#   nearest_collision_particle = nothing
#   for body in compound_body.bodies
#     collided_particle = collide_with_body(particle, body)
#     nearest_collision_particle = earlier(nearest_collision_particle, collided_particle)
#   end
#   return nearest_collision_particle
# end


"""
    IsotropicHelicoid <: Body
    IsotropicHelicoid(cube_side_length::T_L, r::T_r) where {T_L <: Real, T_r <: Real}

Model of Lord Kelvin's isotropic helicoid, constructed from 12 rods. Each rod is a cylinder of radius `r`, capped with two spheres of the same radius. The helicoid is inscribed in a cube of side length `cube_side_length`, centered at ``(0, 0, 0)``. This cube side length indirectly defines the length of each rod.

The rods that form the helicoid can be accessed via `helicoid.elements`, which is a `CompoundBody`.
"""
struct IsotropicHelicoid{T_L <: Real, T_r <: Real, T_el <: CompoundBody} <: Body
  cube_side_length::T_L
  r::T_r
  elements::T_el
end


function IsotropicHelicoid(cube_side_length::T_L, r::T_r) where {T_L <: Real, T_r <: Real}
  L = cube_side_length / 6
  return IsotropicHelicoid(cube_side_length, r,
    CompoundBody(
      CompoundBody(
        CylinderSideWall(SA[ 1,  2,  3] * L,            SA[ 2,  1, -3] * L, r),
                  Sphere(SA[ 1,  2,  3] * L, r), Sphere(SA[ 2,  1, -3] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[-2,  1,  3] * L, SA[-1,  2, -3] * L, r),
                  Sphere(SA[-2,  1,  3] * L, r), Sphere(SA[-1,  2, -3] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[-1, -2,  3] * L,            SA[-2, -1, -3] * L, r),
                  Sphere(SA[-1, -2,  3] * L, r), Sphere(SA[-2, -1, -3] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 2, -1,  3] * L,            SA[ 1, -2, -3] * L, r),
                  Sphere(SA[ 2, -1,  3] * L, r), Sphere(SA[ 1, -2, -3] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 2,  3,  1] * L,            SA[ 1, -3,  2] * L, r),
                  Sphere(SA[ 2,  3,  1] * L, r), Sphere(SA[ 1, -3,  2] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[-1,  3,  2] * L,            SA[-2, -3,  1] * L, r),
                  Sphere(SA[-1,  3,  2] * L, r), Sphere(SA[-2, -3,  1] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 1,  3, -2] * L,            SA[ 2, -3, -1] * L, r),
                  Sphere(SA[ 1,  3, -2] * L, r), Sphere(SA[ 2, -3, -1] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[-2,  3, -1] * L,            SA[-1, -3, -2] * L, r),
                  Sphere(SA[-2,  3, -1] * L, r), Sphere(SA[-1, -3, -2] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 3,  1,  2] * L,            SA[-3,  2,  1] * L, r),
                  Sphere(SA[ 3,  1,  2] * L, r), Sphere(SA[-3,  2,  1] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 3, -2,  1] * L,            SA[-3, -1,  2] * L, r),
                  Sphere(SA[ 3, -2,  1] * L, r), Sphere(SA[-3, -1,  2] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 3, -1, -2] * L,            SA[-3, -2, -1] * L, r),
                  Sphere(SA[ 3, -1, -2] * L, r), Sphere(SA[-3, -2, -1] * L, r)),
      CompoundBody(
        CylinderSideWall(SA[ 3,  2, -1] * L,            SA[-3,  1, -2] * L, r),
                  Sphere(SA[ 3,  2, -1] * L, r), Sphere(SA[-3,  1, -2] * L, r)),
    )
  )
end


function collide_with_body(particle::Particle, helicoid::IsotropicHelicoid)
  return collide_with_body(particle, helicoid.elements)
end


"""
    sample_particle_state_uniform_flow(mean_v, std_v, domain_radius) -> Particle

Sample an initial particle state for a gas in thermal equilibrium flowing past a body:

* The velocity `Particle.v` is drawn from a 3D normal distribution with isotropic standard deviation `std_v` and mean drift `mean_v`. The respective PDF with ``\\boldsymbol v_0``, ``\\sigma_v`` and ``\\langle \\boldsymbol v_0 \\rangle`` for `v`, `std_v` and `mean_v`,  is
```math
f(\\boldsymbol v_0) = \\frac{1}{(\\sqrt{2 \\pi \\sigma_v^2})^3} \\exp \\biggl[- \\frac{(\\boldsymbol v_0 - \\langle \\boldsymbol v_0 \\rangle) \\cdot (\\boldsymbol v_0 - \\langle \\boldsymbol v_0 \\rangle)}{2 \\sigma_v^2}\\biggr].
```
* The initial position `Particle.x` is drawn uniformly on the disk of radius `domain_radius`; the disk is oriented perpendicular to `Particle.v` and is located at a distance `domain_radius` from the domain center.
* The time `Particle.t` is set to zero.
"""
function sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
  T = eltype(mean_v)
  v = std_v * SA[randn(T), randn(T), randn(T)] .+ mean_v
  e_v = normalize(v)
  a = SA[randn(T), randn(T), randn(T)]
  b = normalize(a .- (a ⋅ e_v) * e_v)
  x = -e_v * domain_radius .+ b * √rand(T) * domain_radius
  return Particle(v, x, zero(T))
end


"""
    find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)

Repeatedly sample particles from a gas in thermal equilibrium flowing past `body` (via [`sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)`](@ref)) until one collides with `body`.

Arguments:
- `body`: body past which the gas flows.
- `mean_v`: mean (drift) velocity of gas particles.
- `std_v`: isotropic velocity standard deviation of gas particles.
- `domain_radius`: the radius of the spherical simulation domain centered at ``(0, 0, 0)`` that fully contains `body`.

Returns `(particle_before, particle_after, missed_particles_number)`:
- `particle_before`: particle state before the collision.
- `particle_after`: particle state immediately after the first collision with `body`.
- `missed_particles_number`: number of sampled non-colliding particles before the first hit.
"""
function find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)
  missed_particles_number = 0
  
  while true
    particle_before = sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
    particle_after  = collide_with_body(particle_before, body)
    
    if !isnothing(particle_after)
      return particle_before, particle_after, missed_particles_number
    end
    
    missed_particles_number += 1
  end 
end


tensor4(A::SMatrix{3, 3, T, 9}, B::SMatrix{3, 3, T, 9}) where {T} =
  @SArray [A[i,j] * B[k,l] for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

# tensor4(A::SMatrix{3, 3, T, 9}, B::SMatrix{3, 3, T, 9}) where {T} = begin
#   a = SVector{9, T}(Tuple(A))
#   b = SVector{9, T}(Tuple(B))
#   SArray{Tuple{3, 3, 3, 3}, T, 4, 81}(Tuple(a * transpose(b)))
# end

"""
    force_and_torque_uniform_flow(;
      body,
      torque_point,
      mean_v,
      std_v,
      max_collisions,
      domain_radius,
      sample_size)

Monte Carlo estimate of the kinematic force and torque exerted by a uniform flow of Knudsen gas in thermal equilibrium
on `body`. Kinematic force and the torque are those normalized by the gas mass density.

Arguments:
- `body`: body past which the gas flows.
- `torque_point`: point about which the torque is computed.
- `mean_v`: mean (drift) velocity ``\\langle \\boldsymbol v_0 \\rangle`` of gas particles used by [`sample_particle_state_uniform_flow`](@ref).
- `std_v`: isotropic velocity standard deviation used by [`sample_particle_state_uniform_flow`](@ref).
- `max_collisions`: maximum number of successive collisions tracked for each particle (including the first).
- `domain_radius`: radius of the spherical simulation domain centered at ``(0, 0, 0)`` that fully contains `body`.
- `sample_size`: number of *colliding* particles (hits) used in the Monte Carlo estimator.

Returns a named tuple `(F, M, dF_v, dM_v, P_hit, ΣF, ΣM, ΣdF_v, ΣdM_v, ΣP_hit)` with the following Monte Carlo estimates:
- `F`: kinematic force vector.
- `M`: kinematic torque vector about `torque_point`.
- `dF_v`: Jacobian ``\\partial \\boldsymbol F / \\partial \\langle \\boldsymbol v_0 \\rangle`` of the kinematic force with respect to the mean velocity.
- `dM_v`: Jacobian ``\\partial \\boldsymbol M / \\partial \\langle \\boldsymbol v_0 \\rangle`` of the kinematic torque with respect to the mean velocity.
- `P_hit`: probability that a randomly drawn particle hits `body`.
- `ΣF`: covariance tensor of the Monte Carlo estimator `F` (a `3×3` tensor).
- `ΣM`: covariance tensor of the Monte Carlo estimator `M` (a `3×3` tensor).
- `ΣdF_v`: covariance tensor of the Monte Carlo estimator `dF_v` (a `3×3×3×3` tensor).
- `ΣdM_v`: covariance tensor of the Monte Carlo estimator `dM_v` (a `3×3×3×3` tensor).
- `ΣP_hit`: variance of the Monte Carlo estimator `P_hit`.

Notes:
- Particles are generated by repeatedly sampling from [`sample_particle_state_uniform_flow`](@ref) until a collision occurs;
  thus `sample_size` counts hits, while `P_hit` is inferred from the total number of misses encountered along the way.
- The code is multithreaded with `Base.Threads` to take advantage of multi-core CPUs.

# Example

```julia
force_and_torque_uniform_flow(
    body           = Sphere(SA[0.0, 0.0, 1.0], 0.5),
    torque_point   = SA[0.0, 0.0, 0.0],
    mean_v         = SA[1.0, 0.0, 0.0],
    std_v          = 1.0,
    domain_radius  = 2.0,
    max_collisions = 1,
    sample_size    = 10^7)
```

"""
function force_and_torque_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size)
  T = eltype(mean_v)
  
  F_sums     = [zero(SVector{3, T}) for _ in 1:nthreads()]
  M_sums     = [zero(SVector{3, T}) for _ in 1:nthreads()]
  
  dF_v_sums  = [zero(SMatrix{3, 3, T, 9}) for _ in 1:nthreads()]
  dM_v_sums  = [zero(SMatrix{3, 3, T, 9}) for _ in 1:nthreads()]

  F²_sums = [zero(SMatrix{3, 3, T, 9}) for _ in 1:nthreads()]
  M²_sums = [zero(SMatrix{3, 3, T, 9}) for _ in 1:nthreads()]

  dF_v²_sums = [zero(SArray{Tuple{3, 3, 3, 3}, T, 4, 81}) for _ in 1:nthreads()]
  dM_v²_sums = [zero(SArray{Tuple{3, 3, 3, 3}, T, 4, 81}) for _ in 1:nthreads()]
  
  missed_sums  = zeros(Int, nthreads())

  inv_var_v = 1 / (std_v^2) 

  @threads for i in 1:sample_size
    p_b, p_a, missed = find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)
    tid = threadid()

    Δv = p_a.v - p_b.v
    w  = norm(p_b.v)
    s  = (p_b.v - mean_v) * inv_var_v

    Δf = w * Δv
    Δm = w * (p_a.x - torque_point) × Δv

    missed_sums[tid]  += missed
    
    for _ in 2:max_collisions
      p_b = p_a
      p_a = collide_with_body(p_b, body)
      
      if isnothing(p_a) break end
      
      Δv = p_a.v - p_b.v

      Δf += w * Δv
      Δm += w * (p_a.x - torque_point) × Δv
    end

    Δdf_v = Δf * s'
    Δdm_v = Δm * s'
    
    F_sums[tid] += Δf
    M_sums[tid] += Δm

    dF_v_sums[tid] += Δdf_v
    dM_v_sums[tid] += Δdm_v
  
    F²_sums[tid] += Δf * Δf'
    M²_sums[tid] += Δm * Δm'

    dF_v²_sums[tid] += tensor4(Δdf_v, Δdf_v)
    dM_v²_sums[tid] += tensor4(Δdm_v, Δdm_v)
  end

  sample_size_T = T(sample_size)
  
  P_hit  = sample_size_T / (sample_size_T + sum(missed_sums))
  ΣP_hit = P_hit^2 * (1 - P_hit) / sample_size_T
  
  A  =  P_hit * π * domain_radius^2
  A² = (P_hit^2 + ΣP_hit) * (π * domain_radius^2)^2
  
  F  = - sum(F_sums) / sample_size_T * A
  M  = - sum(M_sums) / sample_size_T * A

  dF_v = - sum(dF_v_sums) / sample_size_T * A
  dM_v = - sum(dM_v_sums) / sample_size_T * A
  
  F² =   sum(F²_sums) / sample_size_T * A²
  M² =   sum(M²_sums) / sample_size_T * A²

  dF_v² = sum(dF_v²_sums) / sample_size_T * A²
  dM_v² = sum(dM_v²_sums) / sample_size_T * A²

  ΣF = (F² - F * F') / (sample_size_T - 1) + ΣP_hit / P_hit^2 * F * F'
  ΣM = (M² - M * M') / (sample_size_T - 1) + ΣP_hit / P_hit^2 * M * M'

  ΣdF_v = (dF_v² - tensor4(dF_v, dF_v)) / (sample_size_T - 1) + ΣP_hit / P_hit^2 * tensor4(dF_v, dF_v)
  ΣdM_v = (dM_v² - tensor4(dM_v, dM_v)) / (sample_size_T - 1) + ΣP_hit / P_hit^2 * tensor4(dM_v, dM_v)

  return ( F =  F,  M =  M,  dF_v =  dF_v,  dM_v =  dM_v,  P_hit = P_hit,
          ΣF = ΣF, ΣM = ΣM, ΣdF_v = ΣdF_v, ΣdM_v = ΣdM_v, ΣP_hit = ΣP_hit)
end


end # module GBGKnudsenGas