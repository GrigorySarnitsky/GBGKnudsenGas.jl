"""
    GBGKnudsenGas

Simulate Knudsen flow of particles past a stationary body, where the particle mean free path is much larger than the body size. In this regime, particles interact with the body but not with each other. The body is assumed to be massive, so particle collisions do not move it.

The module provides building blocks to create your own flow setups:

- geometric body primitives (currently [`Sphere`](@ref) and [`CylinderSideWall`](@ref)),
- composition of bodies via [`CompoundBody`](@ref),
- single-particle collision queries via [`collide_with_body`](@ref),
- tools for modelling a gas in thermal equilibrium flowing past a body.

## Collision model

All particle–body collisions use a [Maxwell reflection law](https://doi.org/10.1146/annurev-fluid-121021-014808) with diffuse fraction ``\\alpha`` and wall thermal velocity scale ``\\sigma_{v_\\text{wall}}``:

- ``\\alpha = 0``: purely specular reflection,
- ``\\alpha = 1``: purely diffuse reflection,
- ``0 < \\alpha < 1``: Bernoulli mixture of specular and diffuse reflection.

A gas particle that reflects diffusely has an outgoing velocity ``\\boldsymbol v`` distributed according to the following probability density function, where ``\\boldsymbol n`` denotes the outward wall normal:

```math
f_\\text{diff}(\\boldsymbol v) = \\left\\{
\\begin{aligned}
&\\frac{\\boldsymbol v \\cdot \\boldsymbol n}{2 \\pi \\sigma_{v_\\text{wall}}^4} \\exp \\biggl(- \\frac{\\boldsymbol v \\cdot \\boldsymbol v}{2 \\sigma_{v_\\text{wall}}^2} \\biggr) && \\text{if \$ \\boldsymbol v \\cdot \\boldsymbol n > 0\$},\\\\
& 0 && \\text{otherwise.}
\\end{aligned}\\right.
```

## Conventions

- Cartesian coordinates are used throughout the package.
- The package is purely kinematic: it works with distances, times, and velocities, but not masses or momenta.
- Particles can reflect only from the exterior side of a body’s surface, so if a particle ever enters a body due to an error, it will freely exit again.
"""
module GBGKnudsenGas

export Particle, Body, Sphere, CylinderSideWall, IsotropicHelicoidRods, CompoundBody, collide_with_body,
 sample_particle_state_uniform_flow, find_colliding_particle_uniform_flow, force_and_torque_uniform_flow


using Random
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


function reflect_velocity(v, n, α, wall_std_v)
  if α ≤ 0
    return reflect_velocity_specularly(v, n)
  elseif α ≥ 1
    return reflect_velocity_diffusively(n, wall_std_v)
  else
    return  rand() < α ? reflect_velocity_diffusively(n, wall_std_v) : reflect_velocity_specularly(v, n)
  end
end


function reflect_velocity_specularly(v, n)
  return v - 2 * n * (v ⋅ n)
end


function reflect_velocity_diffusively(n, wall_std_v)
  T = typeof(wall_std_v)
  gaussian_v = wall_std_v * SA[randn(T), randn(T), randn(T)]
  rayleigh_v = wall_std_v * √(2 * randexp(T))
  return gaussian_v + (rayleigh_v - gaussian_v ⋅ n) * n
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
    Sphere(x::SVector{3}, r; α = zero(r), wall_std_v = one(r))

Sphere of radius `r` centered at point `x`.

Gas particles reflect from the sphere via the Maxwell reflection law with diffuse fraction `α` and wall thermal velocity scale `wall_std_v`.
"""
struct Sphere{T <: Real} <: Body
  x::SVector{3, T}
  r::T
  r²::T
  α::T
  wall_std_v::T
end


function Sphere(x::SVector{3}, r; α = zero(r), wall_std_v = one(r))
  T = promote_type(eltype(x), typeof(r), typeof(α), typeof(wall_std_v))
  Sphere(SVector{3,T}(x), T(r), T(r)^2, T(α), T(wall_std_v))
end


function collide_with_body(particle::Particle, sphere::Sphere)
  Δx = particle.x - sphere.x
  Δx², v² = Δx ⋅ Δx, particle.v ⋅ particle.v
  b = particle.v ⋅ Δx
  d = b^2 - v² * (Δx² - sphere.r²)
  
  if b > 0 || d < 0 || Δx² < sphere.r² || iszero(v²) return nothing end
  
  time_to_collision = (-b - √d) / v²
  x_at_collision  = particle.x + particle.v * time_to_collision
  Δx_at_collision = x_at_collision - sphere.x
  
  n = Δx_at_collision / sphere.r
  v_after_collision = reflect_velocity(particle.v, n, sphere.α, sphere.wall_std_v)
  
  return Particle(v_after_collision, x_at_collision, particle.t + time_to_collision)
end


"""
    CylinderSideWall <: Body
    CylinderSideWall(x_1::SVector{3}, x_2::SVector{3}, r; α = zero(r), wall_std_v = one(r))

Finite circular cylinder defined by its axis endpoints `x_1`, `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. End caps can be added via [`CompoundBody`](@ref), like we do in [`IsotropicHelicoidRods`](@ref).

Gas particles reflect from the side wall via the Maxwell reflection law with diffuse fraction `α` and wall thermal velocity scale `wall_std_v`.
"""
struct CylinderSideWall{T <: Real} <: Body
  x_1::SVector{3, T}
  x_2::SVector{3, T}
  r::T
  axis::SVector{3, T}
  l::T
  r²::T
  α::T
  wall_std_v::T
end


function CylinderSideWall(x_1::SVector{3}, x_2::SVector{3}, r; α = zero(r), wall_std_v = one(r))
  T = promote_type(eltype(x_1), eltype(x_2), typeof(r), typeof(α), typeof(wall_std_v))
  axis = normalize(SVector{3,T}(x_2) - SVector{3,T}(x_1))
  l = T(norm(SVector{3,T}(x_2) - SVector{3,T}(x_1)))
  CylinderSideWall(SVector{3,T}(x_1), SVector{3,T}(x_2), T(r), axis, l, T(r)^2, T(α), T(wall_std_v))
end


function collide_with_body(particle::Particle, cylinder::CylinderSideWall)
  Δx = particle.x - cylinder.x_1
  Δx_⊥ = Δx - (Δx ⋅ cylinder.axis) * cylinder.axis
  v_⊥  = particle.v - (particle.v ⋅ cylinder.axis) * cylinder.axis
  Δx_⊥², v_⊥² = Δx_⊥ ⋅ Δx_⊥, v_⊥ ⋅ v_⊥
  b = v_⊥ ⋅ Δx_⊥
  d = b^2 - v_⊥² * (Δx_⊥² - cylinder.r²)
  
  if b > 0 || d < 0 || Δx_⊥² < cylinder.r² || iszero(v_⊥²) return nothing end
  
  time_to_collision = (-b - √d) / v_⊥²
  x_at_collision  = particle.x + particle.v * time_to_collision
  Δx_at_collision = x_at_collision - cylinder.x_1
  
  if !(0 < Δx_at_collision ⋅ cylinder.axis < cylinder.l) return nothing end

  Δx_⊥_at_collision = Δx_at_collision - (Δx_at_collision ⋅ cylinder.axis) * cylinder.axis
  n = Δx_⊥_at_collision / cylinder.r

  v_after_collision = reflect_velocity(particle.v, n, cylinder.α, cylinder.wall_std_v)
  
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


"""
    IsotropicHelicoidRods <: Body
    IsotropicHelicoidRods(cube_side_length, r; α = zero(r), wall_std_v = one(r))

Model of Lord Kelvin's isotropic helicoid, constructed from 12 rods. Each rod is a cylinder of radius `r`, capped with two spheres of the same radius. The helicoid is inscribed in a cube of side length `cube_side_length + 2r`, centered at ``(0, 0, 0)``. This cube side length indirectly defines the length of each rod.

Gas particles reflect from the helicoid via the Maxwell reflection law with diffuse fraction `α` and wall thermal velocity scale `wall_std_v`.

The rods that form the helicoid can be accessed via `helicoid.elements`, which is a `CompoundBody`.
"""
struct IsotropicHelicoidRods{T <: Real, T_el <: CompoundBody} <: Body
  cube_side_length::T
  r::T
  elements::T_el
end


function IsotropicHelicoidRods(cube_side_length, r; α = zero(r), wall_std_v = one(r))
  L = cube_side_length / 6
  return IsotropicHelicoidRods(cube_side_length, r,
    CompoundBody(
      CompoundBody(
        CylinderSideWall(SA[ 1,  2,  3] * L,                           SA[ 2,  1, -3] * L, r; α, wall_std_v),
                  Sphere(SA[ 1,  2,  3] * L, r; α, wall_std_v), Sphere(SA[ 2,  1, -3] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[-2,  1,  3] * L,                           SA[-1,  2, -3] * L, r; α, wall_std_v),
                  Sphere(SA[-2,  1,  3] * L, r; α, wall_std_v), Sphere(SA[-1,  2, -3] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[-1, -2,  3] * L,                           SA[-2, -1, -3] * L, r; α, wall_std_v),
                  Sphere(SA[-1, -2,  3] * L, r; α, wall_std_v), Sphere(SA[-2, -1, -3] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 2, -1,  3] * L,                           SA[ 1, -2, -3] * L, r; α, wall_std_v),
                  Sphere(SA[ 2, -1,  3] * L, r; α, wall_std_v), Sphere(SA[ 1, -2, -3] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 2,  3,  1] * L,                           SA[ 1, -3,  2] * L, r; α, wall_std_v),
                  Sphere(SA[ 2,  3,  1] * L, r; α, wall_std_v), Sphere(SA[ 1, -3,  2] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[-1,  3,  2] * L,                           SA[-2, -3,  1] * L, r; α, wall_std_v),
                  Sphere(SA[-1,  3,  2] * L, r; α, wall_std_v), Sphere(SA[-2, -3,  1] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 1,  3, -2] * L,                           SA[ 2, -3, -1] * L, r; α, wall_std_v),
                  Sphere(SA[ 1,  3, -2] * L, r; α, wall_std_v), Sphere(SA[ 2, -3, -1] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[-2,  3, -1] * L,                           SA[-1, -3, -2] * L, r; α, wall_std_v),
                  Sphere(SA[-2,  3, -1] * L, r; α, wall_std_v), Sphere(SA[-1, -3, -2] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 3,  1,  2] * L,                           SA[-3,  2,  1] * L, r; α, wall_std_v),
                  Sphere(SA[ 3,  1,  2] * L, r; α, wall_std_v), Sphere(SA[-3,  2,  1] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 3, -2,  1] * L,                           SA[-3, -1,  2] * L, r; α, wall_std_v),
                  Sphere(SA[ 3, -2,  1] * L, r; α, wall_std_v), Sphere(SA[-3, -1,  2] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 3, -1, -2] * L,                           SA[-3, -2, -1] * L, r; α, wall_std_v),
                  Sphere(SA[ 3, -1, -2] * L, r; α, wall_std_v), Sphere(SA[-3, -2, -1] * L, r; α, wall_std_v)),
      CompoundBody(
        CylinderSideWall(SA[ 3,  2, -1] * L,                           SA[-3,  1, -2] * L, r; α, wall_std_v),
                  Sphere(SA[ 3,  2, -1] * L, r; α, wall_std_v), Sphere(SA[-3,  1, -2] * L, r; α, wall_std_v)),
    )
  )
end


function collide_with_body(particle::Particle, helicoid::IsotropicHelicoidRods)
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
  return Particle(v, x, zero(eltype(v)))
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



function force_and_torque_sums_uniform_flow_serial(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size, T)
  F_sum = zero(SVector{3, T})
  M_sum = zero(SVector{3, T})
  
  dF_v_sum = zero(SMatrix{3, 3, T, 9})
  dM_v_sum = zero(SMatrix{3, 3, T, 9})

  F²_sum = zero(SMatrix{3, 3, T, 9})
  M²_sum = zero(SMatrix{3, 3, T, 9})

  dF_v²_sum = zero(SArray{Tuple{3, 3, 3, 3}, T, 4, 81})
  dM_v²_sum = zero(SArray{Tuple{3, 3, 3, 3}, T, 4, 81})
  
  missed_sum = zero(Int)

  for i in 1:sample_size
    p_b, p_a, missed = find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)

    Δv = p_a.v - p_b.v
    w  = norm(p_b.v)
    s  = p_b.v - mean_v

    Δf = w * Δv
    Δm = w * (p_a.x - torque_point) × Δv

    missed_sum += missed
    
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
    
    F_sum += Δf
    M_sum += Δm

    dF_v_sum += Δdf_v
    dM_v_sum += Δdm_v
  
    F²_sum += Δf * Δf'
    M²_sum += Δm * Δm'

    dF_v²_sum += tensor4(Δdf_v, Δdf_v)
    dM_v²_sum += tensor4(Δdm_v, Δdm_v)
  end

  N_hit  = T(sample_size)
  N_miss = T(missed_sum)
  
  return (;F_sum, M_sum, dF_v_sum, dM_v_sum, F²_sum, M²_sum, dF_v²_sum, dM_v²_sum, N_hit, N_miss)
end


function force_and_torque_sums_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size, T, parallel)
  if parallel
    number_of_threads = Threads.nthreads(:default)
    number_of_chunks = number_of_threads
    base_chunk_size, remainder = divrem(sample_size, number_of_chunks)
    
    parallel_force_and_torque_sums = Array{Any, 1}(undef, number_of_chunks)

    @threads for i in collect(eachindex(parallel_force_and_torque_sums))
      chunk_sample_size = base_chunk_size + (i <= remainder)
      parallel_force_and_torque_sums[i] = force_and_torque_sums_uniform_flow_serial(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size = chunk_sample_size, T)
    end
    
    field_names  = keys(parallel_force_and_torque_sums[1])
    return (; (name => sum(d[name] for d in parallel_force_and_torque_sums) for name in field_names)...)
  else
    return force_and_torque_sums_uniform_flow_serial(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size, T)
  end
end

  

"""
    force_and_torque_uniform_flow(;
      body,
      torque_point,
      mean_v,
      std_v,
      max_collisions,
      domain_radius,
      sample_size,
      parallel = false
    )

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
- `parallel`: if `true`, use [Julia multi-threading](https://docs.julialang.org/en/v1/base/multi-threading/) to evaluate samples in parallel.

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
- The covariances (all the `Σ`'s) are computed as described in the [development tips](development.md).

# Example

```julia
force_and_torque_uniform_flow(
  body           = Sphere(SA[0.0, 0.0, 1.0], 0.5),
  torque_point   = SA[0.0, 0.0, 0.0],
  mean_v         = SA[1.0, 0.0, 0.0],
  std_v          = 1.0,
  domain_radius  = 2.0,
  max_collisions = 1,
  sample_size    = 10^7,
  parallel       = true,
)
```
"""
function force_and_torque_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size, parallel = false)
  T = promote_type(eltype(mean_v), eltype(torque_point), typeof(std_v), typeof(domain_radius))

  d = force_and_torque_sums_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size, T, parallel)
  
  P_hit  = (d.N_hit - 1) / (d.N_miss + d.N_hit - 1)
  ΣP_hit = P_hit^2 * (1 - P_hit) / d.N_hit
  
  A = T(P_hit * π * domain_radius^2)
  B = T(P_hit * π * domain_radius^2 / std_v^2)
  
  F = - d.F_sum / d.N_hit * A
  M = - d.M_sum / d.N_hit * A

  dF_v = - d.dF_v_sum / d.N_hit * B
  dM_v = - d.dM_v_sum / d.N_hit * B
  
  F² = d.F²_sum / d.N_hit * A^2
  M² = d.M²_sum / d.N_hit * A^2

  dF_v² = d.dF_v²_sum / d.N_hit * B^2
  dM_v² = d.dM_v²_sum / d.N_hit * B^2

  ΣF = (F² - F * F') / (d.N_hit - 1) + ΣP_hit / P_hit^2 * F * F'
  ΣM = (M² - M * M') / (d.N_hit - 1) + ΣP_hit / P_hit^2 * M * M'

  ΣdF_v = (dF_v² - tensor4(dF_v, dF_v)) / (d.N_hit - 1) + ΣP_hit / P_hit^2 * tensor4(dF_v, dF_v)
  ΣdM_v = (dM_v² - tensor4(dM_v, dM_v)) / (d.N_hit - 1) + ΣP_hit / P_hit^2 * tensor4(dM_v, dM_v)

  return (; F,  M,  dF_v,  dM_v,  P_hit,
           ΣF, ΣM, ΣdF_v, ΣdM_v, ΣP_hit) 
end


end # module GBGKnudsenGas