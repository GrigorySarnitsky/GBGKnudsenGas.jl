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
    Particle(v::SVector{3, T}, x::SVector{3, T}, t::T) where {T<:Real}

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.
"""
struct Particle{T<:Real}
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
collide_with_body(::Particle, ::Body)::Union{Nothing, Particle} = error("No collide_with_body defined for this body")


"""
    Sphere <: Body
    Sphere(x::SVector{3, Float64}, r::Float64)

Sphere of radius `r` centered at point `x`.
"""
struct Sphere{T_x<:Real, T_r<:Real} <: Body
  x::SVector{3,T_x}
  r::T_r
end


function collide_with_body(particle::Particle, sphere::Sphere)::Union{Nothing, Particle}
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
    CylinderSideWall(x_1::SVector{3, Float64}, x_2::SVector{3, Float64}, r::Float64)

Finite circular cylinder defined by its axis endpoints `x_1` and `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. Endcaps can be added via `CompoundBody`, like we do in `IsotropicHelicoid`.
"""
struct CylinderSideWall{T_x<:Real, T_r<:Real} <: Body
  x_1::SVector{3,T_x}
  x_2::SVector{3,T_x}
  r  ::T_r
end


function collide_with_body(particle::Particle, cylinder::CylinderSideWall)::Union{Nothing, Particle}
  Δx = particle.x - cylinder.x_1
  x_2_1 = cylinder.x_2 - cylinder.x_1
  x_2_1² = x_2_1 ⋅ x_2_1

  ###
  iszero(x_2_1²) && return nothing
  
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
    CompoundBody{T<:Tuple{Vararg{<:Body}}}

A compound body whose components are stored in the tuple `bodies`. A collision with a `CompoundBody` is defined as the earliest collision with any of its components.

This type should be the default choice for representing compound bodies, as it is usually faster.

Use the convenience constructors

    CompoundBody(bodies::Vararg{<:Body})
    CompoundBody(bodies::AbstractVector{<:Body})

to build instances.
"""
struct CompoundBody{T<:Tuple{Vararg{<:Body}}} <: Body
  bodies::T
end


CompoundBody(bodies::Vararg{<:Body}) = CompoundBody(bodies)
CompoundBody(bodies::AbstractVector{<:Body}) = CompoundBody(Tuple(bodies))


function earlier(a::Union{Nothing, Particle}, b::Union{Nothing, Particle})
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

# function collide_with_body(particle::Particle, compound_body::CompoundBodyDynamic)::Union{Nothing, Particle}
#   nearest_collision_particle::Union{Nothing, Particle} = nothing
#   for body in compound_body.bodies
#     collided_particle = collide_with_body(particle, body)
#     nearest_collision_particle = earlier(nearest_collision_particle, collided_particle)
#   end
#   return nearest_collision_particle
# end


"""
    IsotropicHelicoid <: Body
    IsotropicHelicoid(cube_side_length::Float64, r::Float64)

A model of Lord Kelvin's isotropic helicoid, constructed from 12 rods. Each rod is a cylinder of radius `r`, capped with two spheres of the same radius. The helicoid is inscribed in a cube of side length `cube_side_length`. This cube side length indirectly defines the length of each rod.

The rods that form the helicoid can be accessed via `helicoid.elements`, which is a `CompoundBody`.
"""
struct IsotropicHelicoid{T_L<:Real, T_r<:Real, T_el<:CompoundBody} <: Body
  cube_side_length::T_L
  r::T_r
  elements::T_el
end


function IsotropicHelicoid(cube_side_length::T_L, r::T_r) where {T_L<:Real, T_r<:Real}
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


function collide_with_body(particle::Particle, helicoid::IsotropicHelicoid)::Union{Nothing, Particle}
  return collide_with_body(particle, helicoid.elements)
end


"""
    sample_particle_state_uniform_flow(mean_v, std_v, domain_radius) -> Particle

Sample an initial particle state for a gas in thermal equilibrium flowing past a body:

* The velocity `Particle.v` is drawn from a 3D normal distribution with isotropic standard deviation `std_v` and mean drift `mean_v`.
* The initial position `Particle.x` is drawn uniformly on the disk of radius `domain_radius`; the disk is oriented perpendicular to `Particle.v` and is located at a distance `domain_radius` from the domain center.
* The time `Particle.t` is set to `0.0`.
"""
function sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
  T = promote_type(eltype(mean_v), typeof(std_v), typeof(domain_radius))
  v = std_v * SA[randn(T), randn(T), randn(T)] .+ mean_v

  ###
  iszero(v ⋅ v) && return sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
  
  e_v = normalize(v)
  a = SA[randn(T), randn(T), randn(T)]

  ###
  a_⊥ = a .- (a ⋅ e_v) * e_v
  iszero(a_⊥ ⋅ a_⊥) && return sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
  b = normalize(a_⊥)

  #b = normalize(a .- (a ⋅ e_v) * e_v)
  x = -e_v * domain_radius .+ b * √rand(T) * domain_radius
  return Particle(v, x, zero(T))
end


"""
    find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)

Repeatedly draw particles from [`sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)`](@ref) until one collides with `body`.

Arguments:
- `body`: body past which the gas flows.
- `mean_v`: mean velocity of gas particles.
- `std_v`: isotropic velocity standard deviation of gas particles.
- `domain_radius`: the radius of a spherical simulation domain centered at ``(0, 0, 0)``, which the body must fully fit into.

Returns tuple `(particle_before, particle_after, missed_particles_number)` where:

- `particle_before::Particle`: particle state before the collision.
- `particle_after::Particle`: particle state after the collision.
- `missed_particles_number::Int64` number of non-colliding particles encountered before a hit was found.
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


"""
    force_and_torque_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size)

Calculate force and torque a flow of a Knudsen gas at thermal equilibrium excerts on a body. Both the force and the torque are normalized by the gas density (``\\boldsymbol F\rho``, ``\\boldsymbol M\rho``).

Arguments:
- `body`: body past which the gas flows.
- `torque_point`: point with respect to which the flow is computed.
- `mean_v`: mean velocity of gas particles.
- `std_v`: isotropic velocity standard deviation of gas particles.
- `max_collisions`: maximum number of collisions of the particle with the body. The first one is always computed, but
- `domain_radius`: the radius of a spherical simulation domain centered at ``(0, 0, 0)``, which the body must fully fit into.
- `sample_size`: number of colliding particles used for the Monte Carlo estimation of the force and torque.

The code is multithreaded with `Base.Threads` to take advantage of multi-core CPUs.
"""
function force_and_torque_uniform_flow(;body, torque_point, mean_v, std_v, max_collisions, domain_radius, sample_size)
  T = eltype(mean_v)
  F_sums      = [SVector{3, T}(0, 0, 0) for _ in 1:nthreads()]
  M_sums      = [SVector{3, T}(0, 0, 0) for _ in 1:nthreads()]
  missed_sums = zeros(Int, nthreads())

  @threads for i in 1:sample_size
    p_b, p_a, missed = find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)
    tid = threadid()

    Δv = p_a.v - p_b.v
    w  = norm(p_b.v)
    F_sums[tid] += w * Δv
    M_sums[tid] += w * (p_a.x - torque_point) × Δv
    missed_sums[tid] += missed
    
    for _ in 2:max_collisions
      p_b = p_a
      p_a = collide_with_body(p_b, body)
      
      if isnothing(p_a) break end
      
      Δv = p_a.v - p_b.v
      F_sums[tid] += w * Δv
      M_sums[tid] += w * (p_a.x - torque_point) × Δv
    end
  end

  collision_probability = T(sample_size) / (T(sample_size) + sum(missed_sums))
  
  F = - sum(F_sums) / T(sample_size) * collision_probability * π * domain_radius^2
  M = - sum(M_sums) / T(sample_size) * collision_probability * π * domain_radius^2

  return (F = F, M = M, P_hit = collision_probability)
end


end # module GBGKnudsenGas