"""
    GBGKnudsenGas

Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a special body (isotropic helicoid) is included.

Currently, all particle–body collisions are assumed to be specular.
"""
module GBGKnudsenGas


export Particle, Body, Sphere, CylinderSideWall, IsotropicHelicoid, CompoundBody, CompoundBodyDynamic,
 collide_with_body, get_initial_particle_state


using StaticArrays
using LinearAlgebra


"""
    Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.
"""
struct Particle
  v::SVector{3, Float64}
  x::SVector{3, Float64}
  t::Float64
end


"""
    Body

Abstract supertype for all bodies that can interact with a `Particle` via `collide_with_body(::Particle, ::Body)`.
"""
abstract type Body end


"""
    collide_with_body(particle::Particle, body::Body) -> Union{Nothing,Particle}

Compute the next collision of `particle` with `body`, assuming straight-line motion between collisions.

Return a new `Particle` located at the collision point, with updated velocity (after collision) and time, or `nothing` if there is no collision.
"""
collide_with_body(::Particle, ::Body)::Union{Nothing,Particle} = error("No collide_with_body defined for this body")


"""
    Sphere <: Body
    Sphere(x::SVector{3, Float64}, x::Float64)

Sphere of radius `r` centered at point `x`.
"""
struct Sphere <: Body
  x::SVector{3, Float64}
  r::Float64
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

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. Endcaps can be added via the `CompoundBody`, like we do in `IsotropicHelicopiod`.
"""
struct CylinderSideWall <: Body
  x_1::SVector{3, Float64}
  x_2::SVector{3, Float64}
  r  ::Float64
end

function collide_with_body(particle::Particle, cylinder::CylinderSideWall)::Union{Nothing, Particle}
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

struct CompoundBodyDynamic <: Body
  bodies::Array{Body, 1}
end

function collide_with_body(particle::Particle, compound_body::CompoundBodyDynamic)::Union{Nothing, Particle}
  nearest_collision_particle::Union{Nothing, Particle} = nothing
  nearest_collision_t::Float64 = Inf
  for body in compound_body.bodies
    collided_particle = collide_with_body(particle, body)
    if collided_particle !== nothing && collided_particle.t < nearest_collision_t
      nearest_collision_t = collided_particle.t
      nearest_collision_particle = collided_particle
    end
  end
  return nearest_collision_particle
end

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

struct IsotropicHelicoid <: Body
  cube_side_length::Float64
  r::Float64
  elements::CompoundBody
end

IsotropicHelicoid(cube_side_length, r) = IsotropicHelicoid(cube_side_length, r,
  CompoundBody(
    CompoundBody(
      CylinderSideWall(SA[ 1,  2,  3] / 6 * cube_side_length, SA[ 2,  1, -3 ]/ 6 * cube_side_length, r),
      Sphere(SA[ 1,  2,  3] / 6 * cube_side_length, r), Sphere(SA[ 2,  1, -3] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[-2,  1,  3] / 6 * cube_side_length, SA[-1,  2, -3] / 6 * cube_side_length, r),
      Sphere(SA[-2,  1,  3] / 6 * cube_side_length, r), Sphere(SA[-1,  2, -3] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[-1, -2,  3] / 6 * cube_side_length, SA[-2, -1, -3] / 6 * cube_side_length, r),
      Sphere(SA[-1, -2,  3] / 6 * cube_side_length, r), Sphere(SA[-2, -1, -3] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 2, -1,  3] / 6 * cube_side_length, SA[ 1, -2, -3] / 6 * cube_side_length, r),
      Sphere(SA[ 2, -1,  3] / 6 * cube_side_length, r), Sphere(SA[ 1, -2, -3] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 2,  3,  1] / 6 * cube_side_length, SA[ 1, -3,  2] / 6 * cube_side_length, r),
      Sphere(SA[ 2,  3,  1] / 6 * cube_side_length, r), Sphere(SA[ 1, -3,  2] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[-1,  3,  2] / 6 * cube_side_length, SA[-2, -3,  1] / 6 * cube_side_length, r),
      Sphere(SA[-1,  3,  2] / 6 * cube_side_length, r), Sphere(SA[-2, -3,  1] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 1,  3, -2] / 6 * cube_side_length, SA[ 2, -3, -1] / 6 * cube_side_length, r),
      Sphere(SA[ 1,  3, -2] / 6 * cube_side_length, r), Sphere(SA[ 2, -3, -1] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[-2,  3, -1] / 6 * cube_side_length, SA[-1, -3, -2] / 6 * cube_side_length, r),
      Sphere(SA[-2,  3, -1] / 6 * cube_side_length, r), Sphere(SA[-1, -3, -2] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 3,  1,  2] / 6 * cube_side_length, SA[-3,  2,  1] / 6 * cube_side_length, r),
      Sphere(SA[ 3,  1,  2] / 6 * cube_side_length, r), Sphere(SA[-3,  2,  1] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 3, -2,  1] / 6 * cube_side_length, SA[-3, -1,  2] / 6 * cube_side_length, r),
      Sphere(SA[ 3, -2,  1] / 6 * cube_side_length, r), Sphere(SA[-3, -1,  2] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 3, -1, -2] / 6 * cube_side_length, SA[-3, -2, -1] / 6 * cube_side_length, r),
      Sphere(SA[ 3, -1, -2] / 6 * cube_side_length, r), Sphere(SA[-3, -2, -1] / 6 * cube_side_length, r)),
    CompoundBody(
      CylinderSideWall(SA[ 3,  2, -1] / 6 * cube_side_length, SA[-3,  1, -2] / 6 * cube_side_length, r),
      Sphere(SA[ 3,  2, -1] / 6 * cube_side_length, r), Sphere(SA[-3,  1, -2] / 6 * cube_side_length, r)),
  ))

function collide_with_body(particle::Particle, helicoid::IsotropicHelicoid)::Union{Nothing, Particle}
  return collide_with_body(particle, helicoid.elements)
end

function get_initial_particle_state(mean_v, std_v, domain_radius)
  ϕ, θ = 2*π*rand(), acos(2*rand() - 1)
  x = domain_radius * SA[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
  v = std_v * SA[randn(), randn(), randn()] + mean_v
  Particle(v, x, 0.0)
end

end # module GBGKnudsenGas