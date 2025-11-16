module GBGKnudsenGas

export Particle, Body, Sphere, CylinderSideWall, IsotropicHelicoid, CompoundBody, collide_with_body, get_initial_particle_state

using StaticArrays
using LinearAlgebra

struct Particle
  t::Float64
  x::SVector{3, Float64}
  v::SVector{3, Float64}
end

abstract type Body end

collide_with_body(::Particle, ::Body) = error("No collide_with_body defined for this body")

struct Sphere <: Body
  x::SVector{3, Float64}
  r::Float64
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
  # heck performance by changing Δx_at_collision ⋅ Δx_at_collision to sphere.r^2
  v_after_collision = particle.v - 2 * Δx_at_collision * (particle.v ⋅ Δx_at_collision) / (Δx_at_collision ⋅ Δx_at_collision)
  
  return Particle(particle.t + time_to_collision, x_at_collision, v_after_collision)
end

struct CylinderSideWall <: Body
  x_1::SVector{3, Float64}
  x_2::SVector{3, Float64}
  r  ::Float64
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
  # check performance by changing Δx_⊥_at_collision ⋅ Δx_⊥_at_collision to cylinder.r^2
  v_after_collision = particle.v - 2 * Δx_⊥_at_collision * (particle.v ⋅ Δx_⊥_at_collision) / (Δx_⊥_at_collision ⋅ Δx_⊥_at_collision)
  
  return Particle(particle.t + time_to_collision, x_at_collision, v_after_collision)
end

struct CompoundBody <: Body
  bodies::Array{Body, 1}
end

function collide_with_body(particle::Particle, compound_body::CompoundBody)
  # change to a loop later
  # best_particle = nothing
  # best_t = Inf
  # for body in compound_body.bodies
  #   p′ = collide_with_body(particle, body)
  #   if p′ !== nothing && p′.t < best_t
  #     best_t = p′.t
  #     best_particle = p′
  #   end
  # end
  # return best_particle
  
  collisions = filter(!isnothing, [collide_with_body(particle, body) for body in compound_body.bodies])
  if isempty(collisions) return nothing end
  return collisions[argmin(p.t for p in collisions)]
end

struct IsotropicHelicoid <: Body
  cube_side_length::Float64
  r::Float64
  elements::CompoundBody
end

IsotropicHelicoid(l, r) = IsotropicHelicoid(l, r,
  CompoundBody([
    CompoundBody([CylinderSideWall(SA[ 1,  2,  3]/6*l, SA[ 2,  1, -3]/6*l, r), Sphere(SA[ 1,  2,  3]/6*l, r), Sphere(SA[ 2,  1, -3]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[-2,  1,  3]/6*l, SA[-1,  2, -3]/6*l, r), Sphere(SA[-2,  1,  3]/6*l, r), Sphere(SA[-1,  2, -3]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[-1, -2,  3]/6*l, SA[-2, -1, -3]/6*l, r), Sphere(SA[-1, -2,  3]/6*l, r), Sphere(SA[-2, -1, -3]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 2, -1,  3]/6*l, SA[ 1, -2, -3]/6*l, r), Sphere(SA[ 2, -1,  3]/6*l, r), Sphere(SA[ 1, -2, -3]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 2,  3,  1]/6*l, SA[ 1, -3,  2]/6*l, r), Sphere(SA[ 2,  3,  1]/6*l, r), Sphere(SA[ 1, -3,  2]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[-1,  3,  2]/6*l, SA[-2, -3,  1]/6*l, r), Sphere(SA[-1,  3,  2]/6*l, r), Sphere(SA[-2, -3,  1]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 1,  3, -2]/6*l, SA[ 2, -3, -1]/6*l, r), Sphere(SA[ 1,  3, -2]/6*l, r), Sphere(SA[ 2, -3, -1]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[-2,  3, -1]/6*l, SA[-1, -3, -2]/6*l, r), Sphere(SA[-2,  3, -1]/6*l, r), Sphere(SA[-1, -3, -2]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 3,  1,  2]/6*l, SA[-3,  2,  1]/6*l, r), Sphere(SA[ 3,  1,  2]/6*l, r), Sphere(SA[-3,  2,  1]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 3, -2,  1]/6*l, SA[-3, -1,  2]/6*l, r), Sphere(SA[ 3, -2,  1]/6*l, r), Sphere(SA[-3, -1,  2]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 3, -1, -2]/6*l, SA[-3, -2, -1]/6*l, r), Sphere(SA[ 3, -1, -2]/6*l, r), Sphere(SA[-3, -2, -1]/6*l, r)]),
    CompoundBody([CylinderSideWall(SA[ 3,  2, -1]/6*l, SA[-3,  1, -2]/6*l, r), Sphere(SA[ 3,  2, -1]/6*l, r), Sphere(SA[-3,  1, -2]/6*l, r)]),
  ]))

function collide_with_body(particle::Particle, helicoid::IsotropicHelicoid)
  return collide_with_body(particle, helicoid.elements)
end

function get_initial_particle_state(mean_v, std_v, domain_radius)
  ϕ, θ = 2*π*rand(), acos(2*rand() - 1)
  x = domain_radius * SA[sin(θ)*cos(ϕ), sin(θ)*sin(ϕ), cos(θ)]
  v = std_v * SA[randn(), randn(), randn()] + mean_v
  Particle(0.0, x, v)
end

end # module GBGKnudsenGas