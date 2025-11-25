# GBGKnudsenGas.jl

Simulate Knudsen flow of particles past a stationary body, where the particle mean free path is much larger than the body size. In this regime, particles interact with the body but not with each other. The body is assumed to be massive, so particle collisions do not move it.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a sphere is included.

All particle–body collisions are assumed to be specular. The package uses a Cartesian coordinate system. Note, `GBGKnudsenGas.jl` is purely kinematic: it works with distances, times, and velocities, but not masses or momenta.

The `GBG` part in the name `GBGKnudsenGas` comes from the city of Göteborg (Gothenburg, Sweden), where the primary development took place.

## Installation

This is an unregistered package, meaning that it is not listed in the `General` Julia registry and to install it you run

```julia
]add https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl
```

or

```julia
using Pkg
Pkg.add(url="https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl")
```

Then you can use it as usual:

```julia
using GBGKnudsenGas
```

## Using GBGKnudsenGas.jl

## Complete example: Epstein drag on a sphere

Now we provide a concrete example of using `GBGKnudsenGas.jl`. We consider the flow of a Knudsen gas past a sphere with specular collisions. Our goal is to calculate the per-collision change in particle velocity, $\langle \Delta \boldsymbol{v} \rangle_\text{c}$. Ensemble averages are denoted by angle brackets, $\langle {\cdot} \rangle$, primes ${\cdot}'$ denote fluctuating part with repsect to the average value. The index $\text{c}$ indicates averaging over the particles that actually collide with the sphere. In standard notation for conditional averages,

```math
\langle \Delta \boldsymbol{v} \rangle_\text{c} = \langle \boldsymbol{v}_\text{after} - \boldsymbol{v}_\text{before} | \text{particle has collieded} \rangle,
```

where $\boldsymbol{v}_\text{before}$ and $\boldsymbol{v}_\text{after}$ are the particle's velocity before and after the collision with the sphere.

We want to check the Eprstein's 1924 calculation that $\langle \Delta \boldsymbol{v} \rangle_\text{c}$ is proportional to the mean flow velocity $\langle \boldsymbol{v} \rangle_\infty$ at infinity:

```math
  \langle \Delta \boldsymbol{v} \rangle_\text{c} = - \delta \langle \boldsymbol{v} \rangle_\infty,
```

with $\delta = 1$ for specular collisions. In this form this relatin can be found in ??? 2003. See also the original paper, Epstein, Paul S. "On the resistance experienced by spheres in their motion through gases." _Physical Review,_ 23.6 (1924), 710.

Full code of this computation is available in `examples/epstein_sphere.jl`, where it is written for optimal performance. Below is a more explainatory version of 

We start by defining our sphere of radius $1.0$, centered at the origin $(0, 0, 0)$ of our Cartesian coordinate system:

```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Next, we need to generate the initial state of the gas particles before they hit the sphere. The gas is at thermal equilibrium, moving past the sphere with the mean flow velocity $\langle \boldsymbol{v} \rangle_\infty$ (denoted `mean_v` in code). Instead of tempreture, we specify the standrad deviation $\sqrt{\langle v'^2 \rangle}$ of particle velocities at equilibrium (`std_v`). We assume that initially all the particles start on the boundary of a spherical domain of radius `domain_radius`.


```julia
function find_single_collision(body, mean_v, std_v, domain_radius)
  while true
    particle_before = get_initial_particle_state(mean_v, std_v, domain_radius)
    particle_after  = collide_with_body(particle_before, body)
    if !isnothing(particle_after)
      return particle_before, particle_after
    end
  end 
end
```

This will give us . But to calculate , we need 

```julia
using Base.Threads

function run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  particles_before_collision = Array{Particle, 1}(undef, sample_size)
  particles_after_collision  = Array{Particle, 1}(undef, sample_size)
  @threads for i in 1:sample_size
    particles_before_collision[i], particles_after_collision[i] = find_single_collision(body, mean_v, std_v, domain_radius)
  end
  return particles_before_collision, particles_after_collision
end
```


## Extending the library

To add another kind of body, one needs to define the body-type, e.g. `Elipsoid` and implement the respective `collide_with_body` function,

```julia
collide_with_body(particle::Particle, elipsoid::Elipsoid)::Union{Nothing, Particle}
```

Existing implementation for `Sphere` can be used for reference.

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

Currently, I haven't thought of adding arbitrary meshes for bodies, but that would be a valuable feature to have.

## API Reference

Below is the API reference generated from the docstrings

### Modules

```@docs
GBGKnudsenGas
```

### Types

```@docs
Particle
Body
Sphere
CylinderSideWall
CompoundBody
CompoundBodyDynamic
IsotropicHelicoid
```

### Functions

```@docs
collide_with_body
get_initial_particle_state
```