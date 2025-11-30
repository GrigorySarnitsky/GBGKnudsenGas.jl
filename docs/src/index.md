# GBGKnudsenGas.jl

This package simulates Knudsen flow of particles past a stationary body, where the particle mean free path is much larger than the body size. In this regime, particles interact with the body but not with each other. The body is assumed to be massive, so particle collisions do not move it.

The package provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). The package includes tools to compute the next collision of a single particle with a body; you can use these to build multi-particle, multi-collision simulations. A simple computational domain for modelling a gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a sphere is included.

All particle–body collisions are assumed to be specular. The package uses a Cartesian coordinate system. Note that `GBGKnudsenGas.jl` is purely kinematic: it works with distances, times, and velocities, but not masses or momenta.

The `GBG` part in the name `GBGKnudsenGas` comes from the city of Göteborg (Gothenburg, Sweden), where the primary development took place.

## Installation

This is an unregistered package, so it is not listed in the `General` Julia registry. To install it, run:

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

The package has three components: the gas particle (type `Particle`), the body past which the gas flows (type `Body`), and the function that calculates the collision between the two (`collide_with_body`). 

The state of a gas particle is represented by the type `Particle` with a constructor

```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

where `v` and `x` are the particle's velocity and position at time `t`. For a given particle `p`, these can be accessed as `p.v`, `p.x`, and `p.t` Note that we use `StaticArrays.jl` to represent the vectors `v` and `x`. For example, let's create a particle `p` that time $t = 0$ has position $(-10, 0, 0)$ and velocity $(1, 0, 0)$:

```julia
using StaticArrays
using GBGKnudsenGas

p = Particle(SA[1.0, 0.0, 0.0], SA[-10.0, 0.0, 0.0], 0.0)
```

Next, bodies immersed in the Knudsen gas are represented by an abstract supertype `Body`. A concrete example is `Sphere`, which represents a sphere centered at point `x` of radius `r`:

```julia
Sphere <: Body
Sphere(x::SVector{3, Float64}, r::Float64)
```

For a given sphere `s`, its center and radius can be accessed as `s.x` and `s.x`. For now, let's make a sphere `s` of radius $1$ and place it at the origin $(0.0, 0.0, 0.0)$ of our Cartesian coordinate system:

```julia
s = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Finally, there is a function `collide_with_body` which, given a `particle` and a `body`, returns the state of the particle immediately after the collision (at the surface of the body). If no collision occurs, it returns `nothing`.

```julia
collide_with_body(particle::Particle, body::Body) -> Union{Nothing, Particle}
```

Let's collide our particle `p` with the sphere `s`:

```julia
julia> collide_with_body(p, s)
Particle([-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 9.0)
```

As expected, the collision occurs at the surface of the sphere, reverses the first velocity component of the particle, and takes `9.0` units of time to happen.

Apart from spheres, `GBGKnudsenGas` provides a constructor for a cylinder without end caps (`CylinderSideWall`), as well as a way to combine simple bodies into a single body (`CompoundBody`). For these and the rest of the features, see the full API reference at the end.

## Complete example: Epstein drag on a sphere (currently wrong)

We provide a concrete example of how to use `GBGKnudsenGas.jl`. We consider the flow of a Knudsen gas past a sphere with specular collisions. Our goal is to calculate the per-collision change in particle velocity, $\langle \Delta \boldsymbol{v} \rangle_\text{c}$. Ensemble averages are denoted by angle brackets, $\langle {\cdot} \rangle$; the corresponding fluctuating parts are denoted with primes ${\cdot}'$. The index $\text{c}$ indicates averaging over the particles that actually collide with the sphere. In standard notation for conditional averages,

```math
\langle \Delta \boldsymbol{v} \rangle_\text{c} = \langle \boldsymbol{v}_\text{after} - \boldsymbol{v}_\text{before} | \text{particle has collided} \rangle,
```

where $`\boldsymbol{v}_\text{before}`$ and $`\boldsymbol{v}_\text{after}`$ are the particle's velocity before and after the collision with the sphere.

We want to check Epstein's 1924 calculation that $`\langle \Delta \boldsymbol{v} \rangle_\text{c}`$ is proportional to the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ at infinity:

```math
  \langle \Delta \boldsymbol{v} \rangle_\text{c} = - \delta \langle \boldsymbol{v} \rangle_\infty,
```

with $\delta = 1$ for specular collisions. In this form this relation can be obtained who knows how. See also the original paper: Epstein, Paul S. "On the resistance experienced by spheres in their motion through gases," _Physical Review,_ 23.6 (1924), 710.

Full code of this computation is available in `examples/epstein_sphere.jl`.

We start by defining our sphere, which we call `epstein_sphere`. We set its radius to $1.0$ and center it in the origin $(0, 0, 0)$ of our Cartesian coordinate system:

```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Next, we need to generate the initial state of the gas particles before they hit the sphere. The gas is at thermal equilibrium, moving past the sphere with the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ (denoted `mean_v` in code). Instead of temperature, we specify the standard deviation $`\sigma_{v, \infty}`$,

```math
 \sigma_{v, \infty} = \sqrt{\frac{1}{3} \langle \boldsymbol v' \cdot \boldsymbol v' \rangle_\infty},
```
 
 of particle velocities far away from the sphere (denoted `std_v` in code). Assuming a monotomic gas, the speed of sound is $c = \sqrt{5/3} \sigma_{v, \infty}$, and the Mach number $\text{Ma}$ of the flow is

```math
\text{Ma} = \sqrt{\frac{3}{5}} \frac{\| \langle \boldsymbol{v} \rangle_\infty \|}{\sigma_{v, \infty}}.
```

We initially place all particles on the boundary of a spherical domain of radius `domain_radius`, which should be sufficiently larger than the unit radius of `epstein_sphere`. Later we set `domain_radius` to `20`.

Thus defined initial state of a gas particle can be obtained with `get_initial_particle_state(mean_v, std_v, domain_radius)`. However, we are interested only in particles that will collide with the sphere. Given that `domain_radius` is much larger than the radius of `epstein_sphere`, most of the generated particles will miss the `epstein_sphere`. Hence we write a simple function `find_single_collision` that keeps drawing particles from the initial distribution until it finds one that collides with the `body`. Later we call `find_single_collision` with `epstein_sphere` for the `body`. The function returns the initial state of the particle `particle_before`, as well as its state after the collision `particle_after`:


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

To accurately estimate ensemble averages, we need a large sample of particles that collided with the body. The function `run_simulation` below finds `sample_size` particles (later we set `sample_size` to `10^6`) that collided with the `body` and returns two arrays, one with the initial state of the particles before the collision (`particles_before_collision`) and one with the particle states after the collision (`particles_after_collision`). The code is multithreaded to take advantage of a multi-core CPU:

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

Finally, we write a function to calculate $`\langle \Delta \boldsymbol{v} \rangle_\text{c}`$ for the specified Mach number `Ma`:

```julia
using Statistics

function mean_Δv_per_collision(body, Ma)
  std_v = 1.0
  mean_v = SA[Ma * std_v * √(5/3), 0.0, 0.0]
  domain_radius = 20.0
  sample_size = 10^6
  particles_before_collision, particles_after_collision = run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  return mean(p_a.v - p_b.v for (p_b, p_a) in zip(particles_before_collision, particles_after_collision))
end
```

We call it for our `epstein_sphere` and 10 different values of the Mach number, and calculate Epstein’s 
$\delta$:

```julia
Ma_s = range(0.0, 3.0, 10)
mean_Δv_s = @time mean_Δv_per_collision.(Ref(epstein_sphere), Ma_s)

for (Δv, Ma) in zip(mean_Δv_s, Ma_s)
  println("Ma = ", round(Ma, digits = 2), ",\t δ = ", -round(Δv[1] / (√(5/3) * Ma), digits = 5))
end
```

This should print something like

```
105.765674 seconds (758.63 k allocations: 1.079 GiB, 0.17% gc time, 1.02% compilation time)
Ma = 0.0,	 δ = Inf
Ma = 0.33,	 δ = 0.9967
Ma = 0.67,	 δ = 1.00221
Ma = 1.0,	 δ = 0.99963
Ma = 1.33,	 δ = 0.99864
Ma = 1.67,	 δ = 1.00079
Ma = 2.0,	 δ = 0.99864
Ma = 2.33,	 δ = 0.99994
Ma = 2.67,	 δ = 1.00015
Ma = 3.0,	 δ = 1.00037
```

As expected, $\delta \approx 1$ and is independent of the Mach number.

## Extending the library

To add another kind of body, you define a new body type (e.g. `Ellipsoid`) and implement the corresponding `collide_with_body` method:

```julia
collide_with_body(particle::Particle, elipsoid::Ellipsoid)::Union{Nothing, Particle}
```

The existing implementation for `Sphere` can be used for reference.

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

Currently, arbitrary meshes for bodies are not supported, but that would be a valuable feature to add.