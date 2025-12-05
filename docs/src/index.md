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

### Basic components

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

### Uniform flow of gas at thermal equilibrium

Figuring out the integral

```math
f_\infty(\boldsymbol v_0) = \exp \frac{\boldsymbol v' \cdot \boldsymbol v'}{}
```

Consider force $\boldsymbol{F}$ that the flowing gas exerts on the body. It can be calculated as

```math
\boldsymbol{F} = \int\limits_{\mathclap{\substack{\text{velocity} \\ \text{space}}}} \biggl[\int\limits_{D(\boldsymbol v)} \rho \Delta \boldsymbol v(\boldsymbol x_0) \, \mathrm{d} \boldsymbol x_0 \biggr] \| \boldsymbol v_0 \| f(\boldsymbol v_0) \, \mathrm{d} \boldsymbol v_0
```

![integration scheme](integration_scheme.svg)

```math
\int\limits_{\mathclap{\substack{\text{velocity} \\ \text{space}}}} \biggl[\int\limits_{D(\boldsymbol v)} H(\boldsymbol v, \boldsymbol x) \, \mathrm{d} \boldsymbol x \biggr] \| \boldsymbol v \| f(\boldsymbol v) \, \mathrm{d} \boldsymbol v
```

## Complete example: Epstein drag

We provide a concrete example of how to use `GBGKnudsenGas.jl` by calculate the drag excerted on a sphere by a uniform flow of Knudsen gas at thermal equilibrium, see the [respective page.](Epstein_drag.md)


## Extending the library

To add another kind of body, you define a new body type (e.g. `Ellipsoid`) and implement the corresponding `collide_with_body` method:

```julia
collide_with_body(particle::Particle, elipsoid::Ellipsoid)::Union{Nothing, Particle}
```

The existing implementation for `Sphere` can be used for reference.

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

Currently, arbitrary meshes for bodies are not supported, but that would be a valuable feature to add.