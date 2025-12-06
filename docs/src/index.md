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

At its core, the package provides functionality to calculate collisions between a single gas particle and a given body. Users are expected to build their own ensembles to model the flow of interest. First, we describe how to use `GBGKnudsenGas.jl` to calculate the dynamics of a single particle; then we demonstrate how to compute ensemble averages for the special case of a uniform gas flow in thermal equilibrium. In addition, a [separate page](Epstein_drag.md) illustrates the calculation of Epstein drag on a sphere.

### Single particle dynamics

The package describes single-particle dynamics through three components: a gas particle (type `Particle`), a body with which the particle interacts (type `Body`), and a function that computes their collision (`collide_with_body`).

The state of a gas particle is represented by the type `Particle` with a constructor

```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

where `v` and `x` are the particle's velocity and position at time `t`. For a given particle `p`, these can be accessed as `p.v`, `p.x`, and `p.t` Note that we use `StaticArrays.jl` to represent the vectors `v` and `x`. For example, let's create a particle `p` that at time $t = 0$ has velocity $(1, 0, 0)$ and position $(-10, 0, 0)$:

```julia
using StaticArrays
using GBGKnudsenGas

p = Particle(SA[1.0, 0.0, 0.0], SA[-10.0, 0.0, 0.0], 0.0)
```

Next, bodies immersed in the Knudsen gas are represented by an abstract supertype `Body`. A concrete example is `Sphere`, which represents a sphere centered at point `x` of radius `r`:

```julia
Sphere(x::SVector{3, Float64}, r::Float64)
```

For a given sphere `s`, its center and radius can be accessed as `s.x` and `s.r`. For now, let's make a sphere `s` of radius $1$ and place it at the origin $(0.0, 0.0, 0.0)$ of our Cartesian coordinate system:

```julia
s = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Finally, the function `collide_with_body` takes a `particle` and a `body`, and returns the state of the particle immediately after its collision with the body (at the surface of the body). If no collision occurs, it returns `nothing`.

```julia
collide_with_body(particle::Particle, body::Body) -> Union{Nothing, Particle}
```

Let's collide our particle `p` with the sphere `s`:

```julia
julia> collide_with_body(p, s)
Particle([-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 9.0)
```

As expected, the collision occurs at the surface of the sphere, reverses the first component of the particle’s velocity, and takes 9.0 units of time.

Apart from spheres, `GBGKnudsenGas` also provides `CylinderSideWall`, which represents a cylinder without end caps, and `CompoundBody`,` which allows multiple simple bodies to be combined into a single composite body.  For details, see the [full API reference.](api.md)


### Simulating flow of gas in thermal equilibrium past a body

Consider a flow of Knudsen gas past a body with specular collisions. The gas is at thermal equilibrium with unperturbed (at infinity) temperature $T$, density $\rho$ and mean velocity $\langle \boldsymbol v_0 \rangle$. Angle brackets $\langle {\cdot} \rangle$ denote ensemble averages, and primes ${\cdot}'$ denote the respective fluctuating parts. Before the collision with the body, the initial particle velocity $\boldsymbol v_0$ has the Maxwell-Boltzmann probability density function $f(\boldsymbol v_0)$:

```math
f(\boldsymbol v_0) = \sqrt{\frac{1}{}} \exp \frac{\boldsymbol v_0' \cdot \boldsymbol v_0'}{\sigma_v^2}
```

with

```math
\sigma_v = 
```

Consider force $\boldsymbol{F}$ that the flowing gas exerts on the body. It can be calculated as

```math
\boldsymbol{F} = \int\limits_{\mathclap{\substack{\text{velocity} \\ \text{space}}}} \biggl[\int\limits_{D(\boldsymbol v)} \rho \Delta \boldsymbol v(\boldsymbol x_0) \, \mathrm{d} \boldsymbol x_0 \biggr] \| \boldsymbol v_0 \| f(\boldsymbol v_0) \, \mathrm{d} \boldsymbol v_0
```

![integration scheme](integration_scheme.svg)

```math
\int\limits_{\mathclap{\substack{\text{velocity} \\ \text{space}}}} \biggl[\int\limits_{D(\boldsymbol v)} H(\boldsymbol v, \boldsymbol x) \, \mathrm{d} \boldsymbol x \biggr] \| \boldsymbol v \| f(\boldsymbol v) \, \mathrm{d} \boldsymbol v
```

### Complete example: Epstein drag

We provide a concrete example of how to use `GBGKnudsenGas.jl` by calculate the drag excerted on a sphere by a uniform flow of Knudsen gas at thermal equilibrium, see the [respective page.](Epstein_drag.md)


## Adding new kinds of bodies

To add another kind of body, you define a new body type (e.g. `Ellipsoid`) and implement the corresponding `collide_with_body` method:

```julia
collide_with_body(particle::Particle, elipsoid::Ellipsoid)::Union{Nothing, Particle}
```

The existing implementation for `Sphere` can be used for reference.

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

Currently, arbitrary meshes for bodies are not supported, but that would be a valuable feature to add.