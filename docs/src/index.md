# GBGKnudsenGas.jl

Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a special body (isotropic helicoid) is included.

Currently, all particle–body collisions are assumed to be specular. Throughout the package, we use a Cartesian coordinate system.

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

## Complete example: Epstein drag

Now we provide a concrete example of using `GBGKnudsenGas.jl`. We consider the flow of a Knudsen gas past a sphere with specular collisions. Our goal is to calculate the per-collision change in particle velocity, $\langle \Delta v \rangle_\text{c}$. Ensemble averages are denoted by angle brackets, $\langle {\cdot} \rangle$, the index $\text{c}$ indicates averaging only over the particles that actually collide with the sphere. In standard notations for conditional averages,

```math
\langle \Delta v \rangle_\text{c} = \langle v_\text{after} - v_\text{before} | \text{particle has collieded} \rangle.
```

Full code of this computation is available in `examples/epstein_sphere.jl`, where it is written for optimal performance. Below is a more explainatory version of 

We start by defining our sphere of radius $1.0$, centered at the origin $(0, 0, 0)$ of our Cartesian coordinate system:

```julia
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Gas particles originate on 


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