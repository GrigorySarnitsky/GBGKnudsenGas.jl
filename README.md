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

## Documentation

Documentation with the full API reference is available from [https://grigorysarnitsky.github.io/GBGKnudsenGas.jl/](https://grigorysarnitsky.github.io/GBGKnudsenGas.jl/).

