# GBGKnudsenGas.jl

Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a special body (an isotropic helicoid) is included.

Currently, all particle–body collisions are assumed to be specular.

The name `GBGKnudsenGas` comes from the city of Göteborg (Gothenburg, Sweden), where its primary development took place.

## Installation

This is an unregistered package, meaning that it is not listed in the `General` Julia registry and to install it you run

    ]add https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl

or

    using Pkg
    Pkg.add(url="hhttps://github.com/GrigorySarnitsky/GBGKnudsenGas.jl")

Then you can use it as usual:

    using GBGKnudsenGas

## API Reference

Below is the API reference generated from the docstrings

```@autodocs
Modules = [GBGKnudsenGas]
```