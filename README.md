
<a id='GBGKnudsenGas.jl'></a>

<a id='GBGKnudsenGas.jl-1'></a>

# GBGKnudsenGas.jl


Simulate Knudsen flow of particles past a body, in the Knudsen regime.


<a id='API'></a>

<a id='API-1'></a>

## API

<a id='GBGKnudsenGas.GBGKnudsenGas' href='#GBGKnudsenGas.GBGKnudsenGas'>#</a>
**`GBGKnudsenGas.GBGKnudsenGas`** &mdash; *Module*.



```julia
GBGKnudsenGas
```

Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a special body (an isotropic helicoid) is included.

Currently, all particle–body collisions are assumed to be specular.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/9d737f4cd1b18def253d15c41166bceba5550a6b/src/GBGKnudsenGas.jl#L1-L11' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Particle' href='#GBGKnudsenGas.Particle'>#</a>
**`GBGKnudsenGas.Particle`** &mdash; *Type*.



```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/9d737f4cd1b18def253d15c41166bceba5550a6b/src/GBGKnudsenGas.jl#L20-L24' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Sphere' href='#GBGKnudsenGas.Sphere'>#</a>
**`GBGKnudsenGas.Sphere`** &mdash; *Type*.



Sphere(x::SVector{3, Float64}, x::Float64)

Sphere centered at position `x` of radius `r`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/9d737f4cd1b18def253d15c41166bceba5550a6b/src/GBGKnudsenGas.jl#L35-L39' class='documenter-source'>source</a><br>

