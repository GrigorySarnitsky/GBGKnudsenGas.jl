
<a id='GBGKnudsenGas.jl'></a>

<a id='GBGKnudsenGas.jl-1'></a>

# GBGKnudsenGas.jl

- [GBGKnudsenGas.jl](index.md#GBGKnudsenGas.jl)
    - [Installation](index.md#Installation)
    - [Using GBGKnudsenGas.jl](index.md#Using-GBGKnudsenGas.jl)
    - [Complete example with visualization](index.md#Complete-example-with-visualization)
    - [Extending the library](index.md#Extending-the-library)
    - [API Reference](index.md#API-Reference)


Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.


The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.


An example demonstrating a complete flow simulation around a special body (isotropic helicoid) is included.


Currently, all particle–body collisions are assumed to be specular.


The `GBG` part in the name `GBGKnudsenGas` comes from the city of Göteborg (Gothenburg, Sweden), where the primary development took place.


<a id='Installation'></a>

<a id='Installation-1'></a>

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


<a id='Using-GBGKnudsenGas.jl'></a>

<a id='Using-GBGKnudsenGas.jl-1'></a>

## Using GBGKnudsenGas.jl


<a id='Complete-example-with-visualization'></a>

<a id='Complete-example-with-visualization-1'></a>

## Complete example with visualization


<a id='Extending-the-library'></a>

<a id='Extending-the-library-1'></a>

## Extending the library


To add another kind of body, one needs to define the body-type, e.g. `Elipsoid` and implement the respective `collide_with_body` function,


```julia
collide_with_body(particle::Particle, elipsoid::Elipsoid)::Union{Nothing, Particle}
```


Existing implementation for `Sphere` can be used for reference.


As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.


Currently, I haven't thought of adding arbitrary meshes for bodies, but that would be a valuable feature to have.


<a id='API-Reference'></a>

<a id='API-Reference-1'></a>

## API Reference


Below is the API reference generated from the docstrings


<a id='Modules'></a>

<a id='Modules-1'></a>

### Modules

<a id='GBGKnudsenGas' href='#GBGKnudsenGas'>#</a>
**`GBGKnudsenGas`** &mdash; *Module*.



```julia
GBGKnudsenGas
```

Simulate Knudsen flow of particles past a body, where the particle mean free path is much larger than the body size. In this regime, particles interact only with the body and not with each other.

The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.

An example demonstrating a complete flow simulation around a special body (isotropic helicoid) is included.

Currently, all particle–body collisions are assumed to be specular.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L1-L11' class='documenter-source'>source</a><br>


<a id='Types'></a>

<a id='Types-1'></a>

### Types

<a id='GBGKnudsenGas.Particle' href='#GBGKnudsenGas.Particle'>#</a>
**`GBGKnudsenGas.Particle`** &mdash; *Type*.



```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L23-L27' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Body' href='#GBGKnudsenGas.Body'>#</a>
**`GBGKnudsenGas.Body`** &mdash; *Type*.



```julia
Body
```

Abstract supertype for all bodies that can interact with a `Particle` via `collide_with_body(::Particle, ::Body)`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L35-L39' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Sphere' href='#GBGKnudsenGas.Sphere'>#</a>
**`GBGKnudsenGas.Sphere`** &mdash; *Type*.



```julia
Sphere <: Body
Sphere(x::SVector{3, Float64}, x::Float64)
```

Sphere of radius `r` centered at point `x`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L53-L58' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CylinderSideWall' href='#GBGKnudsenGas.CylinderSideWall'>#</a>
**`GBGKnudsenGas.CylinderSideWall`** &mdash; *Type*.



```julia
CylinderSideWall <: Body
CylinderSideWall(x_1::SVector{3, Float64}, x_2::SVector{3, Float64}, r::Float64)
```

Finite circular cylinder defined by its axis endpoints `x_1` and `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. Endcaps can be added via `CompoundBody`, like we do in `IsotropicHelicopiod`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L82-L89' class='documenter-source'>source</a><br>


<a id='Functions'></a>

<a id='Functions-1'></a>

### Functions

<a id='GBGKnudsenGas.collide_with_body' href='#GBGKnudsenGas.collide_with_body'>#</a>
**`GBGKnudsenGas.collide_with_body`** &mdash; *Function*.



```julia
collide_with_body(particle::Particle, body::Body) -> Union{Nothing,Particle}
```

Compute the next collision of `particle` with `body`, assuming straight-line motion between collisions.

Return a new `Particle` located at the collision point, with updated velocity (after collision) and time, or `nothing` if there is no collision.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/0b4777f08b730d5f93f234971aede6a956832713/src/GBGKnudsenGas.jl#L43-L49' class='documenter-source'>source</a><br>

