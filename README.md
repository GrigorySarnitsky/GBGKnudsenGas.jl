
<a id='GBGKnudsenGas.jl'></a>

<a id='GBGKnudsenGas.jl-1'></a>

# GBGKnudsenGas.jl


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


<a id='Complete-example:-Epstein-drag'></a>

<a id='Complete-example:-Epstein-drag-1'></a>

## Complete example: Epstein drag


Now we provide a concrete example of using `GBGKnudsenGas.jl`. We consider a flow of Knudsen gas around a sphere with specular collisions. We set our goal to calculate the per-collision change of particle velocity $\langle \Delta v \rangle_\text{c}$:


:$


\langle \Delta v \rangle*\text{c} = \langle v*\text{after collision} - v_\text{before collision} | \text{particle has collieded} \rangle :$


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


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L1-L11' class='documenter-source'>source</a><br>


<a id='Types'></a>

<a id='Types-1'></a>

### Types

<a id='GBGKnudsenGas.Particle' href='#GBGKnudsenGas.Particle'>#</a>
**`GBGKnudsenGas.Particle`** &mdash; *Type*.



```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L23-L27' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Body' href='#GBGKnudsenGas.Body'>#</a>
**`GBGKnudsenGas.Body`** &mdash; *Type*.



```julia
Body
```

Abstract supertype for all bodies that can interact with a `Particle` via `collide_with_body(::Particle, ::Body)`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L35-L39' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Sphere' href='#GBGKnudsenGas.Sphere'>#</a>
**`GBGKnudsenGas.Sphere`** &mdash; *Type*.



```julia
Sphere <: Body
Sphere(x::SVector{3, Float64}, r::Float64)
```

Sphere of radius `r` centered at point `x`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L53-L58' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CylinderSideWall' href='#GBGKnudsenGas.CylinderSideWall'>#</a>
**`GBGKnudsenGas.CylinderSideWall`** &mdash; *Type*.



```julia
CylinderSideWall <: Body
CylinderSideWall(x_1::SVector{3, Float64}, x_2::SVector{3, Float64}, r::Float64)
```

Finite circular cylinder defined by its axis endpoints `x_1` and `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. Endcaps can be added via `CompoundBody`, like we do in `IsotropicHelicopiod`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L82-L89' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CompoundBody' href='#GBGKnudsenGas.CompoundBody'>#</a>
**`GBGKnudsenGas.CompoundBody`** &mdash; *Type*.



```julia
CompoundBody <: Body
CompoundBody{T<:Tuple{Vararg{<:Body}}}
```

A compound body whose components are stored in the tuple `bodies`. A collision with a `CompoundBody` is defined as the earliest collision with any of its components.

This type should be the default choice for representing compound bodies, as it is usually faster.

Use the convenience constructors

```
CompoundBody(bodies::Vararg{<:Body})
CompoundBody(bodies::AbstractVector{<:Body})
```

to build instances.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L146-L160' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CompoundBodyDynamic' href='#GBGKnudsenGas.CompoundBodyDynamic'>#</a>
**`GBGKnudsenGas.CompoundBodyDynamic`** &mdash; *Type*.



```julia
CompoundBodyDynamic <: Body
CompoundBodyDynamic(bodies::Array{Body, 1})
```

A compound body whose components are stored in the array `bodies`. A collision with a `CompoundBodyDynamic` is defined as the earliest collision with any of its components.

This representation might be preferable if the number of component bodies is very large, however the default choice should be tuple-based `CompoundBody`. Performance differences between `CompoundBodyDynamic` and `CompoundBody` should be evaluated on a case-by-case basis.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L121-L128' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.IsotropicHelicoid' href='#GBGKnudsenGas.IsotropicHelicoid'>#</a>
**`GBGKnudsenGas.IsotropicHelicoid`** &mdash; *Type*.



```julia
IsotropicHelicoid <: Body
IsotropicHelicoid(cube_side_length::Float64, r::Float64)
```

A model of Lord Kelvin's isotropic helicoid, constructed from 12 rods. Each rod is a cylinder of radius `r`, capped with two spheres of the same radius. The helicoid is inscribed in a cube of side length `cube_side_length`. This cube side length indirectly defines the length of each rod.

The rods that form the helicoid can be accessed via `helicoid.elements`, which is a `CompoundBody`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L190-L197' class='documenter-source'>source</a><br>


<a id='Functions'></a>

<a id='Functions-1'></a>

### Functions

<a id='GBGKnudsenGas.collide_with_body' href='#GBGKnudsenGas.collide_with_body'>#</a>
**`GBGKnudsenGas.collide_with_body`** &mdash; *Function*.



```julia
collide_with_body(particle::Particle, body::Body) -> Union{Nothing, Particle}
```

Compute the next collision of `particle` with `body`, assuming straight-line motion between collisions.

Return a new `Particle` located at the collision point, with updated velocity (after the collision) and time, or `nothing` if there is no collision.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L43-L49' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.get_initial_particle_state' href='#GBGKnudsenGas.get_initial_particle_state'>#</a>
**`GBGKnudsenGas.get_initial_particle_state`** &mdash; *Function*.



```julia
get_initial_particle_state(mean_v, std_v, domain_radius) -> Particle
```

Sample an initial particle state for a gas in thermal equilibrium flowing past a body:

  * The initial position `Particle.x` is drawn uniformly on the sphere of radius `domain_radius` centered at the origin `[0.0, 0.0, 0.0]`.
  * The velocity `Particle.v` is drawn from a 3D normal distribution with isotropic standard deviation `std_v` and mean drift `mean_v`.
  * The time `Particle.t` is set to `0.0`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ee39465ba20d07fca967f905fc5d3df7905ee5f1/src/GBGKnudsenGas.jl#L251-L259' class='documenter-source'>source</a><br>

