
<a id='GBGKnudsenGas.jl'></a>

<a id='GBGKnudsenGas.jl-1'></a>

# GBGKnudsenGas.jl


Simulate Knudsen flow of particles past a stationary body, where the particle mean free path is much larger than the body size. In this regime, particles interact with the body but not with each other. The body is assumed to be massive, so particle collisions do not move it.


The module provides building blocks to create your own flow setups. The body is built from simple geometries (currently spheres and cylinders). Tools are included to compute the next collision of a single particle with a body, which you can use to build multi-particle, multi-collision simulations. A simple computational domain for modeling gas in thermal equilibrium flowing past a body is also provided.


An example demonstrating a complete flow simulation around a sphere is included.


All particle–body collisions are assumed to be specular. The package uses a Cartesian coordinate system. Note, `GBGKnudsenGas.jl` is purely kinematic: it works with distances, times, and velocities, but not masses or momenta.


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


<a id='Complete-example:-Epstein-drag-on-a-sphere'></a>

<a id='Complete-example:-Epstein-drag-on-a-sphere-1'></a>

## Complete example: Epstein drag on a sphere


Now we provide a concrete example of using `GBGKnudsenGas.jl`. We consider the flow of a Knudsen gas past a sphere with specular collisions. Our goal is to calculate the per-collision change in particle velocity, $\langle \Delta \boldsymbol{v} \rangle_\text{c}$. Ensemble averages are denoted by angle brackets, $\langle {\cdot} \rangle$, primes ${\cdot}'$ denote fluctuating part with repsect to the average value. The index $\text{c}$ indicates averaging over the particles that actually collide with the sphere. In standard notation for conditional averages,


$$
\langle \Delta \boldsymbol{v} \rangle_\text{c} = \langle \boldsymbol{v}_\text{after} - \boldsymbol{v}_\text{before} | \text{particle has collided} \rangle,
$$


where $`\boldsymbol{v}_\text{before}`$ and $`\boldsymbol{v}_\text{after}`$ are the particle's velocity before and after the collision with the sphere.


We want to check the Eprstein's 1924 calculation that $`\langle \Delta \boldsymbol{v} \rangle_\text{c}`$ is proportional to the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ at infinity:


$$
  \langle \Delta \boldsymbol{v} \rangle_\text{c} = - \delta \langle \boldsymbol{v} \rangle_\infty,
$$


with $\delta = 1$ for specular collisions. In this form this relatin can be found in ??? 2003. See also the original paper, Epstein, Paul S. "On the resistance experienced by spheres in their motion through gases." *Physical Review,* 23.6 (1924), 710.


Full code of this computation is available in `examples/epstein_sphere.jl`, where it is written for optimal performance. Below is a more explainatory version of 


We start by defining our sphere, which we call `epstein_sphere`. We set its radius to $1.0$ and center it the origin $(0, 0, 0)$ of our Cartesian coordinate system:


```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```


Next, we need to generate the initial state of the gas particles before they hit the sphere. The gas is at thermal equilibrium, moving past the sphere with the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ (denoted `mean_v` in code). Instead of tempreture, we specify the standrad deviation $`\sigma_{v, \infty}`$,


of particle velocities far away from the sphere (denoted `std_v` in code). Assuming monotomic gas particles, the Mach number $\text{Ma}$ of the flow is


$$
\text{Ma} = \frac{3}{4} \frac{\| \langle \boldsymbol{v} \rangle_\infty \| }{\sigma_v}
$$


We put initially all the particles on the boundary of a spherical domain of radius `domain_radius` which should be sufficiently larger than the unit radius of `epstein_sphere`.


Thus defined initial state of a gas particle can be obtained with `get_initial_particle_state(mean_v, std_v, domain_radius)`. However we are interested only in particles that will collide with the sphere. Given that `domain_radius` is much larger than the radius of `epstein_sphere`, most of the generated particles will miss the `epstein_sphere`. Hence we write a simple function `find_single_collision` that keeps drawing particles from the initial distribution untill it finds one that collides with the `body`. Later we call `find_single_collision` with `epstein_sphere` for the `body`. The function returns the initial state of the particle `particle_before`, as well as its state after the collision `particle_after`:


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


To accuratly estimate ensemble averages, we need a large sample . The function `run_simulation` below finds `sample_size` particles that collided with the `body` and returns two arrays, one with the intitial state of the particles before the collision (`particles_before_collision`) and one with the particle states after the collision (`particles_after_collision`). The code is multithreaded to take advantage of a multi-core CPU:


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


Finally, we want to 


```julia
using Statistics

function mean_Δv_per_collision(body, Ma)
  std_v = 1.0
  mean_v = SA[Ma * std_v * √(8/π), 0.0, 0.0]
  domain_radius = 20.0
  sample_size = 10^6
  
  particles_before_collision, particles_after_collision = run_simulation(body, mean_v, std_v, domain_radius, sample_size)

  return mean(p_a.v - p_b.v for (p_b, p_a) in zip(particles_before_collision, particles_after_collision))
end
```


And call it for `epstein_sphere` and 10 different values of the Mach number:


```julia
Ma_s = 0.1:1:10.0
means_Δv_per_collision = @time mean_Δv_per_collision.(Ma_s, Ref(epstein_sphere))
println(means_Δv_per_collision)
```


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

Currently, all particle–body collisions are assumed to be specular. Throughout the package, we use a Cartesian coordinate system.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L1-L11' class='documenter-source'>source</a><br>


<a id='Types'></a>

<a id='Types-1'></a>

### Types

<a id='GBGKnudsenGas.Particle' href='#GBGKnudsenGas.Particle'>#</a>
**`GBGKnudsenGas.Particle`** &mdash; *Type*.



```julia
Particle(v::SVector{3, Float64}, x::SVector{3, Float64}, t::Float64)
```

Represent a state of a single point-particle (velocity vector `v` and position vector `x`) at time `t`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L23-L27' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Body' href='#GBGKnudsenGas.Body'>#</a>
**`GBGKnudsenGas.Body`** &mdash; *Type*.



```julia
Body
```

Abstract supertype for all bodies that can interact with a `Particle` via `collide_with_body(::Particle, ::Body)`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L35-L39' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.Sphere' href='#GBGKnudsenGas.Sphere'>#</a>
**`GBGKnudsenGas.Sphere`** &mdash; *Type*.



```julia
Sphere <: Body
Sphere(x::SVector{3, Float64}, r::Float64)
```

Sphere of radius `r` centered at point `x`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L53-L58' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CylinderSideWall' href='#GBGKnudsenGas.CylinderSideWall'>#</a>
**`GBGKnudsenGas.CylinderSideWall`** &mdash; *Type*.



```julia
CylinderSideWall <: Body
CylinderSideWall(x_1::SVector{3, Float64}, x_2::SVector{3, Float64}, r::Float64)
```

Finite circular cylinder defined by its axis endpoints `x_1` and `x_2` and radius `r`.

Only the lateral surface (side wall) is present; end caps are not included and do not participate in collisions. Endcaps can be added via `CompoundBody`, like we do in `IsotropicHelicopiod`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L82-L89' class='documenter-source'>source</a><br>

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


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L146-L160' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.CompoundBodyDynamic' href='#GBGKnudsenGas.CompoundBodyDynamic'>#</a>
**`GBGKnudsenGas.CompoundBodyDynamic`** &mdash; *Type*.



```julia
CompoundBodyDynamic <: Body
CompoundBodyDynamic(bodies::Array{Body, 1})
```

A compound body whose components are stored in the array `bodies`. A collision with a `CompoundBodyDynamic` is defined as the earliest collision with any of its components.

This representation might be preferable if the number of component bodies is very large, however the default choice should be tuple-based `CompoundBody`. Performance differences between `CompoundBodyDynamic` and `CompoundBody` should be evaluated on a case-by-case basis.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L121-L128' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.IsotropicHelicoid' href='#GBGKnudsenGas.IsotropicHelicoid'>#</a>
**`GBGKnudsenGas.IsotropicHelicoid`** &mdash; *Type*.



```julia
IsotropicHelicoid <: Body
IsotropicHelicoid(cube_side_length::Float64, r::Float64)
```

A model of Lord Kelvin's isotropic helicoid, constructed from 12 rods. Each rod is a cylinder of radius `r`, capped with two spheres of the same radius. The helicoid is inscribed in a cube of side length `cube_side_length`. This cube side length indirectly defines the length of each rod.

The rods that form the helicoid can be accessed via `helicoid.elements`, which is a `CompoundBody`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L190-L197' class='documenter-source'>source</a><br>


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


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L43-L49' class='documenter-source'>source</a><br>

<a id='GBGKnudsenGas.get_initial_particle_state' href='#GBGKnudsenGas.get_initial_particle_state'>#</a>
**`GBGKnudsenGas.get_initial_particle_state`** &mdash; *Function*.



```julia
get_initial_particle_state(mean_v, std_v, domain_radius) -> Particle
```

Sample an initial particle state for a gas in thermal equilibrium flowing past a body:

  * The initial position `Particle.x` is drawn uniformly on the sphere of radius `domain_radius` centered at the origin `[0.0, 0.0, 0.0]`.
  * The velocity `Particle.v` is drawn from a 3D normal distribution with isotropic standard deviation `std_v` and mean drift `mean_v`.
  * The time `Particle.t` is set to `0.0`.


<a target='_blank' href='https://github.com/GrigorySarnitsky/GBGKnudsenGas.jl/blob/ea4d72434968a26383947202372ff06b944735fa/src/GBGKnudsenGas.jl#L251-L259' class='documenter-source'>source</a><br>

