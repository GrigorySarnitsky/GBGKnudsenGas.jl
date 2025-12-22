# Using GBGKnudsenGas.jl

At its core, the package provides functionality to calculate collisions between a single gas particle and a given body. Users are expected to build their own ensembles to model the flow of interest. First, we describe how to use `GBGKnudsenGas.jl` to calculate the dynamics of a single particle; then we demonstrate how to compute ensemble averages for the special case of a uniform gas flow in thermal equilibrium. Finally, there is a [dedicated page](Epstein_drag.md) illustrating how to build your own flow setup using Epstein drag calculation as an example.


## Single particle dynamics

The package describes single-particle dynamics through three components: a gas particle (type `Particle`), a body with which the particle interacts (type `Body`), and a function that computes their collision (`collide_with_body`).

The state of a gas particle is represented by the type `Particle` with a constructor

```julia
Particle(v::SVector{3, T}, x::SVector{3, T}, t::T) where {T <: Real}
```

where `v` and `x` are the particle's velocity and position at time `t`. For a given particle `p`, these can be accessed as `p.v`, `p.x`, and `p.t` Note that we use `StaticArrays.jl` to represent the vectors `v` and `x`. For example, let's create a particle `p` that at time $t = 0$ has velocity $(1, 0, 0)$ and position $(-10, 0, 0)$:

```julia
using StaticArrays
using GBGKnudsenGas

p = Particle(SA[1.0, 0.0, 0.0], SA[-10.0, 0.0, 0.0], 0.0)
```

Next, bodies immersed in the Knudsen gas are represented by an abstract supertype `Body`. A concrete example is `Sphere`, which represents a sphere centered at point `x` of radius `r`:

```julia
Sphere(x::SVector{3, T_x}, r::T_r) where {T_x <: Real, T_r <: Real}
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

Apart from spheres, `GBGKnudsenGas.jl` also provides `CylinderSideWall`, which represents a cylinder without end caps, and `CompoundBody`, which allows multiple simple bodies to be combined into a single composite body. For details, see the [full API reference.](api.md)


## Simulating stationary flow of a gas in thermal equilibrium past a body

Consider stationary flow of a Knudsen gas past a body with specular collisions. The gas is at thermal equilibrium, with unperturbed (at infinity) temperature $T$, density $\rho$, and mean velocity $\langle \boldsymbol v_0 \rangle$. Angle brackets $\langle {\cdot} \rangle$ denote ensemble averages, and primes ${\cdot}'$ denote the corresponding fluctuations. Before colliding with the body, the initial particle velocity $ \boldsymbol v_0$ follows the Maxwell–Boltzmann probability density function

```math
f(\boldsymbol v_0) = \frac{1}{(\sqrt{2 \pi \sigma_v^2})^3} \exp \biggl[- \frac{(\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle) \cdot (\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle)}{2 \sigma_v^2}\biggr],
```

where

```math
\sigma_v = \sqrt{\frac{k_\text{B} T}{m}},
```

and $m$ is the mass of a single particle.

We now describe a procedure for calculating physical quantities associated with the flow of a Knudsen gas past a body. We use the force $\boldsymbol{F}$ that the flowing gas exerts on the body as an archetypical example. Assume the body fits entirely within a spherical domain of radius $r_\text{dom}$. Let $\boldsymbol v_\infty(\boldsymbol v_0, \boldsymbol x_0)$ be the final velocity of a gas particle that initially had velocity $\boldsymbol v_0$ and position $\boldsymbol x_0$. The particle's velocity increment is $\Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) = \boldsymbol v_\infty - \boldsymbol v_0$. We can express $\boldsymbol{F}$ as an integral over in the velocity and physical spaces:

```math
\boldsymbol{F} = - \int\limits_{\mathclap{\substack{\text{velocity} \\ \text{space}}}} \biggl[\int\limits_{D(\boldsymbol v_0)} \rho  \| \boldsymbol v_0 \| \Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) \, \mathrm{d} \boldsymbol x_0 \biggr] f(\boldsymbol v_0) \, \mathrm{d} \boldsymbol v_0.
```
For a given $\boldsymbol v_0$, the inner integration is taken over the disk $D(\boldsymbol v_0)$ of radius  $r_\text{dom}$, oriented perpendicular to $\boldsymbol v_0$ and located at a distance $r_\text{dom}$ from the domain center, as illustrated below:

![integration scheme](integration_scheme.svg)

We evaluate the integral for $\boldsymbol F$ using a Monte Carlo method. We sample particles whose initial velocities $\boldsymbol v_0$ follow the Maxwell–Boltzmann PDF $f(\boldsymbol v_0)$, and whose initial positions are distributed uniformly over the corresponding disk $D(\boldsymbol v_0)$. With thus defined ensemble of initial particle states, $\boldsymbol F$ is approximated as the following ensemble average:

```math
\boldsymbol F = - \rho  \pi r_\text{dom}^2 \bigl< \| \boldsymbol v_0 \| \Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) \bigr>.
```

To generate this ensemble, `GBGKnudsenGas.jl` provides the function [`sample_particle_state_uniform_flow`](@ref):

```julia
sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
```

Given $\langle \boldsymbol v_0 \rangle$ (`mean_v`), $\sigma_v$ (`std_v`) and $r_\text{dom}$ (`domain_radius`), the function returns a particle `p`, whose fields `p.v` and `p.x` follow the required distributions for $\boldsymbol v_0$ and $\boldsymbol x_0$, with `p.t` set to $0$. Further, there is a function [`find_colliding_particle_uniform_flow`](@ref),

```julia
find_colliding_particle_uniform_flow(body, mean_v, std_v, domain_radius)
```

that repeatedly samples particles via `sample_particle_state_uniform_flow` until it finds the one that collides with body.

Finally, `GBGKnudsenGas.jl` can directly compute the kinematic force $\boldsymbol F / \rho$ and torque $\boldsymbol M / \rho$ about the point that the gas exerts on the immersed body with the function [`force_and_torque_uniform_flow`](@ref). It also computes the Jacobians $\partial (\boldsymbol F / \rho) / \partial \langle \boldsymbol v_0 \rangle$ and $\partial (\boldsymbol M / \rho) / \partial \langle \boldsymbol v_0 \rangle$ with respect ot the mean flow velocity. The Jacobians are computed analytically using the exact gradient of the velocity PDF

```math
\frac{\partial f}{\partial \langle \boldsymbol v_0 \rangle} = \frac{\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle}{\sigma_v^2} f,
```

leading, for example, to the following formula for the Jacobian of force:

```math
\frac{\partial \boldsymbol F}{\partial \langle \boldsymbol v_0 \rangle} = - \frac{\rho  \pi r_\text{dom}^2}{\sigma_v^2} \bigl< \| \boldsymbol v_0 \| \, \Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) \otimes (\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle) \bigr>.
```

## Complete example: Epstein drag

We provide a concrete and complete example of how to use `GBGKnudsenGas.jl` to build your own flow setup using Epstein drag calculation as an example. See the [respective page.](Epstein_drag.md)


## Extending the package

To add another kind of body, you define a new body type (e.g. `Ellipsoid`) and implement the corresponding `collide_with_body` method:

```julia
collide_with_body(particle::Particle, elipsoid::Ellipsoid)::Union{Nothing, Particle}
```

The existing implementation for `Sphere` can be used for reference.

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

Totating bodies, in the limit of small rotational speed, can be added similarly to non-specular collisions. That's because to first order the body rotation merely modifies the collision law.

Currently, arbitrary meshes for bodies are not supported, but that would be a valuable feature to add.