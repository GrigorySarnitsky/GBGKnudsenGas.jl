# Example: Epstein drag calculation


We provide a concrete example of how to use `GBGKnudsenGas.jl`. Namely, we compute the Epstein drag on a sphere and compare it to analytical predictions.

```@contents
Pages = ["Epstein_drag.md"]
```

## Epstein drag

Consider a flow of a Knudsen gas past a sphere with specular collisions. The gas is at thermal equilibrium with unperturbed (at infinity) temperature $T$, density $\rho$ and mean velocity $\langle \boldsymbol v \rangle_\infty$. Its speed of sound $c$ is defined through its adiabatic index $\gamma$ and mass $m$ of its individual particle:

```math
c = \sqrt{\gamma \frac{k_\text{B} T}{m}}.
```

For monatomic gas $\gamma = 5/3$. In addition, we can define the standard deviation of an individual velocity component $\sigma_v$,

```math
  \sigma_v = \frac{1}{3} \langle \boldsymbol v' \cdot \boldsymbol v' \rangle_\infty = \sqrt{\gamma \frac{k_\text{B} T}{m}}.
```

The respective Mach number $\text{Ma}$ is

```math
\text{Ma} = \frac{\| \langle \boldsymbol v \rangle_\infty \|}{c}.
```

The gas flow exerts a drag force on the sphere of magnitude $F_\text{drag}$ and directed opposite direction to $\langle \boldsymbol v \rangle_\infty$. F

Our goal is to calculate the drag excreted on the sphere by the flow.


per-collision change in particle velocity, $\langle \Delta \boldsymbol{v} \rangle_\text{c}$. Ensemble averages are denoted by angle brackets, $\langle {\cdot} \rangle$; the corresponding fluctuating parts are denoted with primes ${\cdot}'$. The index $\text{c}$ indicates averaging over the particles that actually collide with the sphere. In standard notation for conditional averages,

```math
\langle \Delta \boldsymbol{v} \rangle_\text{c} = \langle \boldsymbol{v}_\text{after} - \boldsymbol{v}_\text{before} | \text{particle has collided} \rangle,
```

where $`\boldsymbol{v}_\text{before}`$ and $`\boldsymbol{v}_\text{after}`$ are the particle's velocity before and after the collision with the sphere.

We want to check Epstein's 1924 calculation that $`\langle \Delta \boldsymbol{v} \rangle_\text{c}`$ is proportional to the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ at infinity:

```math
  \langle \Delta \boldsymbol{v} \rangle_\text{c} = - \delta \langle \boldsymbol{v} \rangle_\infty,
```

with $\delta = 1$ for specular collisions. In this form this relation can be obtained who knows how. See also the original paper: Epstein, Paul S. "On the resistance experienced by spheres in their motion through gases," _Physical Review,_ 23.6 (1924), 710.

## Code for Monte Carlo calcualtion of the Epstein drag

Full code of this computation is available in `examples/epstein_sphere.jl`.

We start by defining our sphere, which we call `epstein_sphere`. We set its radius to $1.0$ and center it in the origin $(0, 0, 0)$ of our Cartesian coordinate system:

```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Next, we need to generate the initial state of the gas particles before they hit the sphere. The gas is at thermal equilibrium, moving past the sphere with the mean flow velocity $`\langle \boldsymbol{v} \rangle_\infty`$ (denoted `mean_v` in code). Instead of temperature, we specify the standard deviation $`\sigma_{v, \infty}`$,

```math
 \sigma_{v, \infty} = \sqrt{\frac{1}{3} \langle \boldsymbol v' \cdot \boldsymbol v' \rangle_\infty},
```
 
 of particle velocities far away from the sphere (denoted `std_v` in code). Assuming a monotomic gas, the speed of sound is $c = \sqrt{5/3} \sigma_{v, \infty}$, and the Mach number $\text{Ma}$ of the flow is

```math
\text{Ma} = \sqrt{\frac{3}{5}} \frac{\| \langle \boldsymbol{v} \rangle_\infty \|}{\sigma_{v, \infty}}.
```

We initially place all particles on the boundary of a spherical domain of radius `domain_radius`, which should be sufficiently larger than the unit radius of `epstein_sphere`. Later we set `domain_radius` to `20`.

Thus defined initial state of a gas particle can be obtained with `get_initial_particle_state(mean_v, std_v, domain_radius)`. However, we are interested only in particles that will collide with the sphere. Given that `domain_radius` is much larger than the radius of `epstein_sphere`, most of the generated particles will miss the `epstein_sphere`. Hence we write a simple function `find_single_collision` that keeps drawing particles from the initial distribution until it finds one that collides with the `body`. Later we call `find_single_collision` with `epstein_sphere` for the `body`. The function returns the initial state of the particle `particle_before`, as well as its state after the collision `particle_after`:


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

To accurately estimate ensemble averages, we need a large sample of particles that collided with the body. The function `run_simulation` below finds `sample_size` particles (later we set `sample_size` to `10^6`) that collided with the `body` and returns two arrays, one with the initial state of the particles before the collision (`particles_before_collision`) and one with the particle states after the collision (`particles_after_collision`). The code is multithreaded to take advantage of a multi-core CPU:

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

Finally, we write a function to calculate $`\langle \Delta \boldsymbol{v} \rangle_\text{c}`$ for the specified Mach number `Ma`:

```julia
using Statistics

function mean_Δv_per_collision(body, Ma)
  std_v = 1.0
  mean_v = SA[Ma * std_v * √(5/3), 0.0, 0.0]
  domain_radius = 20.0
  sample_size = 10^6
  particles_before_collision, particles_after_collision = run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  return mean(p_a.v - p_b.v for (p_b, p_a) in zip(particles_before_collision, particles_after_collision))
end
```

We call it for our `epstein_sphere` and 10 different values of the Mach number, and calculate Epstein’s 
$\delta$:

```julia
Ma_s = range(0.0, 3.0, 10)
mean_Δv_s = @time mean_Δv_per_collision.(Ref(epstein_sphere), Ma_s)

for (Δv, Ma) in zip(mean_Δv_s, Ma_s)
  println("Ma = ", round(Ma, digits = 2), ",\t δ = ", -round(Δv[1] / (√(5/3) * Ma), digits = 5))
end
```

This should print something like

```
105.765674 seconds (758.63 k allocations: 1.079 GiB, 0.17% gc time, 1.02% compilation time)
Ma = 0.0,	 δ = Inf
Ma = 0.33,	 δ = 0.9967
Ma = 0.67,	 δ = 1.00221
Ma = 1.0,	 δ = 0.99963
Ma = 1.33,	 δ = 0.99864
Ma = 1.67,	 δ = 1.00079
Ma = 2.0,	 δ = 0.99864
Ma = 2.33,	 δ = 0.99994
Ma = 2.67,	 δ = 1.00015
Ma = 3.0,	 δ = 1.00037
```

As expected, $\delta \approx 1$ and is independent of the Mach number.