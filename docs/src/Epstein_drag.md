# Example: Epstein drag calculation


We provide a concrete example of how to use `GBGKnudsenGas.jl`. Namely, we compute the Epstein drag on a sphere and compare it to analytical predictions.

```@contents
Pages = ["Epstein_drag.md"]
```

## Epstein drag

Consider stationary flow of a Knudsen gas past a sphere with specular collisions. The gas is at thermal equilibrium, with unperturbed (at infinity) temperature $T$, density $\rho$, and mean velocity $\langle \boldsymbol v_0 \rangle$. Angle brackets $\langle {\cdot} \rangle$ denote ensemble averages, and primes ${\cdot}'$ denote the corresponding fluctuations. Before colliding with the body, the initial particle velocity $ \boldsymbol v_0$ follows the Maxwell–Boltzmann probability density function

```math
f(\boldsymbol v_0) = \frac{1}{(\sqrt{2 \pi \sigma_v^2})^3} \exp \biggl[- \frac{(\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle) \cdot (\boldsymbol v_0 - \langle \boldsymbol v_0 \rangle)}{2 \sigma_v^2}\biggr],
```
where

```math
\sigma_v = \sqrt{\frac{k_\text{B} T}{m}},
```
and $m$ is the mass of a single particle. The speed of sound $c$ in the gas is defined through the adiabatic index $\gamma$:

```math
c = \sqrt{\gamma \frac{k_\text{B} T}{m}} = \sqrt \gamma \sigma_v.
```
For a monatomic gas $\gamma = 5/3$. The Mach number $\text{Ma}$ is

```math
\text{Ma} = \frac{\| \langle \boldsymbol v_0 \rangle \|}{c} = \frac{\| \langle \boldsymbol v_0 \rangle \|}{\sigma_v \sqrt \gamma}.
```

Let $\boldsymbol F$ be the force that the flowing gas exerts on the sphere and $r$ be the radius of the sphere. In the limit of small Mach numbers, $\text{Ma} \ll 1$, the Epstein (1924) expression for $F_\text{drag}$ is

```math
\boldsymbol F_\text{Epstein (1924)}
= - \frac{8}{3} \sqrt{2 \pi} \sigma_v \rho r^2 \langle \boldsymbol v_0 \rangle,
```

see the papers

* Epstein, Paul S. (1924). "On the resistance experienced by spheres in their motion through gases." _Physical Review,_ 23(6), 710, [10.1103/PhysRev.23.710](https://doi.org/10.1103/PhysRev.23.710).
* Li, Zhigang, and Wang, Hai (2003). "Drag force, diffusion coefficient, and electric mobility of small particles. I. Theory applicable to the free-molecule regime." _Physical Review E,_ 68(6), 061206, [10.1103/PhysRevE.68.061206](https://doi.org/10.1103/PhysRevE.68.061206)

For arbitrary Mach numbers, a good approximation for $F_\text{drag}$ was given by Kwok (1975):
```math
\boldsymbol F_\text{Kwok (1975)} = \sqrt{1 + \frac{9\pi}{128} \biggr( \frac{\| \langle \boldsymbol v_0 \rangle \|}{\sigma_v} \biggr)^{\!2} } \, \boldsymbol F_\text{Epstein (1924)},
```
see the papers
* Laibe, Guillaume, and Price, Daniel J. (2012). "Dusty gas with smoothed particle hydrodynamics–II. Implicit timestepping and astrophysical drag regimes." _Monthly Notices of the Royal Astronomical Society,_ 420(3), 2365, [10.1111/j.1365-2966.2011.20201.x](https://doi.org/10.1111/j.1365-2966.2011.20201.x)
* Kwok, Sun (1975). "Radiation pressure on grains as a mechanism for mass loss in red giants." _Astrophysical Journal,_ 198(6??), 583, [10.1086/153637](https://doi.org/10.1086/153637).

Kwok's formula is accurate within 1% of the exact expression. The exact formula for $\boldsymbol F$ was given by Baines et al. (1965), but we won't use it, see
* Baines, M. J., Williams, I. P., Asebiomo, A. S., & Agacy, R. L. (1965). "Resistance to the motion of a small sphere moving through a gas." _Monthly Notices of the Royal Astronomical Society,_ 130(1), 63, [10.1093/mnras/130.1.63](https://doi.org/10.1093/mnras/130.1.63).


## Monte Carlo calculation of Epstein drag with `GBGKnudsenGas.jl`

The full code of this computation is available in `examples/epstein_drag.jl`. For a more convinient analysis, we calculate the norm $\| \boldsymbol F \|$ of the force vector $\boldsymbol F$.

### Computing $\| \boldsymbol F \|$

We start by defining our sphere, `epstein_sphere`, with radius $1.0$ and center at the origin $(0, 0, 0)$ of the Cartesian coordinate system:

```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

The calculation of $\| \boldsymbol F \|$ is based on the statistics of the velocity increments of the gas particles. Let $\boldsymbol v_\infty(\boldsymbol v_0, \boldsymbol x_0)$ be the final velocity of a gas particle that initially had velocity $\boldsymbol v_0$ and position $\boldsymbol x_0$. The velocity increment is

```math
\Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) = \boldsymbol v_\infty - \boldsymbol v_0
```
Recall from the [introduction,](index.md) that

```math
\| \boldsymbol F \| = \rho  \pi r_\text{dom}^2 \Bigl\| \bigl< \Delta \boldsymbol v(\boldsymbol v_0, \boldsymbol x_0) \| \boldsymbol v_0 \| \bigr> \Bigr\|.
```

The ensemble average is taken over particles whose initial states ($\boldsymbol v_0$ and $\boldsymbol x_0$) are sampled using `sample_particle_state_uniform_flow`. Let $P(\text{hit})$ denote the probability that a particle collides with the body. Then we can rewrite the expression for $|\boldsymbol F|$ using the average conditioned on a collision:

```math
\| \boldsymbol F \| = \rho  \pi r_\text{dom}^2 P(\text{hit}) \Bigl\| \bigl< \Delta \boldsymbol v(\boldsymbol x_0) \| \boldsymbol v_0 \| \bigm| \text{hit} \bigr> \Bigr\|.
```

To evaluate this conditional average, we first define a helper function `find_single_collision`, which repeatedly draws particles from the initial distribution until one collides with the body. It accepts the following arguments: the body (`body`), the mean velocity $\langle \boldsymbol v_0 \rangle$ (`mean_v`), the velocity standard deviation $\sigma_v$ (`std_v`), and the simulation domain radius $r_\text{dom}$ (`domain_radius`). We later call this function with `epstein_sphere` as the body and $r_\text{dom} = 2$.

The function returns the particle’s initial state (`particle_before`), its state after the collision (`particle_after`), and the number of non-colliding samples encountered before a hit was found (`missed_particles_number`):

```julia
function find_single_collision(body, mean_v, std_v, domain_radius)
  missed_particles_number = 0
  
  while true
    particle_before = sample_particle_state_uniform_flow(mean_v, std_v, domain_radius)
    particle_after  = collide_with_body(particle_before, body)
    
    if !isnothing(particle_after)
      return particle_before, particle_after, missed_particles_number
    end
    
    missed_particles_number += 1
  end 
end
```

To accurately estimate ensemble averages, we need a large sample of particles that collided with the body. The function `run_simulation` below finds `sample_size` particles (later we set `sample_size` to `10^6`) that collided with the `body` and returns two arrays, one with the initial state of the particles before the collision (`particles_before_collision`) and one with the particle states after the collision (`particles_after_collision`). The code is multithreaded to take advantage of a multi-core CPU:

```julia
using Base.Threads

function run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  T = eltype(mean_v)

  particles_before_collision = Vector{Particle{T}}(undef, sample_size)
  particles_after_collision  = Vector{Particle{T}}(undef, sample_size)
  missed_particles_numbers   = Vector{Int}(undef, sample_size)

  @threads for i in 1:sample_size
    particles_before_collision[i], particles_after_collision[i], missed_particles_numbers[i] = find_single_collision(body, mean_v, std_v, domain_radius)
  end

  return particles_before_collision, particles_after_collision, missed_particles_numbers
end
```

Finally, we write a function to calculate $\| \boldsymbol F \|$ for the specified Mach number $\| \langle \boldsymbol v_0 \rangle \|$:

```julia
using Statistics
using LinearAlgebra

function drag_MC(body, mean_v_norm, std_v, domain_radius, sample_size)
  mean_v = SA[mean_v_norm, 0.0, 0.0]
  
  particles_before_collision, particles_after_collision, missed_particles_numbers = run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  
  collision_probability = sample_size / (sum(missed_particles_numbers) + sample_size)
  
  F = -mean(norm(p_b.v) * (p_a.v - p_b.v) for (p_b, p_a) in zip(particles_before_collision, particles_after_collision)) * collision_probability * π * domain_radius^2
  
  return F[1]
end
```

We call it for our `epstein_sphere` and 20 different values $\| \langle \boldsymbol v_0 \rangle \|$. Then we plot the results along with the analytical predictions of Epstein (1924) and Kwok (1975):

```julia
using Plots
using LaTeXStrings

default(;fontfamily = "Computer Modern")

mean_v_norms = range(0.0, 3.0, 20)
drag_values = drag_MC.(Ref(epstein_sphere), mean_v_norms, 1.0, 2.0, 10^6)

drag_Epstein(mean_v_norm, std_v) = 8/3 * √(2 * π) * std_v * epstein_sphere.r^2 * mean_v_norm
drag_Kwok(mean_v_norm, std_v) = drag_Epstein(mean_v_norm, std_v) * √(1 + 9 * π / 128 * (mean_v_norm / std_v)^2)


plot(
  mean_v_norms, [drag_Epstein.(mean_v_norms, 1.0), drag_Kwok.(mean_v_norms, 1.0), drag_values] ./ (8/3 * √(2 * π) * 1.0 * epstein_sphere.r^2),
  label = ["Epstein (1924)" "Kwok (1975)" "GBGKnudsenGas.jl"],
  color = [:black 1 1], linestyle = [:solid :solid :dash],
  ylabel = L"‖\! 𝑭\ ‖ / (\frac{8}{3} \sqrt{2 \pi} \sigma_v \rho r^2)",
  xlabel = L"‖\! 𝒗\ ‖ / σ_v", xticks = 0:0.5:3
)

display(plot!())

savefig("Epstein_drag.svg")
```

![Epstein drag plot](Epstein_drag.svg)

### Computing Epstein coefficient

```julia
using ForwardDiff

epstein_coefficient_MC = ForwardDiff.derivative(mean_v_norm -> drag_MC(epstein_sphere, mean_v_norm, 1.0, 2.0, 10^6), mean_v_norm)

epstein_coefficient = 8/3 * √(2 * π) *  epstein_sphere.r^2

println(epstein_coefficient_MC / epstein_coefficient)
```