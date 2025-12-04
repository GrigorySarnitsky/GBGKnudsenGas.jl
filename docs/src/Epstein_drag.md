# Example: Epstein drag calculation


We provide a concrete example of how to use `GBGKnudsenGas.jl`. Namely, we compute the Epstein drag on a sphere and compare it to analytical predictions.

```@contents
Pages = ["Epstein_drag.md"]
```

## Epstein drag

Consider a flow of a Knudsen gas past a sphere with specular collisions. The gas is at thermal equilibrium with unperturbed (at infinity) temperature $T$, density $\rho$ and mean velocity $\langle \boldsymbol v \rangle_\infty$. Angle brackets $\langle {\cdot} \rangle$ denote ensemble averages, and primes ${\cdot}'$ denote the respective fluctuating parts. The speed of sound $c$ in the gas is defined through the adiabatic index $\gamma$ and mass $m$ of an individual gas particle:
```math
c = \sqrt{\gamma \frac{k_\text{B} T}{m}}.
```
For a monatomic gas $\gamma = 5/3$. In addition, we can define the standard deviation of an individual velocity component $\sigma_v$.
```math
  \sigma_v = \frac{1}{3} \langle \boldsymbol v' \cdot \boldsymbol v' \rangle_\infty = \sqrt{\frac{k_\text{B} T}{m}} = \frac{c}{\sqrt{\gamma}}.
```
The Mach number $\text{Ma}$ is
```math
\text{Ma} = \frac{\| \langle \boldsymbol v \rangle_\infty \|}{c} = \frac{\| \langle \boldsymbol v \rangle_\infty \|}{\sigma_v \sqrt \gamma}.
```
The gas flow exerts a drag force on the sphere of magnitude $F_\text{drag}$ directed opposite to $\langle \boldsymbol v \rangle_\infty$. Let the sphere radius be $r$. In the limit of small Mach numbers, $\text{Ma} \ll 1$, the Epstein (1924) expression for $F_\text{drag}$ is
```math
F_\text{drag, Epstein (1924)}
= \frac{8}{3} \sqrt{2 \pi} \sigma_v \rho r^2 \| \langle \boldsymbol v \rangle_\infty \|
= \frac{8}{3} \sqrt{2 \pi} \sigma_v^2 \rho r^2 \sqrt{\gamma} \text{Ma},
```
see the papers

* Epstein, Paul S. (1924). "On the resistance experienced by spheres in their motion through gases." _Physical Review,_ 23(6), 710, [10.1103/PhysRev.23.710](https://doi.org/10.1103/PhysRev.23.710).
* Li, Zhigang, and Wang, Hai (2003). "Drag force, diffusion coefficient, and electric mobility of small particles. I. Theory applicable to the free-molecule regime." _Physical Review E,_ 68(6), 061206, [10.1103/PhysRevE.68.061206](https://doi.org/10.1103/PhysRevE.68.061206)

For arbitrary Mach numbers, a good approximation for $F_\text{drag}$ was given by Kwok (1975):
```math
F_\text{drag, Kwok (1975)} = \sqrt{1 + \frac{9\pi}{128} \gamma \text{Ma}^2} \, F_\text{drag, Epstein (1924)},
```
see the papers
* Laibe, Guillaume, and Price, Daniel J. (2012). "Dusty gas with smoothed particle hydrodynamics–II. Implicit timestepping and astrophysical drag regimes." _Monthly Notices of the Royal Astronomical Society,_ 420(3), 2365, [10.1111/j.1365-2966.2011.20201.x](https://doi.org/10.1111/j.1365-2966.2011.20201.x)
* Kwok, Sun (1975). "Radiation pressure on grains as a mechanism for mass loss in red giants." _Astrophysical Journal,_ 198(6??), 583, [10.1086/153637](https://doi.org/10.1086/153637).

Kwok's formula is accurate within 1% of the exact expression. The exact formula for $F_\text{drag}$ by Baines et al. (1965), but we won't use it, see
* Baines, M. J., Williams, I. P., Asebiomo, A. S., & Agacy, R. L. (1965). "Resistance to the motion of a small sphere moving through a gas." _Monthly Notices of the Royal Astronomical Society,_ 130(1), 63, [10.1093/mnras/130.1.63](https://doi.org/10.1093/mnras/130.1.63).


## Monte Carlo calculations of the Epstein drag with `GBGKnudsenGas.jl`

Full code of this computation is available in `examples/epstein_drag.jl`.

```math
\int
```

We start by defining our sphere, which we call `epstein_sphere`. We set its radius to $1.0$ and center it in the origin $(0, 0, 0)$ of our Cartesian coordinate system:

```julia
using StaticArrays
using GBGKnudsenGas

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)
```

Next, we need to generate the initial state of the gas particles before they hit the sphere. We need to sample the ...

This initial state of a gas particle can be obtained with `get_initial_particle_state(mean_v, std_v, domain_radius)`. However, we are interested only in particles that will collide with the sphere. Given that `domain_radius` is much larger than the radius of `epstein_sphere`, most of the generated particles will miss the `epstein_sphere`. Hence we write a simple function `find_single_collision` that keeps drawing particles from the initial distribution until it finds one that collides with the `body`. Later we call `find_single_collision` with `epstein_sphere` for the `body`. The function returns the initial state of the particle `particle_before`, as well as its state after the collision `particle_after`:


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

![Epstein drag plot](Epstein_drag.svg)