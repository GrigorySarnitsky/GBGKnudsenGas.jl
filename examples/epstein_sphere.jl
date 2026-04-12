using StaticArrays
using GBGKnudsenGas
using LinearAlgebra
using SpecialFunctions
using Plots
using LaTeXStrings
using ForwardDiff

pgfplotsx()

for command in ["\\usepackage{amsmath}"]
  if command ∉ PGFPlotsX.CUSTOM_PREAMBLE
    push!(PGFPlotsX.CUSTOM_PREAMBLE, command)
  end
end

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


function force_MC(body, mean_v, std_v, domain_radius, sample_size)
  T = eltype(mean_v)
  F_sum = SVector{3, T}(0, 0, 0)
  missed_sum = zero(Int)

  for i in 1:sample_size
    p_b, p_a, missed = find_single_collision(body, mean_v, std_v, domain_radius)
    F_sum += norm(p_b.v) * (p_a.v - p_b.v)
    missed_sum += missed
  end

  collision_probability = sample_size / (sample_size + missed_sum)
  force = - F_sum / sample_size * collision_probability * π * domain_radius^2

  return force
end


epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)

const σ_v   = 1.0
const r_dom = 2.0

mean_v_values = [SA[mean_v, 0.0, 0.0] for mean_v in range(0.0, 3.0, 20)]
force_MC_values = force_MC.(Ref(epstein_sphere), mean_v_values, σ_v, r_dom, 10^6)

force_Epstein(mean_v, std_v) = 8/3 * √(2 * π) * std_v * epstein_sphere.r^2 * mean_v

function force_Baines(mean_v, std_v)
  ψ = norm(mean_v) / (√2 * std_v)
  return 3/8 * ψ * ((1/ψ + 1/(2ψ^3)) * exp(-ψ^2) + (1 + 1/ψ^2 - 1/(4ψ^4)) * √π * erf(ψ)) * force_Epstein(mean_v, std_v)
end

force_Kwok(mean_v, std_v) = force_Epstein(mean_v, std_v) * √(1 + 9 * π / 128 * (norm(mean_v) / std_v)^2)


normalize_force(F) = norm(F) / (8/3 * √(2 * π) * σ_v * epstein_sphere.r^2)

plot(
  norm.(mean_v_values) / σ_v,
  [
    normalize_force.(force_Epstein.(mean_v_values, σ_v)),
    normalize_force.(force_Baines.(mean_v_values, σ_v)),
    normalize_force.(force_Kwok.(mean_v_values, σ_v)),
    normalize_force.(force_MC_values)
  ],
  label = ["Epstein (1924)" "Baines et al. (1965)" "Kwok (1975)" "GBGKnudsenGas.jl"], legend = :topleft,
  color = [:black 1 2 :black], linestyle = [:solid :solid :dash :dash],
  ylabel = L"\| \boldsymbol{F} \| \bigm/ (\frac{8}{3} \sqrt{2 \pi} \sigma_v \rho r^2)",
  xlabel = L"\| \langle \boldsymbol v_0 \rangle \| / \sigma_v", xticks = 0:0.5:3
)

display(plot!())

savefig("Epstein_drag.svg")


force_linear_response_MC = ForwardDiff.jacobian(mean_v -> force_MC(epstein_sphere, mean_v, σ_v, r_dom, 10^7), SA[0.0, 0.0, 0.0])

epstein_coefficient = 8/3 * √(2 * π) * σ_v * epstein_sphere.r^2

println("\n(dF / dv) / epstein_coefficient = ")
display(round.(force_linear_response_MC / epstein_coefficient, digits = 5))