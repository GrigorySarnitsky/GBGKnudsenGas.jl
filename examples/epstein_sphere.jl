using StaticArrays
using Statistics
using GBGKnudsenGas
using Base.Threads

function find_single_collision(body, mean_v, std_v, domain_radius)
  while true
    particle_before = get_initial_particle_state(mean_v, std_v, domain_radius)
    particle_after  = collide_with_body(particle_before, body)
    if !isnothing(particle_after)
      return particle_before, particle_after
    end
  end 
end

function run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  particles_before_collision = Array{Particle, 1}(undef, sample_size)
  particles_after_collision  = Array{Particle, 1}(undef, sample_size)
  @threads for i in 1:sample_size
    particles_before_collision[i], particles_after_collision[i] = find_single_collision(body, mean_v, std_v, domain_radius)
  end
  return particles_before_collision, particles_after_collision
end

function mean_Δv_per_collision(body, Ma)
  std_v = 1.0
  mean_v = SA[Ma * std_v * √(5/3), 0.0, 0.0]
  domain_radius = 20.0
  sample_size = 10^6
  particles_before_collision, particles_after_collision = run_simulation(body, mean_v, std_v, domain_radius, sample_size)
  return mean(p_a.v - p_b.v for (p_b, p_a) in zip(particles_before_collision, particles_after_collision))
end

epstein_sphere = Sphere(SA[0.0, 0.0, 0.0], 1.0)

Ma_s = range(0.0, 3.0, 10)
mean_Δv_s = @time mean_Δv_per_collision.(Ref(epstein_sphere), Ma_s)

for (Δv, Ma) in zip(mean_Δv_s, Ma_s)
  println("Ma = ", round(Ma, digits = 2), ",\t δ = ", -round(Δv[1] / (√(5/3) * Ma), digits = 5))
end