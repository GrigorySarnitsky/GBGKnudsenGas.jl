# Development tips for GBGKnudsenGas.jl

## Estimators for means and their standard errors

When studying the flow of a Knudsen gas past a body, it is natural to collect statistics only from gas particles that actually hit the body, since particles that miss it neither change their state nor affect the body. In a Monte Carlo simulation, this corresponds to collecting data based on a fixed, predetermined number $N_\text{hit}$ of collided particles (and stopping once we reach $N_\text{hit}$). If $N_\text{miss}$ denotes the number of particles that miss the body, then the total number of simulated particles is $N_\text{hit} + N_\text{miss}$. Unlike $N_\text{hit}$, $N_\text{miss}$ is a random variable, and consequently so is the total number of simulated particles. This randomness requires care when computing statistics across all of the particles.

Let $\boldsymbol F$ be a (tensor) quantity to be averaged over the particles, for example, the change in particle momentum. Crucially, for non-colliding particles, $\boldsymbol F = 0$. Our goal is to compute the average $\langle \boldsymbol F \rangle$ across all of the particles, both that hit the body and that missed it. A naive Monte Carlo estimator one can come up with is

```math
  \langle \boldsymbol F \rangle_\text{MC, naive} = \frac{1}{N_\text{hit} + N_\text{miss}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i,
```

with $\boldsymbol F_i$ being the value of $\boldsymbol F$ for $i$-th colliding particle. This would be a correct and unbiased estimator for $\langle \boldsymbol F \rangle$ if both $N_\text{hit}$ and $N_\text{miss}$ were deterministic quantities. However, since $N_\text{miss}$  is random, we instead take a different approach.

Let $P(\text{hit})$ denote the probability that a particle hits the body, and let $\langle \boldsymbol F \mid \text{hit} \rangle$ denote the conditional average of $\boldsymbol F$ over particles that hit the body. Then $\langle \boldsymbol F \rangle$ can be written as

```math
\langle \boldsymbol F \rangle = \langle \boldsymbol F \mid \text{hit} \rangle P(\text{hit}).
```

Since $N_\text{hit}$ is deterministic, the unbiased Monte Carlo estimator for $\langle \boldsymbol F \mid \text{hit} \rangle$ is

```math
\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} = \frac{1}{N_\text{hit}} \sum_{i = 1}^{N_\text{hit}} \boldsymbol F_i.
```

In our setup, we continue the simulation until the predetermined number of collisions $N_\text{hit}$ is reached, at which point we immediately stop. This implies that $N_\text{miss}$ follows the [negative binomial distribution](https://en.wikipedia.org/wiki/Negative_binomial_distribution). The corresponding unbiased Monte Carlo estimator for $P(\text{hit})$ is

```math
P_\text{MC}(\text{hit}) = \frac{N_\text{hit} - 1}{N_\text{hit} + N_\text{miss} - 1}.
```

Then the estimator for $\langle \boldsymbol F \rangle$ can be written down as

```math
\langle \boldsymbol F \rangle_\text{MC} = \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} P_\text{MC}(\text{hit}) = \frac{1}{N_\text{hit} + N_\text{miss} - 1} \frac{N_\text{hit} - 1}{N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i.
```
Compare this to the naive estimator $\langle \boldsymbol F \rangle_\text{MC, naive}$. Next, we want to quantify the accuracy of $\langle \boldsymbol F \rangle_\text{MC}$, i.e. find its covariance. For this we first find the covariances of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} $ and $P_\text{MC}(\text{hit})$.

### Conditional estimator covariances

#### Covariance of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$

To quantify the accuracy of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$, we need to calculate its covariance:

```math
  \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
  = \Bigl< \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \Bigm| \text{hit} \Bigr>.
```

To calculate this covariance, we use the unbiased property of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$ (its average coincides with the true average):

```math
\bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigm| \text{hit} \bigr> = \langle \boldsymbol F \mid \text{hit} \rangle,
```

and that two measurements, $i$ and $j$, are statistically independent, leading to

```math
\langle \boldsymbol F_i \otimes \boldsymbol F_j \mid \text{hit} \rangle = \langle \boldsymbol F_i  \mid \text{hit} \rangle \otimes \langle \boldsymbol F_j \mid \text{hit} \rangle
\quad\text{if $i \neq j$}.
```


Next,
```math
  \begin{align*}
    \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
    &= \Bigl< \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \Bigm| \text{hit} \Bigr> \\
    &= \bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigm| \text{hit} \bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \Bigl< \frac{1}{N_\text{hit}^2} \sum_{i = 1}^{N_\text{hit}} \sum_{j = 1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_j \Bigm| \text{hit} \Bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \frac{1}{N_\text{hit}^2} \Bigl[ N_\text{hit} \langle \boldsymbol F \otimes \boldsymbol F \mid \text{hit} \rangle + N_\text{hit} (N_\text{hit} - 1) \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle  \boldsymbol F \mid \text{hit} \rangle \Bigr] - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \frac{1}{N_\text{hit}} \bigl[\langle \boldsymbol F \otimes \boldsymbol F \mid \text{hit} \rangle - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \bigr].
  \end{align*}
```
Thus we showed
```math
  \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
  = \frac{1}{N_\text{hit}} \operatorname{cov}(\boldsymbol F \mid \text{hit}).
```

Using unbiased estimator for the covariance
```math
  \begin{align*}
    \operatorname{cov}_\text{MC}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
    &= \frac{1}{N_\text{hit} (N_\text{hit} - 1)} \sum_{i=1}^{N_\text{hit}} \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr] \otimes \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr] \\
    &= \frac{1}{N_\text{hit} (N_\text{hit} - 1)} \biggl( \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_i - N_\text{hit} \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \biggr) \\
    &= \frac{1}{N_\text{hit} - 1} \biggl[ \frac{1}{N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}\biggr].
  \end{align*}
```

#### Covariance of $P_\text{MC}(\text{hit})$


The accuracy. Using the unbiaseness of $P_\text{MC}(\text{hit})$, we calculate
```math
  \begin{align*}
    \operatorname{var}(P_\text{MC}(\text{hit}))
    &= \bigl< \bigl[P_\text{MC}(\text{hit}) - P(\text{hit})\bigr]^2 \bigr> \\
    &= \langle P_\text{MC}^2 \rangle - P^2(\text{hit}) \\
    &= \frac{N_\text{hit} - 1}{N_\text{hit} + N_\text{miss} - 1},
  \end{align*}
```



### Covariance of the unconditional estimator $\langle \boldsymbol F \rangle$

## Extending the package

### New body kinds

To add another kind of body, you define a new body type (e.g. `Ellipsoid`) and implement the corresponding `collide_with_body` method:

```julia
collide_with_body(particle::Particle, elipsoid::Ellipsoid)::Union{Nothing, Particle}
```

The existing implementation for `Sphere` can be used for reference.

### Non-specular collisions

As for non-specular collisions, I envision the following architecture. Each body type is extended with the property `collision_type`, which tells how the particles interact with the body. Then the function `collide_with_body` should have a branching that chooses the correct collision law based on the body's `collision_type`.

### Rotating bodies

Rotating bodies, in the limit of small rotational speed, can be added similarly to non-specular collisions. That's because to first order the body rotation merely modifies the collision law.

### Arbitrary meshes

Currently, arbitrary meshes for bodies are not supported, but that would be a valuable feature to add.