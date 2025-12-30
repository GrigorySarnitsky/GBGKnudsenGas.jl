# Development tips for GBGKnudsenGas.jl

Here we assmeble tips on the advanced usage, extension and development of  `GBGKnudsenGas.jl`. It is written a bit loosly, but if you got to the point you need to read it, you should be abgle to get through.

## Estimators for means and their standard errors

When studying the flow of a Knudsen gas past a body, it is natural to collect statistics only from gas particles that actually hit the body, since particles that miss it neither change their state nor affect the body. In a Monte Carlo simulation, this corresponds to collecting data based on a fixed, predetermined number $N_\text{hit}$ of collided particles (and stopping once we reach $N_\text{hit}$). If $N_\text{miss}$ denotes the number of particles that miss the body, then the total number of simulated particles is $N_\text{miss} + N_\text{hit}$. Unlike $N_\text{hit}$, $N_\text{miss}$ is a random variable, and consequently so is the total number of simulated particles. This randomness requires care when computing statistics across all of the particles.

Let $\boldsymbol F$ be a (tensor) quantity to be averaged over the particles, for example, the change in particle momentum. Crucially, for non-colliding particles, $\boldsymbol F = 0$. Our goal is to compute the average $\langle \boldsymbol F \rangle$ across all of the particles, both that hit the body and that missed it. A naive Monte Carlo estimator one can come up with is

```math
  \langle \boldsymbol F \rangle_\text{MC, naive} = \frac{1}{N_\text{miss} + N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i,
```

with $\boldsymbol F_i$ being the value of $\boldsymbol F$ for $i$-th colliding particle. The respective naive covariance of $\langle \boldsymbol F \rangle_\text{MC, naive}$ (its measure of accuracy) would be

```math
  \begin{align*}
    \operatorname{cov}(\langle \boldsymbol F \rangle_\text{MC, naive})
    &= \bigl< \bigl[ \langle \boldsymbol F \rangle_\text{MC, naive} - \langle \boldsymbol F \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \rangle_\text{MC, naive} - \langle \boldsymbol F \rangle\bigr] \bigr> \\
    &= \frac{1}{(N_\text{miss} + N_\text{hit}) (N_\text{miss} + N_\text{hit} - 1)} \biggl[ \sum_{i = 1}^{N_\text{hit}} \bigl[ \boldsymbol F_i - \langle \boldsymbol F \rangle_\text{MC, naive} \bigr] \otimes \bigl[ \boldsymbol F_i - \langle \boldsymbol F \rangle_\text{MC, naive} \bigr] + N_\text{miss} \langle \boldsymbol F \rangle_\text{MC, naive} \otimes \langle \boldsymbol F \rangle_\text{MC, naive} \biggr].
  \end{align*}
```

These would be correct and unbiased estimators if both $N_\text{hit}$ and $N_\text{miss}$ were deterministic quantities. However, since $N_\text{miss}$  is random, we instead take a different approach.

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
P_\text{MC}(\text{hit}) = \frac{N_\text{hit} - 1}{N_\text{miss} + N_\text{hit} - 1}.
```

Then the estimator for $\langle \boldsymbol F \rangle$ can be written down as

```math
\langle \boldsymbol F \rangle_\text{MC} = \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} P_\text{MC}(\text{hit}) = \frac{1}{N_\text{miss} + N_\text{hit} - 1} \frac{N_\text{hit} - 1}{N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i.
```
Compare this to the naive estimator $\langle \boldsymbol F \rangle_\text{MC, naive}$.

Next, we want to quantify the accuracy of $\langle \boldsymbol F \rangle_\text{MC}$, i.e. find its covariance. For this we first find the covariances of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} $ and $P_\text{MC}(\text{hit})$.

### Conditional estimator covariances

#### Covariance of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$

To quantify the accuracy of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$, we need to calculate its covariance:

```math
  \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} )
  = \Bigl< \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \Bigr>.
```

It helps to remember that $\langle {\cdot} \rangle$ here denotes averaging over the ensemble of Monte Carlo simulations, each of each having exactly $N_\text{hit}$ collided particles.

To calculate this covariance, we use the unbiased property of $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$ (its average coincides with the true average):

```math
\bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigm| \text{hit} \bigr> = \langle \boldsymbol F \mid \text{hit} \rangle,
```

and that two measurements, $i$ and $j$, are statistically independent, leading to

```math
\langle \boldsymbol F_i \otimes \boldsymbol F_j \mid \text{hit} \rangle = \langle \boldsymbol F_i  \mid \text{hit} \rangle \otimes \langle \boldsymbol F_j \mid \text{hit} \rangle
\quad\text{if $i \neq j$}.
```

With that, we proceed:
```math
  \begin{align*}
    \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} )
    &= \Bigl< \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \Bigr> \\
    &= \bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \Bigl< \frac{1}{N_\text{hit}^2} \sum_{i = 1}^{N_\text{hit}} \sum_{j = 1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_j \Bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \frac{1}{N_\text{hit}^2} \sum_{i = 1}^{N_\text{hit}} \sum_{j = 1}^{N_\text{hit}} \langle \boldsymbol F_i \otimes \boldsymbol F_j \mid \text{hit} \rangle - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \frac{1}{N_\text{hit}^2} \Bigl[ N_\text{hit} \langle \boldsymbol F \otimes \boldsymbol F \mid \text{hit} \rangle + N_\text{hit} (N_\text{hit} - 1) \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle  \boldsymbol F \mid \text{hit} \rangle \Bigr] - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \frac{1}{N_\text{hit}} \bigl[\langle \boldsymbol F \otimes \boldsymbol F \mid \text{hit} \rangle - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \bigr].
  \end{align*}
```
Thus we confirmed the usual formula for the [standard error:](https://en.wikipedia.org/wiki/Standard_error)
```math
  \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC})
  = \frac{1}{N_\text{hit}} \operatorname{cov}(\boldsymbol F \mid \text{hit}).
```

Using unbiased estimator for the covariance

```math
  \operatorname{cov}_\text{MC}(\boldsymbol F \mid \text{hit}) = \frac{1}{N_\text{hit} - 1} \sum_{i=1}^{N_\text{hit}} \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr] \otimes \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr]
```

we finally get
```math
  \begin{align*}
    \operatorname{cov}_\text{MC}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} )
    &= \frac{1}{N_\text{hit} (N_\text{hit} - 1)} \sum_{i=1}^{N_\text{hit}} \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr] \otimes \bigl[ \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigr] \\
    &= \frac{1}{N_\text{hit} (N_\text{hit} - 1)} \biggl( \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_i - N_\text{hit} \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \biggr) \\
    &= \frac{1}{N_\text{hit} - 1} \biggl[ \frac{1}{N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}\biggr].
  \end{align*}
```

#### Variance of $P_\text{MC}(\text{hit})$

Now we calculate the variance of $P_\text{MC}(\text{hit})$,

```math
  \operatorname{var}(P_\text{MC}(\text{hit})) = \bigl< \bigl[P_\text{MC}(\text{hit}) - P(\text{hit})\bigr]^2 \bigr>.
```

For this we need to use the probability mass function $f$ of $N_\text{miss}$:

```math
  f(N_\text{miss}) = \frac{(N_\text{miss} + N_\text{hit} - 1)!}{(N_\text{hit} - 1)! N_\text{miss}!} P^{N_\text{hit}}(\text{hit}) \bigl[1 - P(\text{hit}) \bigr]^{N_\text{miss}}.
```

Now we proceed with calculating $P_\text{MC}(\text{hit})$, studying in the limit of large $N_\text{hit}$:
```math
  \begin{align*}
    \operatorname{var}(P_\text{MC}(\text{hit}))
    &= \bigl< \bigl[P_\text{MC}(\text{hit}) - P(\text{hit})\bigr]^2 \bigr> \\
    &= \sum_{N_{\text{miss}} = 0}^\infty \biggl( \frac{N_\text{hit} - 1}{N_\text{miss} + N_\text{hit} - 1} - P(\text{hit})\biggr)^{\!2} \frac{(N_\text{miss} + N_\text{hit} - 1)!}{(N_\text{hit} - 1)! N_\text{miss}!} P^{N_\text{hit}}(\text{hit}) \bigl[1 - P(\text{hit}) \bigr]^{N_\text{miss}} \\
    &= P^2(\text{hit}) \bigl[ {}_2F_1(1, 1, N_\text{hit}, 1 - P(\text{hit})) - 1 \bigr] \\
    &= P^2(\text{hit}) \biggl[ \frac{1 - P(\text{hit})}{N_\text{hit}} + \frac{2(1 - P(\text{hit}))^2}{N_\text{hit}(N_\text{hit} + 1)} + \frac{6 (1 - P(\text{hit}))^3}{N_\text{hit}(N_\text{hit} + 1)(N_\text{hit} + 2)} + \ldots \biggr],
  \end{align*}
```

where is ${}_2 F_1$ the [hypergeometric function](https://en.wikipedia.org/wiki/Hypergeometric_function). The above can be checked with the following Mathematica code:

```mathematica
f = (m + n - 1)!/((n - 1)! m!) (1 - p)^m p^n;
pmc = (n - 1)/(m + n - 1);
s = FullSimplify[Sum[(pmc - p)^2 f, {m, 0, Infinity}]]
Sum[SeriesCoefficient[s, {(1 - p), 0, k}] (1 - p)^k, {k, 0, 4}]
```

with an output

```mathematica
p^2 (-1 + Hypergeometric2F1[1, 1, n, 1 - p])

((1 - p) p^2)/n + (2 (1 - p)^2 p^2)/(n (1 + n)) + (6 (1 - p)^3 p^2)/(n (1 + n) (2 + n)) + (24 (1 - p)^4 p^2)/(n (1 + n) (2 + n) (3 + n))
```

Thus to estimate the accuracy of $P_\text{MC}(\text{hit})$ we use the below estimator for its variance:

```math
  \operatorname{var}_\text{MC} (P_\text{MC}(\text{hit})) = \frac{1}{N_\text{hit}} P_\text{MC}^2(\text{hit}) \bigl[1 - P_\text{MC}(\text{hit})\bigr].
```

### Covariance of the unconditional estimator $\langle \boldsymbol F \rangle$

With the known covariances of $P_\text{MC}$ and $\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}$ we are well equiped for finding the covariacne of $\langle \boldsymbol F \rangle_\text{MC} = \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} P_\text{MC}(\text{hit})$. Since we a dealing with the product of two independent random variables, its coavariacne is

```math
  \begin{align*}
  \operatorname{cov}(\langle \boldsymbol F \rangle_\text{MC})
  &= \Bigl< \bigl[ \langle \boldsymbol F \rangle_\text{MC} - \langle \boldsymbol F \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \rangle_\text{MC} - \langle \boldsymbol F \rangle \bigr] \Bigr> \\
  &= \operatorname{cov}(\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} P_\text{MC}(\text{hit})) \\
  &= \Bigl< P_\text{MC}(\text{hit})^2 \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \Bigr> - \langle \boldsymbol F \rangle \otimes \langle \boldsymbol F \rangle \\
  &= \bigl< P_\text{MC}(\text{hit})^2 \bigr> \Bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \Bigr> - \langle \boldsymbol F \rangle \otimes \langle \boldsymbol F \rangle \\
  &= \Bigl[ \operatorname{var} (P_\text{MC}(\text{hit})) + P(\text{hit})^2 \Bigr] \Bigl[\operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}) + \langle \boldsymbol F \mid \text{hit} \rangle\otimes \langle \boldsymbol F \mid \text{hit} \rangle \Bigr] - \langle \boldsymbol F \rangle \otimes \langle \boldsymbol F \rangle \\
  &= \bigl[ P(\text{hit})^2 + \operatorname{var} (P_\text{MC}(\text{hit})) \bigr] \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}) + \operatorname{var} (P_\text{MC}(\text{hit})) \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle.
  \end{align*}
```

So the Monte Carlo estimator for $\operatorname{cov}(\langle \boldsymbol F \rangle_\text{MC})$ is

```math
  \operatorname{cov}_\text{MC}(\langle \boldsymbol F \rangle_\text{MC})
  = \bigl[ P_\text{MC}(\text{hit})^2 + \operatorname{var} (P_\text{MC}(\text{hit})) \bigr] \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}) + \operatorname{var} (P_\text{MC}(\text{hit})) \, \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}.
```

### Summary of the estimators

To summarize, consider a Monte Carlo simulation of a flow of a Knudsen gas past a body. We run it untill we find $N_\text{hit}$ particles that collided with the body and then stop. During the simulatin we also count $N_\text{miss}$ particles that hit the body. Let $\boldsymbol F$ be a (tensor) quantity to be averaged over the particles, for example, the change in particle momentum. Crucially, for non-colliding particles, $\boldsymbol F = 0$. Our goal is to compute the average $\langle \boldsymbol F \rangle$ across all of the particles, both that hit the body and that missed it. Porper Monte carlo estimators for $\langle \boldsymbol F \rangle$ and its covariance are:

```math
\boxed{
  \begin{align*}
    \langle \boldsymbol F \rangle_\text{MC}
    &=\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} P_\text{MC}(\text{hit}) \\[2ex]
    \operatorname{cov}(\langle \boldsymbol F \rangle_\text{MC})
    &= \Bigl< \bigl[ \langle \boldsymbol F \rangle_\text{MC} - \langle \boldsymbol F \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \rangle_\text{MC} - \langle \boldsymbol F \rangle \bigr] \Bigr> \\
    &= P_\text{MC}^2(\text{hit}) \, \operatorname{cov}_\text{MC}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit}) + \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \, \operatorname{var}_\text{MC}(P_\text{MC}(\text{hit})) \\[1ex]
    \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}
    &= \frac{1}{N_\text{hit}} \sum_{i = 1}^{N_\text{hit}} \boldsymbol F_i, \\[2ex]
    P_\text{MC}(\text{hit})
    &= \frac{N_\text{hit} - 1}{N_\text{miss} + N_\text{hit} - 1}, \\[2ex]
    \operatorname{cov}_\text{MC}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
    &= \frac{1}{N_\text{hit} - 1} \biggl[ \frac{1}{N_\text{hit}} \sum_{i=1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_i - \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC}\biggr], \\[2ex]
    \operatorname{var}_\text{MC} (P_\text{MC}(\text{hit}))
    &= \frac{1}{N_\text{hit}} P_\text{MC}^2(\text{hit}) \bigl[1 - P_\text{MC}(\text{hit})\bigr].
  \end{align*}
}
```

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