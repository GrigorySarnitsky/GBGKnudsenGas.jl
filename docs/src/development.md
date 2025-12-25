# Development tips for GBGKnudsenGas.jl

## Estimators for means and their standard errors

```math
\langle \boldsymbol F \rangle = P(\text{hit}) \langle \boldsymbol F \mid \text{hit} \rangle.
```

Unbiased estimators $\langle \boldsymbol F \mid \text{hit} \rangle$ for

```math
P_\text{MC}(\text{hit}) = \frac{N_\text{hit} - 1}{N_\text{hit} + N_\text{miss} - 1}, \qquad
\langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} = \frac{1}{N_\text{hit}} \sum_{i = 1}^{N_\text{hit}} \boldsymbol F_i.
```

Ubiased means that their averages are coincide with true averages:
```math
\langle P_\text{MC}(\text{hit}) \rangle = P(\text{hit}), \qquad
\bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigm| \text{hit} \bigr> = \langle \boldsymbol F \mid \text{hit} \rangle.
```

The accuracy
```math
\operatorname{var}(P_\text{MC}(\text{hit})) = 
```

```math
  \begin{align}
    \operatorname{cov}(\boldsymbol \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \mid \text{hit})
    &= \Bigl< \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \otimes \bigl[ \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} - \langle \boldsymbol F \mid \text{hit} \rangle \bigr] \Bigm| \text{hit} \Bigr> \\
    &= \bigl< \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \otimes \langle \boldsymbol F \mid \text{hit} \rangle_\text{MC} \bigm| \text{hit} \bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\
    &= \Bigl< \frac{1}{N_\text{hit}^2} \sum_{i = 1}^{N_\text{hit}} \sum_{j = 1}^{N_\text{hit}} \boldsymbol F_i \otimes \boldsymbol F_j \Bigm| \text{hit} \Bigr> - \langle \boldsymbol F \mid \text{hit} \rangle \otimes \langle \boldsymbol F \mid \text{hit} \rangle \\

  \end{align}
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