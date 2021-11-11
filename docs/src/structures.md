```@meta
CurrentModule = EpiSiming
```

# Structures

## Epidemic phases for an individual

We define a type to hold the epidemic phase of an individual (e.g. asymptotic, exposed, infected, recovered, and so on). This type is defined as a primitive type which uses a single byte (8 bits), in order to reduce the memory usage. That is the smallest possible size for a primitive type. Even Bools are 8 bits. The only information entity that uses less memory is BitVector, but that is a Vector, not a singleton, and is not very performant for random access; it is more suitable for sequential access.

We could have used straight `Int8` to represent the phases of the disease, but defining a new type has the advantage of allowing us to overload `Base.show` to display the phase in a more friendly and recognizable way. So, the SUSCEPTIBLE state can be displayed as "SUSCEPTIBLE" or "S", instead of say "0". Imagine figuring out whether "5" is RECOVERED or DECEASED.

```@docs
Phase
```

## Residences

```@docs
Residences
```

## Population

```@docs
Population
```
## Clusters

```@docs
Clusters
```

## Networks

```@docs
Networks
```

## Scenario

```@docs
Scenario
```
