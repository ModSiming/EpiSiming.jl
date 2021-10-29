# Types

## Types for the stages of the epidemics

We define a primitive type which uses a single byte (8 bits) to reduce the memory usage. That is the smallest possible size for a primitive type. Even Bools are 8 bits. The only thing that uses less memory is BitVector, but that is a ... Vector, not a singleton.

We could also use straight `Int8` to represent the states, but the nice thing about defining a new type is that we can define how it is displayed. So, the SUSCEPTIBLE state can be displayed as "SUSCEPTIBLE" or "S", instead of say "0". Imagine figuring out whether "5" is RECOVERED, QUARENTINED or DECEASED.

We make the primitive type a subtype of Number to make it easier for insertion on Vectors (e.g. if `u` is a vector of State, we can do `u .= SUSCEPTIBLE`).

```@docs
EpiSiming.State
```