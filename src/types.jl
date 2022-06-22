# Types for Episiming.jl
#= 
"""
    primitive type Phase 8 end

`Phase` is a primitive type with 8 bits, representing the different phases, or stages,
of a disease, for a given individual. The following phases are defined as (global)
constants:
* `SUSCEPTIBLE = Phase(1)`
* `EXPOSED = Phase(2)`
* `INFECTED = Phase(3)`
* `ASYMPTOMATIC = Phase(4)`
* `RECOVERED = Phase(5)`
* `DECEASED = Phase(6)`
"""
primitive type Phase 8 end

Phase(x::Int) = reinterpret(Phase, UInt8(x))
Base.Int(x::Phase) = convert(Int, reinterpret(UInt8, x)) # needed for Base.show(io::IO, x::Phase)
Base.Broadcast.broadcastable(x::Phase) = Ref(x)
Base.zero(::Phase) = Phase(1)
Base.zero(::Type{Phase}) = Phase(1) =#

"""
    struct Phase

`Phase` is a composite type with a single 8 bits UInt8 value, representing the different phases, or stages,
of a disease, for a given individual. 

The following phases are defined as (global)
constants:
# `NULL = Phase(0)`
* `SUSCEPTIBLE = Phase(1)`
* `EXPOSED = Phase(2)`
* `INFECTED = Phase(3)`
* `ASYMPTOMATIC = Phase(4)`
* `RECOVERED = Phase(5)`
* `DECEASED = Phase(6)`
* `UNKNOWN = Phase(7)`

PS: I have tried with a primitive type Phase 8 end but the footprint and performance are the same, so a simple struct is easier and, in fact,
recommended.
"""
struct Phase
    val::UInt8
end

Phase(x::Int) = Phase(UInt8(x))
Base.Int(x::Phase) = convert(Int, x.val)
Base.Broadcast.broadcastable(x::Phase) = Ref(x)
Base.zero(::Phase) = Phase(0)
Base.zero(::Type{Phase}) = Phase(0)

PHASE_LIST = (
    :NULL,
    :SUSCEPTIBLE,
    :EXPOSED,
    :INFECTED,
    :ASYMPTOMATIC,
    :RECOVERED,
    :DECEASED,
    :UNKNOWN
)

for (i, n) in enumerate(PHASE_LIST)
    @eval const $n = Phase($i-1)
end

function Base.show(io::IO, x::Phase)
    if x in eval.(PHASE_LIST)
        if get(io, :compact, true)
            print(io, first(string(PHASE_LIST[Int(x)+1])))
        else
            print(io, string(PHASE_LIST[Int(x)+1]))
        end
    else
        print(io, "Undefined phase")
    end
end

Base.show(io::IO, ::MIME"text/plain", x::Phase) =
    x in eval.(PHASE_LIST) ?
        print(io, "Epidemic phase:\n   ", x) :
        print(io, "Undefined phase")

 #' Some colors for displaying the epidemics phase of the population

const phase_colors = Dict(
    NULL => :black,
    SUSCEPTIBLE => :lightgray,
    EXPOSED => :orange,
    INFECTED => :red,
    ASYMPTOMATIC => :violet,
    RECOVERED => :blue,
    DECEASED => :brown,
    UNKNOWN => :white
)

#' ## Defining the main composite types
#'
#' ### Population

"""
    struct Population{R, S, T, U, V, W} <: AbstractVector{Tuple{{R, S, T, U, V, W}}}

Population structure with the following fields:
* `phase::Vector{R}`
* `previous_transition::Vector{S}`
* `next_transition::Vector{S}`
* `next_phase::Vector{R}`
* `residence::Vector{T}`
* `position::Vector{Tuple{U, U}}`
* `age::Vector{V}`
* `susceptibility::Vector{W}`
* `infectivity::Vector{W}`

Usually, `R=Phase`, `S=Int`, `T=U=Int`, `U=W=Float64`, `X=Symbol`.
If is a huge population and memory is critical, one can set `S=Int16`, `U=Float16`,
`V=Int8`, `W=Float16`.
"""
struct Population{R, S, T, U, V, W} <: AbstractVector{Tuple{R, S, T, U, V, W}}
    phase::Vector{R}
    previous_transition::Vector{S}
    next_transition::Vector{S}
    next_phase::Vector{R}
    residence::Vector{T}
    position::Vector{Tuple{U, U}}
    age::Vector{V}
    susceptibility::Vector{W}
    infectivity::Vector{W}
end

Base.size(population::Population) = size(population.phase)
Base.IndexStyle(::Type{<:Population}) = IndexLinear()
Base.getindex(population::Population, n::Int) = (
    population.phase[n],
    population.previous_transition[n],
    population.next_transition[n],
    population.next_phase[n],
    population.residence[n],
    population.position[n],
    population.age[n],
    population.susceptibility[n],
    population.infectivity[n]
)

function Base.setindex!(
    population::Population{R, S, T, U, V, W}, 
    v::Tuple{
        Vector{R},
        Vector{S},
        Vector{S},
        Vector{R},
        Vector{T},
        Vector{Tuple{U, U}},
        Vector{V},
        Vector{W},
        Vector{W}
    },
    n::Int
) where {R, S, T, U, V, W}
    population.phase[n] = v[1]
    population.previous_transition[n] = v[2]
    population.next_transition[n] = v[3]
    population.next_phase[n] = v[4]
    population.residence[n] = v[5]
    population.position[n] = v[6]
    population.age[n] = v[7]
    population.susceptibility[n] = v[8]
    population.infectivity[n] = v[9]
    return v
end

Base.similar(population::Population{R, S, T, U, V, W}) where {R, S, T, U, V, W} = Population{R, S, T, U, V, W}(
    Vector{R}(undef, length(population.phase)),
    Vector{S}(undef, length(population.previous_transition)),
    Vector{S}(undef, length(population.next_transition)),
    Vector{R}(undef, length(population.next_phase)),
    Vector{T}(undef, length(population.residence)),
    Vector{Tuple{U, U}}(undef, length(population.position)),
    Vector{V}(undef, length(population.age)),
    Vector{W}(undef, length(population.susceptibility)),
    Vector{W}(undef, length(population.infectivity))
)

# It shouldn't be needed to extend copy; it should work automatically
# from the above, but something is not okay
Base.copy(population::EpiSiming.Population) = EpiSiming.Population(population.phase, population.previous_transition, population.next_transition, population.next_phase, population.residence, population.position, population.age, population.susceptibility, population.infectivity)

#' ### Residences

#' We parametrize it with types `T` and `S`, but with `Int` and `Float64` in mind.
#'

"""
    Residences{T, S} <: AbstractVector{Tuple{T, S}}

Residences structure with the following fields:
* `block::Vector{T}`
* `position::Vector{Tuple{S, S}}`
* `num_residents::Vector{T}`
* `residents::Vector{Vector{T}}`

Usually, `T` is an `Int` and `S` is a floating point type.
"""
struct Residences{T, S} <: AbstractVector{Tuple{T, S}}
    block::Vector{T}
    position::Vector{Tuple{S, S}}
    num_residents::Vector{T}
    residents::Vector{Vector{T}}
end

Base.size(residence::Residences) = size(residence.block)
Base.IndexStyle(::Type{<:Residences}) = IndexLinear()
Base.getindex(residence::Residences, n::Int) = (
    residence.block[n],
    residence.position[n],
    residence.num_residents[n],
    residence.residents[n]
)

function Base.setindex!(
    residence::Residences{T, S}, 
    v::Tuple{Vector{T}, Vector{Tuple{S, S}}, Vector{T}, Vector{Vector{T}}},
    n::Int
) where {T, S}
    residence.block[n] = v[1]
    residence.position[n] = v[2]
    residence.num_residents[n] = v[3]
    residence.residents[n] = v[4]
    return v
end

Base.similar(residences::Residences{T, S}) where {T, S} = Residences{T, S}(
    Vector{T}(undef, length(residences.block)),
    Vector{Tuple{T, S}}(undef, length(residences.positions)),
    Vector{T}(undef, length(residences.num_residents)),
    Vector{Vector{T}}(undef, length(residences.residents))
)


#' ### Clusters

"""
    struct Clusters{R, S, T, U}

Clusters structure with the following fields:
* `id::R`
* `name::S`
* `contact_rate::T`
* `clusters::Vector{Vector{U}}`
"""
struct Clusters{R, S, T, U}
    id::R
    name::S
    contact_rate::T
    clusters::Vector{Vector{U}}
end

#' ### networks

"""
struct Networks

Structure for the set of networks
"""
struct Networks
end
#' ### Scenario

"""
Structure for the Scenario, with the following fields:

* `name::String`
* `info::String`
* `num_population::Int`
* `population::Population`
* `residences::Residences`
* `res_size_distribution::StatsBase.Weights{Float64, Float64, Vector{Float64}}`
* `pop_pyramid::StatsBase.Weights{Float64, Float64, Vector{Float64}}`
* `Γ_susceptibility::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}`
* `Γ_infectivity::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}`
* `contact_rate::NamedTuple{(:residences, :general), Tuple{Float64, Float64}}`
* `recovery_rate::NamedTuple{
        (:exposed, :infected, :asymptomatic),
        Tuple{Float64, Float64, Float64}
  }`
"""
struct Scenario
    name::String
    info::String
    num_population::Int
    population::Population
    residences::Residences
    res_size_distribution::StatsBase.Weights{Float64, Float64, Vector{Float64}}
    pop_pyramid::StatsBase.Weights{Float64, Float64, Vector{Float64}}
    Γ_susceptibility::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
    Γ_infectivity::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
    contact_rate::NamedTuple{(:residences, :general), Tuple{Float64, Float64}}
    recovery_rate::NamedTuple{
        (:exposed, :infected, :asymptomatic),
        Tuple{Float64, Float64, Float64}
    }
end
