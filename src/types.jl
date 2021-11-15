# Types for Episiming.jl

"""
    primitive type Phase 8 end

`Phase` is a primitive type with 8 bits, representing the different phases, or stages,
of a disease, for a given individual. The following phases are defined as (global)
constants:
* `SUSCEPTIBLE = Phase(0)`
* `EXPOSED = Phase(1)`
* `INFECTED = Phase(2)`
* `ASYMPTOMATIC = Phase(3)`
* `RECOVERED = Phase(4)`
* `DECEASED = Phase(5)`
"""
primitive type Phase 8 end

Phase(x::Int) = reinterpret(Phase, UInt8(x))
Base.Int(x::Phase) = convert(Int, reinterpret(UInt8, x)) # needed for Base.show(io::IO, x::Phase)
Base.Broadcast.broadcastable(x::Phase) = Ref(x)
Base.zero(::Phase) = Phase(0)
Base.zero(::Type{Phase}) = Phase(0)

phaselist = (
    :SUSCEPTIBLE,
    :EXPOSED,
    :INFECTED,
    :ASYMPTOMATIC,
    :RECOVERED,
    :DECEASED
)

for (i, n) in enumerate(phaselist)
    @eval const $n = Phase($i - 1)
end

function Base.show(io::IO, x::Phase)
    if x in eval.(phaselist)
        if get(io, :compact, true)
            print(io, first(string(phaselist[Int(x) + 1])))
        else
            print(io, string(phaselist[Int(x) + 1]))
        end
    else
        print(io, "Undefined phase $x")
    end
end

Base.show(io::IO, ::MIME"text/plain", x::Phase) =
    x in eval.(phaselist) ?
        print(io, "Epidemic phase:\n   ", x) :
        print(io, "Undefined phase $x")

 #' Some colors for displaying the epidemics phase of the population

const phase_colors = Dict(
    SUSCEPTIBLE => :lightgray,
    EXPOSED => :orange,
    INFECTED => :red,
    ASYMPTOMATIC => :violet,
    RECOVERED => :blue,
    DECEASED => :black
)

#' ## Defining the main composite types
#'
#' ### Population

"""
    struct Population{R, S, T, U, V, W, X} <: AbstractVector{Tuple{{R, S, T, U, V, W, X}}}

Population structure with the following fields:
* `phase::Vector{R}`
* `past_transition::Vector{S}`
* `next_transition::Vector{Tuple{R, S}}`
* `residence::Vector{T}`
* `position::Vector{Tuple{U, U}}`
* `age::Vector{V}`
* `susceptibility::Vector{W}`
* `infectivity::Vector{W}`
* `clusters::Vector{Dict{X, T}}`
* `networks::Vector{Dict{X, T}}`

Usually, `R=Phase`, `S=Int`, `T=U=Int`, `U=W=Float64`, `X=Symbol`.
If is a huge population and memory is critical, one can set `S=Int16`, `U=Float16`,
`V=Int8`, `W=Float16`, `X=Int8`.
"""
struct Population{R, S, T, U, V, W, X} <: AbstractVector{Tuple{R, S, T, U, V, W, X}}
    phase::Vector{R}
    past_transition::Vector{S}
    next_transition::Vector{Tuple{R, S}}
    residence::Vector{T}
    position::Vector{Tuple{U, U}}
    age::Vector{V}
    susceptibility::Vector{W}
    infectivity::Vector{W}
    clusters::Vector{Dict{X, T}}
    networks::Vector{Dict{X, T}}
end

Base.size(population::Population) = size(population.phase)
Base.getindex(population::Population, n::Int) = (
    population.phase[n],
    population.past_transition[n],
    population.next_transition[n],
    population.residence[n],
    population.position[n],
    population.age[n],
    population.susceptibility[n],
    population.infectivity[n],
    population.clusters[n],
    population.networks[n]
)

function Base.setindex!(
    population::Population{R, S, T, U, V, W, X}, 
    v::Tuple{
        Vector{R},
        Vector{S},
        Vector{Tuple{R, S}},
        Vector{T},
        Vector{Tuple{U, U}},
        Vector{V},
        Vector{W},
        Vector{W},
        Vector{Dict{X, T}},
        Vector{Dict{X, T}}
    },
    n::Int
) where {R, S, T, U, V, W, X}
    population.phase[n] = v[1]
    population.past_transition[n] = v[2]
    population.next_transition[n] = v[3]
    population.residence[n] = v[4]
    population.position[n] = v[5]
    population.age[n] = v[6]
    population.susceptibility[n] = v[7]
    population.infectivity[n] = v[8]
    population.clusters[n] = v[9]
    population.networks[n] = v[10]
    return v
end

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
