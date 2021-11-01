# Types for Episiming.jl

"""
    primitive type Phase <: Number 8 end

`Phase` is a primitive subtype of `Number` with 8 bits, representing the different
phases, or stages, of a disease, for a given individual. The following phases are
defined as (global) constants:
* SUSCEPTIBLE = Phase(0)
* EXPOSED = Phase(1)
* INFECTED = Phase(2)
* ASYMPTOMATIC = Phase(3)
* RECOVERED = Phase(4)
* DECEASED = Phase(5)
"""
primitive type Phase <: Number 8 end

Phase(x::Int) = reinterpret(Phase, Int8(x))
Base.Int8(x::Phase) = reinterpret(Int8, x)
Base.Int(x::Phase) = convert(Int, Int8(x))

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
    struct Population{T, S} <: AbstractVector{Tuple{T, S}}

Population structure with the following fields:
* phase::Vector{S}
* event_history::Vector{T}
* residence::Vector{Int}
* position::Vector{Tuple{Float64, Float64}}
* age::Vector{Int8}
* susceptibility::Vector{Float64}
* infectivity::Vector{Float64}
* clusters::Vector{Dict{Symbol, Int}}
* networks::Vector{Dict{Symbol, Int}}
"""
struct Population{T, S} <: AbstractVector{Tuple{T, S}}
    phase::Vector{S}
    event_history::Vector{T}
    residence::Vector{Int}
    position::Vector{Tuple{Float64, Float64}}
    age::Vector{Int8}
    susceptibility::Vector{Float64}
    infectivity::Vector{Float64}
    clusters::Vector{Dict{Symbol, Int}}
    networks::Vector{Dict{Symbol, Int}}
end

Base.size(population::Population) = size(population.phase)
Base.getindex(population::Population, n::Int) = (
    population.phase[n],
    population.event_history[n],
    population.residence[n],
    population.position[n],
    population.age[n],
    population.susceptibility[n],
    population.infectivity[n],
    population.clusters[n],
    population.networks[n]
)

function Base.setindex!(
    population::Population{T, S}, 
    v::Tuple{
        Vector{S}, Vector{T}, Vector{Int}, Vector{Tuple{Float64, Float64}},
        Vector{Int8}, Vector{Float64}, Vector{Float64},
        Vector{Dict{Symbol, Int}}, Vector{Dict{Symbol, Int}}
    },
    n::Int
) where {T, S}
    population.phase[n] = v[1]
    population.event_history[n] = v[2]
    population.residence[n] = v[3]
    population.position[n] = v[4]
    population.age[n] = v[5]
    population.susceptibility[n] = v[6]
    population.infectivity[n] = v[7]
    population.clusters[n] = v[8]
    population.networks[n] = v[9]
    return v
end

#' ### Residences

#' We parametrize it with types `T` and `S`, but with `Int` and `Float64` in mind.
#'

"""
    Residences{T, S} <: AbstractVector{Tuple{T, S}}

Residences structure with the following fields:
* block::Vector{T}
* position::Vector{Tuple{S, S}}
* num_residents::Vector{T}
* residents::Vector{Vector{T}}
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
    struct Clusters

Clusters structure with the following fields:
* id::Symbol
* name::String
* contact_rate::Float64
* clusters::Vector{Vector{Int}}
"""
struct Clusters
    id::Symbol
    name::String
    contact_rate::Float64
    clusters::Vector{Vector{Int}}
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

* name::String
* info::String
* num_population::Int
* population::Population
* residences::Residences
* res_size_distribution::StatsBase.Weights{Float64, Float64, Vector{Float64}}
* pop_pyramid::StatsBase.Weights{Float64, Float64, Vector{Float64}}
* Γ_susceptibility::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
* Γ_infectivity::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
* contact_rate::NamedTuple{(:residences, :general), Tuple{Float64, Float64}}
* recovery_rate::NamedTuple{
      (:exposed, :infected, :asymptomatic),
      Tuple{Float64, Float64, Float64}
  }
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
