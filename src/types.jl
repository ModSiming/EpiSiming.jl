#' # Types for Episiming.jl

#' ## Defining the individual epidemics state constants
#'
#' We define a primitive type which uses a single byte (8 bits) to reduce the memory usage.
#' That is the smallest possible size for a primitive type. Even Bools are 8 bits.
#' The only thing that uses less memory is BitVector, but that is a ... Vector, not a 
#' singleton.
#'
#' We could also use straight `Int8` to represent the states, but the nice thing about
#' defining a new type is that we can define how it is displayed. So, the SUSCEPTIBLE state
#' can be displayed as "SUSCEPTIBLE" or "S", instead of say "0". Imagine figuring out
#' whether "5" is RECOVERED, QUARENTINED or DECEASED.
#'
#' We make the primitive type a subtype of Number to make it easier for insertion on
#' Vectors (e.g. if `u` is a vector of State, we can do `u .= SUSCEPTIBLE`).

primitive type State <: Number 8 end

State(x::Int) = reinterpret(State, Int8(x))
Base.Int8(x::State) = reinterpret(Int8, x)
Base.Int(x::State) = convert(Int, Int8(x))

statelist = (
    :SUSCEPTIBLE,
    :EXPOSED,
    :INFECTED,
    :ASYMPTOMATIC,
    :RECOVERED,
    :QUARENTINED,
    :DECEASED
)

for (i, n) in enumerate(statelist)
    @eval const $n = State($i - 1)
end

Base.show(io::IO, x::State) =
    x in eval.(statelist) ?
        print(io, string(statelist[Int(x) + 1])) :
        print(io, "Undefined state $x")

function Base.show(io::IO, x::State)
    if x in eval.(statelist)
        if get(io, :compact, true)
            print(io, first(string(statelist[Int(x) + 1])))
        else
            print(io, string(statelist[Int(x) + 1]))
        end
    else
        print(io, "Undefined state $x")
    end
end

Base.show(io::IO, ::MIME"text/plain", x::State) =
    x in eval.(statelist) ?
        print(io, "Epidemic state:\n   ", x) :
        print(io, "Undefined state $x")



 #' Some colors for displaying the epidemics state of the population

const state_colors = Dict(
    SUSCEPTIBLE => :lightgray,
    EXPOSED => :orange,
    INFECTED => :red,
    ASYMPTOMATIC => :violet,
    RECOVERED => :blue,
    QUARENTINED => :yellow,
    DECEASED => :black
)

#' ## Defining the main composite types
#'
#' ### Population

struct Population{T, S} <: AbstractVector{Tuple{T, S}}
    state::Vector{S}
    event_history::Vector{T}
    residence::Vector{Int}
    position::Vector{Tuple{Float64, Float64}}
    age::Vector{Int8}
    susceptibility::Vector{Float64}
    infectivity::Vector{Float64}
    clusters::Vector{Dict{Symbol, Int}}
    networks::Vector{Dict{Symbol, Int}}
end

Base.size(population::Population) = size(population.state)
Base.getindex(population::Population, n::Int) = (
    population.state[n],
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
    population.state[n] = v[1]
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

struct Clusters
    id::Symbol
    name::String
    contact_rate::Float64
    clusters::Vector{Vector{Int}}
end

#' ### networks

#' ### Scenario

struct Scenario
    name::String
    info::String
    num_population::Int
    population::Population
    residences::Residences
    res_size_distribution::StatsBase.Weights{Float64, Float64, Vector{Float64}}
    pop_pyramid::StatsBase.Weights{Float64, Float64, Vector{Float64}}
    Γ_suscepibility::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
    Γ_infectivity::NamedTuple{(:shape, :scale), Tuple{Float64, Float64}}
    contact_rate::NamedTuple{(:residences, :general), Tuple{Float64, Float64}}
    recovery_rate::NamedTuple{
        (:exposed, :infected, :asymptomatic),
        Tuple{Float64, Float64, Float64}
    }
end
