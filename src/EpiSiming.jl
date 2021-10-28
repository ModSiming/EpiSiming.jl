module EpiSiming

using Random

using SparseArrays

using StatsBase
using Distributions

using ProgressMeter

export evolve!
export SUSCEPTIBLE, EXPOSED, INFECTED, ASYMPTOMATIC, RECOVERED, QUARENTINED, DECEASED

#=
using Graphs
using MetaGraphs
=#

include("types.jl")
include("evolution.jl")
include("scenarios/random_scenario.jl")

end
