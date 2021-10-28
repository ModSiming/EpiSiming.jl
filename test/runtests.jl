using EpiSiming

using Test
using Random
using StatsBase

Random.seed!(0)

@testset "Utils" begin
  include("random_scenario.jl")
end