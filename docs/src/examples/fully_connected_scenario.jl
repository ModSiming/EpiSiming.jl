#+ echo = false
#
# Weave it to markdown, for showing on github and in Documenter, with
# ```julia
# ]activate .
# using Weave
# cd("docs/src/examples")
# weave("fully_connected_scenario.jl", fig_path = "fully_connected_scenario_img", doctype = "github")
# ```
#
# Weave it to html with 
# `weave("fully_connected_scenario.jl", fig_path = "fully_connected_scenario_img")`

#' # Epidemic simulation via a discrete-time, agent-based stochastic model with a fully-connected scenario

#' ## Introduction
#'
#' We construct, here, a fully-connected scenario for an epidemics, in which the population is composed
#' of a number of agents (or individuals), aggregated into individual residences, and interacting with
#' each other through a fully-connected newtork:

#' ## Loading packages

using Random

using SparseArrays

using StatsBase
using Distributions

using Plots
using StatsPlots: groupedbar

using EpiSiming

#=
using Graphs
using Graphs
using GraphPlot
=#

#' ## Building the Scenario

#' ### Set global parameters

#' #### Population number

num_population = 1_000

#' #### Region size

region_size = (1, 1)

#' #### Distribution of residences per size

res_size_distribution = [1.0]

#' #### Susceptibility and infectivity parameters

#' Parameters of the Gamma distribution for the susceptibility and infectivity factors

Γ = (
    sus_shape = 2.0,
    sus_scale = 0.5,
    inf_shape = 4.0,
    inf_scale= 0.25
)

#' #### Age pyramid

pop_pyramid = let age_max = 100, p = 2, pyramid_func(a, age_max, p) = ( a + 1 ) * (age_max - a)^p
    Weights(
        pyramid_func.(0:age_max, age_max, p) / sum(pyramid_func.(0:age_max, age_max, p))
    )
end

#' #### Contact rates

τ = Dict{Symbol, Float64}(:complete => 1.0)

#' #### Random number generator
#'
#' Set random number generator for repeatability in testings

rng = MersenneTwister(123)

#' ### Generating scenario

#' #### Generating the single-block with the whole population

pop_blocks = [num_population]

#' #### Generating residences per residents per block

res_blocks_distrib = fill([num_population], 1, 1)
#' #### Generating the list of residences with all the required information

residences = EpiSiming.gen_residences(rng, res_blocks_distrib)

#' #### Generating population

population = EpiSiming.gen_population(
    rng,
    residences,
    res_blocks_distrib,
    Γ,
    pop_pyramid
)

#' ### Generate complete cluster

clusters = Dict{Symbol, Vector{Vector{Int}}}(:complete => [collect(1:num_population)])

#' ## Evolution of the epidemics

#' ### Evolution parameters

num_steps = 90
time_step = 1

#' ### Initial infection

num_exposed_at_time_0 = 2

exposed_at_time_0 = sample(rng, 1:num_population, num_exposed_at_time_0, replace = false)

population.phase .= SUSCEPTIBLE
for n in exposed_at_time_0
    population.phase[n] = EXPOSED
    next_phase, next_change = EpiSiming.transition_rules(rng, EXPOSED, 1)
    population.next_transition[n] = next_change
    population.next_phase[n] = next_phase
end

#' ### Time evolution

@time evolution = EpiSiming.evolve!(
    rng,
    population,
    residences,
    clusters,
    τ,
    num_steps,
    time_step,
    verbose_step = 10
)

@time summary = EpiSiming.get_summary(evolution)

#' Visualizations

display(
    plot(
        summary[:, 1:end],
        labels = string.(reduce(hcat, EpiSiming.PHASE_LIST[2:end-1])),
        legend = :left,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
)

#= display(
    plot(
        summary[:, [3, 4, 5, 7]],
        labels = string.(reduce(hcat, EpiSiming.PHASE_LIST[[3, 4, 5, 7]])),
        legend = :right,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
)
 =#
#' Solving with transitions

#' Reinitializing population

population.phase .= SUSCEPTIBLE
for n in exposed_at_time_0
    population.phase[n] = EXPOSED
    next_phase, next_change = EpiSiming.transition_rules(rng, EXPOSED, 1)
    population.next_transition[n] = next_change
    population.next_phase[n] = next_phase
end

@time transitions = EpiSiming.evolve_for_transitions!(
    rng,
    population,
    residences,
    clusters,
    τ,
    num_steps,
    time_step,
    verbose_step = 10
)

@time summary_trsn = EpiSiming.get_summary_from_transitions(transitions)

#' Visualizations

display(
    plot(
        summary_trsn[:, 1:end],
        labels = string.(reduce(hcat, EpiSiming.PHASE_LIST[2:end-1])),
        legend = :left,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
)
#= 
display(
    plot(
        summary_trsn[:, [3, 4, 5, 7]],
        labels = string.(reduce(hcat, EpiSiming.PHASE_LIST[[3, 4, 5, 7]])),
        legend = :right,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
) =#
