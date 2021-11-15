#+ echo = false
#
# Weave it to markdown, for showing on github and in Documenter, with
# ```julia
# ]activate .
# using Weave
# cd("docs/src/examples")
# weave("random_scenario.jl", fig_path = "random_scenario_img", doctype = "github")
# ```
#
# Weave it to html with 
# `weave("random_scenario.jl", fig_path = "random_scenario_img")`

#' # Epidemic simulation via a discrete-time, agent-based stochastic model with a random scenario

#' ## Introduction
#'
#' We construct, here, the scenario for an epidemics, in which the population is composed
#' of a number of agents (or individuals), aggregated into residences, and interacting with
#' each other within the residences and through the following structures:
#' 1. clusters of school and workplaces,
#' 2. networks of social and commerce connections; and
#' 3. causual encounters with the general population.
#' 
#' In the current example, we build a random scenario, in a rectangular region.

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

num_population = 10_000

#' #### Region size

region_size = (6, 12)

#' #### Distribution of residences per size

res_size_distribution = [10, 22, 33, 22, 5, 5, 2, 1] |> v -> Weights(v / sum(v))

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

τ = Dict{Symbol, Float64}(
    :residences => 0.3,
    :work_places => 0.1,
    :school_places => 0.2
)

#' #### Recovery rates

γ = (
    rate_expos = 0.25, # rate out of exposed
    rate_infec = 0.1, # rate out of infected
    rate_asymp = 0.2 # rate out of asymptomatic
)

#' #### Fate probabilities

prob = (
    asymp = 0.6, # probability of becoming asymptomatic (vs. symptomatic = infected)
    decease = 0.02 # probability of deceasing
)


#' #### Random number generator
#'
#' Set random number generator for repeatability in testings

rng = MersenneTwister(123)

#' ### Generating scenario

#' #### Generating a random matrix of blocks with the population divided per block

pop_blocks = EpiSiming.gen_random_pop_blocks(rng, num_population, region_size)

#' Heatmap

display(
    heatmap(
        pop_blocks,
        c = cgrad([:white, :orange, :red]),
        annotate = [
            (j, i, text(pop_blocks[i, j], 8, :black, :center))
            for i in 1:size(pop_blocks, 1) for j in 1:size(pop_blocks, 2)
        ],
        xlabel = "km",
        ylabel = "km",
        title = "Distribution per block of the $num_population area residents",
        titlefont = 10
    )
)

#' #### Generating residences per block

res_blocks_distrib = EpiSiming.gen_res_blocks(pop_blocks, res_size_distribution)

#' Visualizations

display(
    groupedbar(
        [res_size_distribution sum(res_blocks_distrib) /
            sum(Iterators.flatten(res_blocks_distrib))],
        xticks = eachindex(res_size_distribution),
        bar_position = :dodge,
        bar_width = 0.7,
        label = ["target" "generated"],
        xlabel = "residence size",
        ylabel = "fraction of residences",
        title = "distribution of residences per size",
        titlefont = 10
    )
)

display(
    histogram(
        reshape(pop_blocks ./ map(sum, res_blocks_distrib), :, 1),
        bins = firstindex(res_size_distribution):0.1:lastindex(res_size_distribution),
        label = nothing,
        xlabel = "ratio of residents per residence",
        ylabel = "number of blocks",
        title = "Distribution of the number of blocks per average residence size" *
            "\nMean of the blocks average size: " * 
            "$(round(mean(pop_blocks ./ map(sum, res_blocks_distrib)), digits = 2))",
        titlefont = 10
    )
)

#' #### Generating the list of residences with all the required information

residences = EpiSiming.gen_residences(rng, res_blocks_distrib)

#' Visualizations

display(
    scatter(
        residences.position,
        markersize = 2,
        markerstrokecolor = :auto,
        label = false,
        xlabel = "km",
        ylabel = "km",
        title = "Residences",
        titlefont = 10
    )
)

display(
    scatter(
        filter( p -> (2.0 ≤ first(p) ≤ 3.0) && (2.0 ≤ last(p) ≤ 3.0), residences.position),
        markersize = 2,
        markerstrokecolor = :auto,
        label = false,
        xlabel = "km",
        ylabel = "km",
        title = "Residences",
        titlefont = 10
    )
)

display(
    scatter(
        filter( p -> (2.0 ≤ first(p) ≤ 2.2) && (2.0 ≤ last(p) ≤ 2.2), residences.position),
        markersize = 2,
        markerstrokecolor = :auto,
        label = false,
        xlabel = "km",
        ylabel = "km",
        title = "Residences",
        titlefont = 10
    )
)

#' #### Generating population

population = EpiSiming.gen_population(
    rng,
    residences,
    res_blocks_distrib,
    Γ,
    pop_pyramid
)

#' #### Susceptibility distribution

display(
    scatter(
        population.susceptibility,
        title = "Susceptibility (gamma distribution with mean $(Γ.sus_shape * Γ.sus_scale))",
        titlefont = 10,
        label = nothing,
        xlabel = "sample of population",
        ylabel = "susceptibility"
    )
)

let bins = 100, interval = extrema(population.susceptibility) |>  u -> u[2] - u[1]
    plt = histogram(
        population.susceptibility,
        bins = bins,
        label = "generated",
        title = "susceptibility histogram",
        titlefont = 10,
        xlabel = "susceptibility",
        ylabel = "population"
    )
    plot!(
        plt,
        x -> (num_population * interval / bins + 1) *
            Distributions.pdf(Distributions.Gamma(Γ.sus_shape, Γ.sus_scale), x),
            label = "target (Gamma distribution)"
    )
    display(plt)
end

#' #### Infectivity

display(
    scatter(
        population.infectivity,
        title = "Infectivity (gamma distribution with mean $(Γ.inf_shape * Γ.inf_scale))",
        titlefont = 10,
        label = nothing,
        xlabel = "sample of population",
        ylabel = "infectivity"
    )
)

let bins = 100, interval = extrema(population.infectivity) |>  u -> u[2] - u[1]
    plt = histogram(
        population.infectivity,
        bins = bins,
        label = "generated",
        title = "infectivity histogram",
        titlefont = 10,
        xlabel = "infectivity",
        ylabel = "population"
    )

    plot!(
        plt,
        x -> (num_population * interval / bins + 1) *
            Distributions.pdf(Distributions.Gamma(Γ.inf_shape, Γ.inf_scale), x),
        label = "target (Gamma distribution)"
    )
    display(plt)
end

#' #### Age distribution

plt = histogram(
    population.age,
    bins = length(pop_pyramid) + 1,
    label = "generated",
    title = "Age pyramid",
    titlefont = 10,
    xlabel = "age",
    ylabel = "individuals"
)

plot!(plt, 0:length(pop_pyramid) - 1, num_population * pop_pyramid, label = "target")
display(plt)

#' #### Agents location

display(
    scatter(
        population.position,
        c = map(s -> EpiSiming.phase_colors[s], population.phase),
        markersize = 1,
        markerstrokecolor = :auto,
        label = false,
        xlabel = "km",
        ylabel = "km",
        title = "Population",
        titlefont = 10
    )
)

display(
    scatter(
        filter(p -> (3.0 < first(p) < 4.0) & (1.0 < last(p) < 2.0), population.position),
        c = map(
            n -> EpiSiming.phase_colors[population.phase[n]],
            filter(
                n -> (3.0 < first(population.position[n]) < 4.0) &
                    (1.0 < last(population.position[n]) < 2.0),
                    eachindex(population.position)
            )
        ),
        markersize = 4,
        markerstrokecolor = :auto,
        xlabel = "km",
        ylabel = "km",
        label = false,
        title = "Population",
        titlefont = 10
    )
)

display(
    scatter(
        filter(p -> (3.2 < first(p) < 3.8) & (1.2 < last(p) < 1.8), population.position),
        c = map(
            n -> EpiSiming.phase_colors[population.phase[n]],
            filter(
                n -> (3.2 < first(population.position[n]) < 3.8) &
                    (1.2 < last(population.position[n]) < 1.8),
                    eachindex(population.position)
            )
        ),
        markersize = 4,
        markerstrokecolor = :auto,
        xlabel = "km",
        ylabel = "km",
        label = false,
        title = "Population",
        titlefont = 10
    )
)

display(
    scatter!(
        filter( p -> (3.2 ≤ first(p) ≤ 3.8) && (1.2 ≤ last(p) ≤ 1.8), residences.position),
        markersize = 4,
        markerstrokecolor = :auto,
        color = :green,
        label = false,
        xlabel = "km",
        ylabel = "km",
        title = "Residences",
        titlefont = 10
    )
)

#' ### Generate clusters

max_size = 100
α = 1.8

clusters = Dict{Symbol, Vector{Vector{Int}}}()

push!(
    clusters,
    :work_places => EpiSiming.gen_clusters(
        rng,
        filter(n -> population.age[n] ≥ 18, 1:num_population),
        max_size,
        α
    )
)

push!(
    clusters,
    :school_places => EpiSiming.gen_clusters(
        rng,
        filter(n -> population.age[n] < 20, 1:num_population),
        max_size,
        α
    )
)

for (k, v) in clusters
    for (i, r) in enumerate(v)
        for n in r
            push!(population.clusters[n], k => i)
        end
    end
end

#' Visualizations

plt = histogram(
    length.(clusters[:work_places]),
    bins = max_size,
    xlabel = "work size",
    ylabel = "number of work places",
    label = "generated",
    title = "Distribution of work places according to size",
    titlefont = 10
)
plot!(
    plt,
    length(clusters[:work_places]) * EpiSiming.decay.(1:max_size, α) / sum(EpiSiming.decay.(1:max_size, α)),
    linewidth = 2,
    label = "target"
)
display(plt)

#' ## Evolution of the epidemics

#' ### Initial infection

num_exposed_at_time_0 = div(num_population, 500) # ( 1/500 = 0.002 = 0.2%) 20 for 10_000

exposed_at_time_0 = sample(rng, 1:num_population, num_exposed_at_time_0, replace = false)

# population.phase[exposed_at_time_0] .= EXPOSED
for n in exposed_at_time_0
    population.phase[n] = EXPOSED
    next_change, next_phase = EpiSiming.transition_rules(rng, EXPOSED, 1)
    population.next_transition[n] = (next_phase, next_change)
end

#' ### Evolution parameters

num_steps = 360
time_step = 1

#' ### Time evolution

@time evolution = EpiSiming.evolve!(
    rng,
    population,
    residences,
    clusters,
    τ,
    γ,
    prob,
    num_steps,
    time_step,
    verbose_step = 10
)

@time summary = EpiSiming.get_summary(evolution)

#' Visualizations

display(
    plot(
        summary[:, 2:end],
        labels = string.(reduce(hcat, EpiSiming.phaselist[2:end])),
        legend = :left,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
)

display(
    plot(
        summary[:, [2, 3, 4, 6]],
        labels = string.(reduce(hcat, EpiSiming.phaselist[[2, 3, 4, 6]])),
        legend = :right,
        xlabel = "day",
        ylabel = "cases",
        title = "time evolution of cases",
        titlefont = 10
    )
)
