#' # Evolution

#' ## Force of infection

"""
Compute the force of infection on each individual, i.e. return a vector `λ` where each
`λ[i]` is the force of infection on the individual with index `i`.

The force of infection `λ[i]` is computed as the sum of the force of infection in each
group of direct contacts of the individual, be it their residence, their clusters
(e.g. school, workplace), their network neighbors (e.g. friends, neighbors,
local shops), and random contact (e.g. commuting, sporadical contacts).

For each group say `i_1, ..., i_n`, the associated force of contact is computed as
the sum of the infectivity rate for each of the infected individuals, multiplied by the
transmission rate in that group and divided by the number `n-1` of co-members in the group.
"""
function force_of_infection!(λ, population, residences, clusters, τ)

    for n in eachindex(population)
        λ[n] = 0.0
        if population.phase[n] == SUSCEPTIBLE
            res = population.residence[n]
            number_of_coresidents = residences.num_residents[res] - 1
            if number_of_coresidents ≥ 1
                for k in residences.residents[res]
                    if population.phase[k] == INFECTED || 
                        population.phase[k] == ASYMPTOMATIC
                        λ[n] += population.infectivity[k]
                    end
                end
                λ[n] *= τ[:residences] / number_of_coresidents
            end
        end
    end
    for cluster in keys(clusters)
        for ind in eachindex(clusters[cluster])
            cluster_members = clusters[cluster][ind]
            num_of_cluster_comembers = length(cluster_members) - 1
            if num_of_cluster_comembers ≥ 1 # should not even include clusters with a single member
                cocluster_infectivity = 0.0
                for n in cluster_members
                    if population.phase[n] == INFECTED || 
                        population.phase[n] == ASYMPTOMATIC
                        cocluster_infectivity += population.infectivity[n]
                    end
                end
                for n in cluster_members
                    λ[n] += τ[cluster] * cocluster_infectivity / num_of_cluster_comembers
                end
            end
        end
    end
    return λ
end

#' ## Single step forward

"""
step_forward!(rng, population, chances, λ, k)

A single step forward, based on the given force of infection `λ` (already computed with
[`force_of_infection!`](@ref)).

`chances` is preallocated and is mutated each step with new random numbers.

The phase transitions from susceptible to exposed depends on the success of

    `chances[n] ≤ 1 - exp(-λ[n] * population.susceptibility[n])`

The transition from the other states are given by [`transition_rules`](@ref).
"""
function step_forward!(rng, population, chances, λ, k)
    rand!(rng, chances)
    for n in eachindex(population)
        phase = population.phase[n]
        if phase == SUSCEPTIBLE
            if @fastmath chances[n] ≤ 1 - exp(-λ[n] * population.susceptibility[n])
                population.previous_transition[n] = k
                population.phase[n] = EXPOSED
                next_phase, next_change = transition_rules(rng, EXPOSED, k)
                population.next_transition[n] = next_change
                population.next_phase[n] = next_phase
            end
        elseif k ≥ population.next_transition[n]
            population.previous_transition[n] = k
            population.phase[n] = population.next_phase[n]
            next_phase, next_change = transition_rules(rng, phase, k)
            population.next_transition[n] = next_change
            population.next_phase[n] = next_phase
        end
    end
    nothing
end

#' ## Time evolution

#' Main evolution functions

"""
    evolve!(
        rng, population, residences, clusters, τ, num_steps, time_step;
        verbose_step::Integer = 0
    )

It essentially loops through [`step_forward!`](@ref) as many as `num_steps` times, and
returns an array with the transition phases ocurring for each individual.
"""
function evolve!(
    rng, population, residences, clusters, τ, num_steps, time_step;
    verbose_step::Integer = 0
)

    num_population = length(population)
    evolution = Matrix{Phase}(undef, num_population, num_steps)
    evolution[:, 1] .= population.phase
    λ = Vector{Float64}(undef, num_population)
    chances = Vector{Float64}(undef, num_population)
    
    for k in 2:num_steps
        force_of_infection!(λ, population, residences, clusters, τ)
        step_forward!(rng, population, chances, λ, k)
        evolution[:, k] .= population.phase
        if verbose_step > 0 && mod(k, verbose_step) == 0
            @info "Done time step $k (day $(k * time_step))"
        end
    end
    return evolution
end

#' Summary with the total population in each compartment

"""
    get_summary(evolution)

With `m` phases (SUSCEPTIBLE, EXPOSED, etc.) and `num_steps` iterations,
`get_summary(evolution)` returns a `num_steps x m` matrix of Integers with the population
count at each iteration, on each phase.
"""
function get_summary(evolution)
    num_population, num_steps = size(evolution)
    summary = zeros(Int, num_steps, 7)
    for k in 1:num_steps
        for n in 1:num_population
            summary[k, Int(evolution[n, k])] += 1
        end
    end
    return summary
end

"""
    evolve_for_transitions!(
        rng, population, residences, clusters, τ, num_steps, time_step;
        verbose_step::Integer = 0
    )

It essentially loops through [`step_forward!`](@ref) as many as `num_steps` times, and
returns a sparse array with the transition phases ocurring for each individual.
"""
function evolve_for_transitions!(
    rng, population, residences, clusters, τ, num_steps, time_step;
    verbose_step::Integer = 0
)

    num_population = length(population)
    transitions = spzeros(Phase, num_population, num_steps)
    transitions[:, 1]  .= population.phase
    current = copy(population.phase)
    λ = Vector{Float64}(undef, num_population)
    chances = Vector{Float64}(undef, num_population)
    
    for k in 2:num_steps
        force_of_infection!(λ, population, residences, clusters, τ)
        step_forward!(rng, population, chances, λ, k)
        for n in 1:num_population
            phase = population.phase[n]
            if phase != current[n]
                transitions[n, k] = phase
                current[n] = phase
            end
        end
        if verbose_step > 0 && mod(k, verbose_step) == 0
            @info "Done time step $k (day $(k * time_step))"
        end
    end
    return transitions
end

#' Summary with the total population in each compartment from
#' a sparse matrix of transitions (the transition matrix is yet 
#' to be implemented, but should be the way to go).

"""
    get_summary_from_transitions(transitions)

With `m` phases (SUSCEPTIBLE, EXPOSED, etc.) and `num_steps` iterations,
`get_summary(transitions)` returns a `num_steps x m` matrix of Integers with the population
count at each iteration, on each phase.
"""
function get_summary_from_transitions(transitions::SparseMatrixCSC)
    num_population, num_steps = size(transitions)
    summary = zeros(Int, num_steps, 7)
    current_state = Vector(@view(transitions[:, 1]))
    for kj in 1:num_steps
        m = transitions.colptr[kj]
        for ni in 1:num_population
            if m < transitions.colptr[kj+1] && ni == transitions.rowval[m] && !iszero(transitions.nzval[m])
                current_state[ni] = transitions.nzval[m]    
                m += 1
            end
            summary[kj, Int(current_state[ni])] += 1
        end
    end
    return summary
end
