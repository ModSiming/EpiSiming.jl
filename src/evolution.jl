#' # Evolution

#' ## Force of infection

"""
Compute the force of infection for each agent, i.e. return a vector `λ` where each
`λ[i]`` is the force of infection on the agent with index `i`.

Since we expect the number of susceptibles to be nearly the whole population in most cases,
we don't care to filter by susceptible, just leaving it zero for non-susceptibles.
"""
function force_of_infection!(λ, population, residences, clusters, τ)
    for n in eachindex(population)
        if population.state[n] == SUSCEPTIBLE
            res = population.residence[n]
            number_of_coresidents = residences.num_residents[res] - 1
            λ[n] = 0.0
            if number_of_coresidents > 0
                for k in residences.residents[res]
                    if population.state[k] == INFECTED || 
                        population.state[k] == ASYMPTOMATIC
                        λ[n] += population.infectivity[k]
                    end
                end
                λ[n] *= τ[:residences] / number_of_coresidents
            end
            for (cluster, index) in population.clusters[n]
                cluster_members = clusters[cluster][index]
                num_of_cluster_members = length(cluster_members) - 1
                if num_of_cluster_members > 1 # should not even include clusters with a single member
                    cocluster_infectivity = 0.0
                    for k in cluster_members
                        if population.state[k] == INFECTED || 
                            population.state[k] == ASYMPTOMATIC
                            cocluster_infectivity += population.infectivity[k]
                        end
                    end
                    λ[n] += τ[cluster] * cocluster_infectivity / num_of_cluster_members
                end
            end
        end

    end
    return λ
end

#' ## Single step forward

function step_foward!(rng, population, chances, λ, γ, prob, k)
    rand!(rng, chances)
    for n in eachindex(population)
        state = population.state[n]
        if state == SUSCEPTIBLE && chances[n] ≤ 1 - exp(-λ[n])
            population.event_history[n] = k
            population.state[n] = EXPOSED
        elseif state == EXPOSED && chances[n] ≤ γ.rate_expos
            population.event_history[n] = k
            if rand(rng) ≤ prob.asymp
                population.state[n] = ASYMPTOMATIC
            else
                population.state[n] = INFECTED
            end
        elseif state == ASYMPTOMATIC && chances[n] ≤ γ.rate_asymp
            population.event_history[n] = k
            population.state[n] = RECOVERED
        elseif state == INFECTED && chances[n] ≤ γ.rate_infec
            population.event_history[n] = k
            if rand(rng) ≤ prob.decease
                population.state[n] = DECEASED
            else
                population.state[n] = RECOVERED
            end
        end
    end
    nothing
end

#' ## Time evolution

#' Main evolution function

function evolve!(rng, population, residences, clusters, τ, γ, prob, num_steps, time_step;
    verbose::Integer = false
)
    num_population = length(population)
    evolution = spzeros(State, num_population, num_steps)
    evolution[:, 1]  .= population.state # getfield(population, :state)
    λ = Vector{Float64}(undef, num_population)
    chances = Vector{Float64}(undef, num_population)
    for k in 2:num_steps
        force_of_infection!(λ, population, residences, clusters, τ)
        step_foward!(rng, population, chances, λ, γ, prob, k)
        for n in 1:num_population
            state = population.state[n]
            if state != SUSCEPTIBLE 
                evolution[n, k]  = state
            end
        end
        if verbose != false && mod(k, verbose) == 0
            @info "Done time step $k (day $(k * time_step))"
        end
    end
    return evolution
end

#' Summary with the total population in each compartiment

function get_summary(evolution)
    num_population, num_steps = size(evolution)
    summary = zeros(Int, num_steps, 7)
    for k in 1:num_steps
        for n in 1:num_population
            summary[k, Int(evolution[n, k]) + 1] += 1
        end
    end
    return summary
end