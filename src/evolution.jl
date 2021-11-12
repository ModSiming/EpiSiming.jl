#' # Evolution

#' ## Force of infection

"""
Compute the force of infection on each individual, i.e. return a vector `λ` where each
`λ[i]` is the force of infection on the individual with index `i`.

The force of infection `λ[i]` is computed as the sum of the force of infection in each
group of direct contacts of the individual, be it their residence, their clusters
(e.g. school, workplace), their network neighbors (e.g. friends, neighbors,
local shops), and random contact (e.g. commuting, sporadical contacts).

For each group say `i₁, …, iₙ`, the associated force of contact is computed as
the sum of the infectivity rate for each of the infected individuals, multiplied by the
transmission rate in that group and divided by the number `n-1` of co-members in the group.
"""
function force_of_infection!(λ, population, residences, clusters, τ)
    for n in eachindex(population)
        if population.phase[n] == SUSCEPTIBLE
            res = population.residence[n]
            number_of_coresidents = residences.num_residents[res] - 1
            λ[n] = 0.0
            if number_of_coresidents ≥ 1
                for k in residences.residents[res]
                    if population.phase[k] == INFECTED || 
                        population.phase[k] == ASYMPTOMATIC
                        λ[n] += population.infectivity[k]
                    end
                end
                λ[n] *= τ[:residences] / number_of_coresidents
            end
            for (cluster, index) in population.clusters[n]
                cluster_members = clusters[cluster][index]
                num_of_cluster_members = length(cluster_members) - 1
                if num_of_cluster_members ≥ 1 # should not even include clusters with a single member
                    cocluster_infectivity = 0.0
                    for k in cluster_members
                        if population.phase[k] == INFECTED || 
                            population.phase[k] == ASYMPTOMATIC
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
        phase = population.phase[n]
        if phase == SUSCEPTIBLE &&
                chances[n] ≤ 1 - exp(-λ[n] * population.susceptibility[n])
            population.event_history[n] = k
            population.phase[n] = EXPOSED
        elseif phase == EXPOSED && chances[n] ≤ γ.rate_expos
            population.event_history[n] = k
            if rand(rng) ≤ prob.asymp
                population.phase[n] = ASYMPTOMATIC
            else
                population.phase[n] = INFECTED
            end
        elseif phase == ASYMPTOMATIC && chances[n] ≤ γ.rate_asymp
            population.event_history[n] = k
            population.phase[n] = RECOVERED
        elseif phase == INFECTED && chances[n] ≤ γ.rate_infec
            population.event_history[n] = k
            if rand(rng) ≤ prob.decease
                population.phase[n] = DECEASED
            else
                population.phase[n] = RECOVERED
            end
        end
    end
    nothing
end

#' ## Time evolution

#' Main evolution function

function evolve!(rng, population, residences, clusters, τ, γ, prob, num_steps, time_step;
    verbose_step::Integer = 0
)
    num_population = length(population)
    evolution = spzeros(Phase, num_population, num_steps)
    evolution[:, 1]  .= population.phase # getfield(population, :phase)
    λ = Vector{Float64}(undef, num_population)
    chances = Vector{Float64}(undef, num_population)
    for k in 2:num_steps
        force_of_infection!(λ, population, residences, clusters, τ)
        step_foward!(rng, population, chances, λ, γ, prob, k)
        for n in 1:num_population
            phase = population.phase[n]
            if phase != SUSCEPTIBLE 
                evolution[n, k]  = phase
            end
        end
        if verbose_step > 0 && mod(k, verbose_step) == 0
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
