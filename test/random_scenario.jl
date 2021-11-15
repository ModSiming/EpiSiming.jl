@testset "Random Scenario" begin
    num_population = 10_000

    region_size = (6, 12)

    #' Distribution of residences per size
    res_size_distribution = [10, 22, 33, 22, 5, 5, 2, 1] |> v -> Weights(v / sum(v))

    #' Parameters of the Gamma distribution for the susceptibility and infectivity factors
    Γ = (
        sus_shape = 2.0,
        sus_scale = 0.5,
        inf_shape = 4.0,
        inf_scale= 0.25
    )

    #' Age pyramid
    pop_pyramid = let age_max = 100, p = 2, pyramid_func(a, age_max, p) = ( a + 1 ) * (age_max - a)^p
        Weights(
            pyramid_func.(0:age_max, age_max, p) / sum(pyramid_func.(0:age_max, age_max, p))
        )
    end

    #' Contact rates
    τ = Dict{Symbol, Float64}(
        :residences => 0.3,
        :work_places => 0.1,
        :school_places => 0.2
    )

    #' Recovery rates
    γ = (
        rate_expos = 0.25, # rate out of exposed
        rate_infec = 0.1, # rate out of infected
        rate_asymp = 0.2 # rate out of asymptomatic
    )

    #' Fate probabilities
    prob = (
        asymp = 0.6, # probability of becoming asymptomatic (vs. symptomatic = infected)
        decease = 0.02 # probability of deceasing
    )


    #' Set random number generator for repeatability in testings
    rng = MersenneTwister(123)

    #' generate population

    pop_blocks = EpiSiming.gen_random_pop_blocks(rng, num_population, region_size)

    @test size(pop_blocks) == region_size
    
    res_blocks_distrib = EpiSiming.gen_res_blocks(pop_blocks, res_size_distribution)

    residences = EpiSiming.gen_residences(rng, res_blocks_distrib)

    population = EpiSiming.gen_population(
        rng,
        residences,
        res_blocks_distrib,
        Γ,
        pop_pyramid
    )

    #' ### Initial infection

    num_exposed_at_time_0 = div(num_population, 500) # ( 1/500 = 0.002 = 0.2%) 20 for 10_000

    exposed_at_time_0 = sample(rng, 1:num_population, num_exposed_at_time_0, replace = false)
    
    for n in exposed_at_time_0
        population.phase[n] = EXPOSED
        next_phase, next_change = EpiSiming.transition_rules(rng, EXPOSED, 1)
        population.next_transition[n] = (next_phase, next_change)
    end

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

    #' ### Example

    num_steps = 360
    time_step = 1
    @time evolution = evolve!(
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
end