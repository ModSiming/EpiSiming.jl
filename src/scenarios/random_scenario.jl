#' # Random Scenario

#' ## Generating a random matrix of blocks with the population divided per block

"""
Given a number for the total population and a tuple with a rectangular area made of
square blocks, randomly distribute the population among the blocks.
"""
function gen_random_pop_blocks(rng, num_population, region_size)
    mean_per_block = div(num_population, prod(region_size))
    blocks = zeros(Int, region_size...)
    for k in eachindex(blocks)
        let s = sum(blocks)
            if s < num_population
                blocks[k] = min(num_population - s, rand(rng, 1:2mean_per_block))
            end
        end
    end
    let s = sum(blocks)
        if s < num_population
            k = rand(rng, eachindex(blocks))
            blocks[k] += num_population - s
        end
    end
    return blocks
end

gen_random_pop_blocks(num_population, region_size) = gen_random_pop_blocks(
    Random.default_rng(), num_population, region_size
)

#' ### Generating residences per block
#'
#' We are given a matrix `blocks` with the population in each block, and we want to
#' distribute them between residences according to a distribution of the number
#' of residences per size, `residence_size_distribution`.
#'
#' If $n$ is the (known) population of a block; $m$ is the (unknown) number of residences
#' in that block, and $p_j$ is the proportion of residences of size $j$, for
#' $j = 1, \ldots, J$, then there are $p_j m$ residences of size $j$ and $j p_j m$ residents
#' in a residence of size $j$. Thus,
#'
#' $$ n = p_1 m + 2 p_2 m + \ldots + J p_J m
#' $$
#'
#' Therefore, we find $m$ according to
#'
#' $$ m = \frac{n}{p_1 + 2p_2 + \ldots + Jp_J}.
#' $$
#'
#' But $p_1, \ldots, p_J$ are probabilities, so $n$ may not be divisible by the above
#' denominator and $m$ may not be an integer. Even if $m$ were an integer, some $m p_j$
#' may not be. So we need to adjust for round-off errors.
#'
#' So, we don't care to round $m$ (leave it as float) and take $m_j = floor(Int, m * p_j)$.
#' Then we distribute the remaining population...

function gen_res_blocks(rng, pop_blocks::Matrix, res_size_distribution::AbstractWeights)
    rsd = res_size_distribution |> v -> Weights(v / sum(v))
    res_blocks_distrib = fill(zeros(Int, length(rsd)), size(pop_blocks)...)
    for (k, n) in enumerate(pop_blocks)
        m = n / sum(rsd .* eachindex(rsd))
        res_blocks_distrib[k] = floor.(Int, m .* rsd)
        pop_block_remaining =
            n - sum(res_blocks_distrib[k] .* eachindex(rsd))
        while pop_block_remaining > 0
            res_max_size = min(pop_block_remaining, lastindex(rsd))
            temp_dist =
                rsd[begin:res_max_size] |> v -> Weights(v / sum(v))
            rs = sample(rng, eachindex(temp_dist), temp_dist)
            res_blocks_distrib[k][rs] += 1
            pop_block_remaining =
                n - sum(res_blocks_distrib[k] .* eachindex(rsd))
        end
    end
    @info ifelse(
        sum(
            abs,
            map(
                r -> sum(r .* eachindex(rsd)),
                res_blocks_distrib) - pop_blocks
            ) == 0,
        "population successfully distributed within blocks",
        "discrepancy in the distribution of the population within blocks"
    )
    return res_blocks_distrib
end

gen_res_blocks(pop_blocks::Matrix, res_size_distribution::AbstractWeights) =
gen_res_blocks(Random.default_rng(), pop_blocks, res_size_distribution)

#' ### Generating the list of residences with all the required information

#block_grid_num(r) = r # spread out
#block_grid_num(r) = isqrt(r) + 1 # this is as tightly as it can get
block_grid_num(r) = 2 * isqrt(r) # this is if we don't want too tight
#block_grid_num(r) = div(r, 2) + 1 # even less tight
# If `√r = s + ε`, with an integer `s` and `0 ≤ ε < 1`, then `r = (a + ε)^2 < (a + 1)^2`.

function gen_residences(rng, res_blocks_distrib)

    region_size = size(res_blocks_distrib)

    num_residences = sum(Iterators.flatten(res_blocks_distrib))
    residences = Residences(
        Vector{Int}(undef, num_residences),
        Vector{Tuple{Float64, Float64}}(undef, num_residences),
        Vector{Int}(undef, num_residences),
        Vector{Vector{Int}}(undef, num_residences)
    )

    res_index = 0
    pop_index = 0
    @showprogress "Generating $num_residences residences... " for j in 1:region_size[2],
            i in 1:region_size[1]
        num_res_block = sum(res_blocks_distrib[i, j])
        bgn = block_grid_num(num_res_block) # grid size
        pos = sample(
            rng,
            collect(
                Iterators.product(i .+ (0:bgn - 1) ./ bgn, j .+ (0:bgn - 1) ./ bgn)
            ),
            num_res_block,
            replace = false
        )
        block_res_index = 0
        for (k, m) in enumerate(res_blocks_distrib[i, j])
            for r in 1:m
                residences.block[res_index + r] = i + (j - 1) * region_size[1]
                residences.position[res_index + r] = pos[block_res_index + r]
                residences.num_residents[res_index + r] = k
                residences.residents[res_index + r] = collect(pop_index + 1 : pop_index + k)
                pop_index += k
            end
            res_index += m
            block_res_index += m
        end
    end

    return residences
end

gen_residences(res_blocks_distrib) = gen_residences(
    Random.default_rng(),
    res_blocks_distrib
)

#' ### Generating population

#' ### Generating agents
#'

# infectivity and susceptibility (mean = shape * scale)
# didn't see anything about correlation between the two, so I am assuming none

function gen_population(
    rng,
    residences,
    res_blocks_distrib,
    Γ,
    pop_pyramid
)

    num_population = residences.residents[end][end]

    pop_susceptibility = rand(
        rng,
        Distributions.Gamma(Γ.sus_shape, Γ.sus_scale),
        num_population
    )

    pop_infectivity = rand(
        rng,
        Distributions.Gamma(Γ.inf_shape, Γ.inf_scale),
        num_population
    )

    ages = Int8.(sample(rng, 0 : length(pop_pyramid) - 1, pop_pyramid, num_population))

    population = Population(
        fill(SUSCEPTIBLE, num_population), # phase
        fill(1, num_population), # event_history
        fill((EXPOSED, typemax(Int)), num_population), # transition
        Vector{Int}(undef, num_population), # residences
        Vector{Tuple{Float64, Float64}}(undef, num_population), # positions
        ages, # ages
        pop_susceptibility, # susceptibility
        pop_infectivity, # infectivity
        Vector{Dict{Symbol, Int}}(undef, num_population), # clusters
        Vector{Dict{Symbol, Int}}(undef, num_population) # networks
    )

    @showprogress "Generating $num_population agents... " for m in eachindex(residences)
        radius = 0.25 / block_grid_num(sum(res_blocks_distrib[residences.block[m]]))
        res_position = residences.position[m]
        for (k, n) in enumerate(residences.residents[m])
            θ = 2π * (k - 1) / residences.num_residents[m]
            population.position[n] = res_position .+ radius .* (cos(θ), sin(θ))
            population.residence[n] = m
            population.clusters[n] = Dict{Symbol, Int}()
            population.networks[n] = Dict{Symbol, Int}()
        end
    end
    return population
end

function gen_population(
    residences,
    res_blocks_distrib,
    Γ,
    pop_pyramid
)
    return gen_population(
        Random.default_rng(),
        residences,
        res_blocks_distrib,
        Γ,
        pop_pyramid
    )
end


#' ### Generate clusters

decay(r, α) = 1 / (1 + r^α)

function gen_clusters(rng, ind_pop_available, max_size, α)
    ind_pop_shuffled = shuffle(ind_pop_available)
    clusters = Vector{Vector{Int}}()
    assigned = 0
    while assigned < length(ind_pop_available)
        cluster_size = sample(rng, 1:max_size, Weights(decay.(1:max_size, α)))
        push!(clusters, ind_pop_shuffled[assigned + 1 : min(end, assigned + cluster_size)])
        assigned += cluster_size
    end
    return clusters
end

gen_clusters(ind_pop_available, max_size, α) =
    gen_clusters(Random.default_rng(), ind_pop_available, max_size, α)

function transition_rules(rng, phase::Phase, k)
    if phase == EXPOSED
        next_phase = ifelse(rand(rng) < 0.35, ASYMPTOMATIC, INFECTED)
        if next_phase == ASYMPTOMATIC
#=             next_change = k + sample(
                1:5,
                Weights([0.1, 0.3, 0.2, 0.1, 0.1])
            ) =#
            next_change = k + sample(
                1:8,
                Weights([0.04, 0.08, 0.16, 0.31, 0.28, 0.08, 0.04, 0.01])
            )
        else
            next_change = k + sample(
                1:8,
                Weights([0.03, 0.07, 0.14, 0.28, 0.30, 0.10, 0.06, 0.02])
            )
        end
    elseif phase == ASYMPTOMATIC
        next_phase = RECOVERED
        next_change = k + sample(
                1:16,
                Weights([0.01, 0.02, 0.03, 0.04, 0.08, 0.06, 0.15, 0.16, 0.15, 0.14, 0.06, 0.04, 0.02, 0.02, 0.01, 0.01])
            )
    elseif phase == INFECTED
        next_phase = ifelse(rand(rng) < 0.97, RECOVERED, DECEASED)
        if next_phase == RECOVERED
            next_change = k + sample(
                1:16,
                Weights([0.01, 0.02, 0.03, 0.04, 0.08, 0.06, 0.15, 0.16, 0.15, 0.14, 0.06, 0.04, 0.02, 0.02, 0.01, 0.01])
            )
        else
            next_change = k + sample(
                1:16,
                Weights([0.01, 0.02, 0.03, 0.04, 0.08, 0.06, 0.15, 0.16, 0.15, 0.14, 0.06, 0.04, 0.02, 0.02, 0.01, 0.01])
            )
        end
    elseif phase == RECOVERED
        next_phase = RECOVERED
        next_change = typemax(k) # should be T from Population{T, S}
    elseif phase == DECEASED
        next_phase = DECEASED
        next_change = typemax(k) # should be T from Population{T, S}
    else
        throw(
            ArgumentError(
                "Evolution rule not implemented for phase $phase"
            )
        )
    end
    return next_change, next_phase
end

transition_rules(phase::Phase, k) = transition_rules(Random.default_rng(), phase::Phase, k)
