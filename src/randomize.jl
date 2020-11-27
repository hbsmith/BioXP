using Random

"""

    sort_biosystem_compounds(compound_structs,
        biosystem_compounds,
        sortkey,
        zero_mass_behavior)

Return a vector of (compound,mass) sorted by increasing mass.
Always shuffles the order of the zero-mass compounds.

...
# Arguments
- `compound_structs::Compounds`: compounds with metadata
- `biosystem_compounds::IDs`: compounds to sort
- `sortkey::Symbol=:exact_mass`: how to sort compounds. Can also sort by `:mol_weight`
- `zero_mass_behavior::String="end"`: what to do with compounds that have no mass. Can also be `"random"`
...
"""
function sort_biosystem_compounds(
    compound_structs::Compounds,
    biosystem_compounds::IDs,
    sortkey::Symbol=:exact_mass,
    zero_mass_behavior::String="end")

    cpd_masses = [(i,getproperty(compound_structs[i],sortkey)) for i in biosystem_compounds if getproperty(compound_structs[i],sortkey)!=0]
    cpd_zero_masses = [(i,getproperty(compound_structs[i],sortkey)) for i in biosystem_compounds if getproperty(compound_structs[i],sortkey)==0]

    sort!(cpd_masses, by= x->x[2])

    if zero_mass_behavior=="end"
        return vcat(cpd_masses,shuffle(cpd_zero_masses))
    
    elseif zero_mass_behavior=="random"
        for tup in cpd_zero_masses
            insert!(cpd_masses,rand(1:length(cpd_masses)),tup) # this will permit insertion at longer indices as cpd_masses grows
        end
        return cpd_masses
    end
end

"""
    mix_it_up(tuples,beta,n_swaps)

Swap random tuples within `tuples`.
Based on formula from Handorf et al. 2008.
"""
function mix_it_up!(
    tuples::Vector{<:Tuple{String,Real}},
    beta::Real,
    n_swaps::Int,
    rng::AbstractRNG=Random.default_rng()) #where {R <: Real}

    for _ in 1:n_swaps
        i,j = rand(1:length(tuples),2) ## 2 random indices 
        if swap_random(tuples[i],tuples[j],beta,rng) == true
            tuples[i], tuples[j] = tuples[j], tuples[i] ## changes them simultaneously
        end 
    end 
end

"""
    swap_random(ti,tj,beta)

Return `true` if tuples should be flipped.
"""
function swap_random(
    ti::Tuple{String,<:Real},
    tj::Tuple{String,<:Real},
    beta::Real,
    rng::AbstractRNG=Random.default_rng())
    
    p = getflipprobability(ti,tj,beta)

    rand(rng,Float64)<p ? true : false
end

"""
    getflipprobability(ti,tj,beta)

Get the probability of flipping tuples ti and tj.
Based on formula from Handorf et al. 2008.
"""
function getflipprobability(
    ti::Tuple{String,Real},
    tj::Tuple{String,Real},
    beta::Real)

    diff = ti[2] - tj[2]
    
    ## Determine probability of flipping (another way of writing below)
    # ti[1]==0 || tj[1]==0 ? p=.5 : diff>0 ? p = math.e**(-diff/float(beta)) : p=1.0

    if ti[2]==0 || tj[2]==0
        p = .5
    elseif diff > 0
        p = exp(-diff/beta)
    else
        p = 1.0
    end

    return p
    
end

"""
    randomizecompounds(biosystem_compounds,
        compound_structs,
        n_runs,
        n_swaps,
        beta,
        sortkey,
        zero_mass_behavior)

Return a vector of run results. 
Each run result is a randomized list of biosystem_compounds.

...
# Arguments
- `biosystem_compounds::IDs`: compounds to sort
- `compound_structs::Compounds`: compounds with metadata
- `n_runs::Int`: number of runs
- `n_swaps::Int=1000`: number of compounds to within biosystem_compounds to swap when randomizing
- `beta::Real=20`: Mixing coefficient. See Handorf et al. 2008 for details.
- `sortkey::Symbol=:exact_mass`: how to sort compounds. Can also sort by `:mol_weight`
- `zero_mass_behavior::String="end"`: what to do with compounds that have no mass. Can also be `"random"`
...
"""
function randomizecompounds(
    biosystem_compounds::IDs,
    compound_structs::Compounds,
    n_runs::Int,
    n_swaps::Int=1000,
    beta::Real=20,
    sortkey::Symbol=:exact_mass,
    zero_mass_behavior::String="end",
    rng::AbstractRNG=Random.default_rng())

    randomized_cpd_lists = [Vector{String}(undef,length(biosystem_compounds)) for _ in 1:n_runs] # allocate output memory
    ## Generate seeds and RNG objects based on rng
    rng_seeds = [rand(rng,1:typemax(Int)) for i in 1:n_runs]
    rng_list = [MersenneTwister(i) for i in rng_seeds]
    for r in 1:n_runs
        cpds_masses = sort_biosystem_compounds(compound_structs, biosystem_compounds, sortkey, zero_mass_behavior)
        mix_it_up!(cpds_masses,beta,n_swaps,rng_list[r])
        randomized_cpd_lists[r] = [c[1] for c in cpds_masses]
    end
    return randomized_cpd_lists
end