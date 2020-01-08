using Random

# export randomizecompounds

"""

    sort_biosystem_compounds(compound_structs,
        biosystem_compounds,
        sortkey,
        zero_mass_behavior)

Sorts compounds by increasing mass.
Always shuffles the order of the zero-mass compounds.

...
# Arguments
- `compound_structs::Compounds`: compounds with metadata
- `biosystem_compounds::IDs`: compounds to sort
- `sortkey::Symbol=:exact_mass`: how to sort compounds. Can also sort by `:molecular_weight`
- `zero_mass_behavior::String="end"`: what to do with compounds that have no mass. Can also be `"random"`
...
"""
function sort_biosystem_compounds(
    compound_structs::Compounds,
    biosystem_compounds::IDs,
    sortkey::Symbol=:exact_mass,
    zero_mass_behavior::String="end")

    cpd_masses = [(i,compound_structs[i][sortkey]) for i in biosystem_compounds if compound_structs[i][sortkey]!=0]
    cpd_zero_masses = [(i,compound_structs[i][sortkey]) for i in biosystem_compounds if compound_structs[i][sortkey]==0]

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

function mix_it_up!(
    tuples::Vector{Tuple{String,Float64}},
    beta::Float64,
    n_swaps::Int)

    for _ in 1:n_swaps
        i,j = rand(1:length(tuples),2) ## 2 random indices 
        if swap_random(tuples[i],tuples[j],beta) == true
            tuples[i], tuples[j] = tuples[j], tuples[i] ## changes them simultaneously
        end 
    end 
end

function swap_random(
    ti::Tuple{String,Float64},
    tj::Tuple{String,Float64},
    beta::Float64)
    
    diff = ti[1] - tj[1]
    
    ## Determine probability of flipping (another way of writing below)
    # ti[1]==0 || tj[1]==0 ? p=.5 : diff>0 ? p = math.e**(-diff/float(beta)) : p=1.0

    if ti[1]==0 || tj[1]==0
        p = .5
    elseif diff > 0
        p = exp(-diff/beta)
    else
        p = 1.0
    end

    rand(Float64)<p ? true : false
end

function randomizecompounds(
    biosystem_compounds::IDs,
    compound_structs::Compounds,
    n_runs::Int,
    n_swaps::Int=1000,
    beta::Float64=20,
    sortkey::Symbol=:exact_mass,
    zero_mass_behavior::String="end"
    )

    randomized_cpd_lists = [Vector{String}(undef,length(biosystem_compounds)) for _ in 1:n_runs] # allocate output memory
    for r in 1:n_runs
        cpds_masses = sort_biosystem_compounds(compound_structs, biosystem_compounds, sortkey, zero_mass_behavior)
        mix_it_up!(cpds_masses,beta,n_swaps)
        randomized_cpd_lists[r] = [c[1] for c in cpds_masses]
    end
    return randomized_cpd_lists
end