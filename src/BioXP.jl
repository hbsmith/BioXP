module BioXP

include("structs.jl")
include("randomize.jl")
include("parser.jl")
include("format.jl")
include("dg.jl")

import JSON
# using .Structs, .Randomize, .Parser

export 
# Functions
readexpansion,
readmaster,
readcompounds,
readids,
readkeyedids,
randomizecompounds,
expand,
find_minimal_seed_set,
simple_write_out,
formatbioxpoutput,
formatexpansion,
filter_reactions_by_dg,

# Structs
System,
Reaction,
Compound,
Expansion,

# Consts
Reactions,
Compounds,
IDs

## reaction_structs only give detailed info about reactions.
##   they are never supposed to be used as a comprehensive set of reactions
##   for the expansion process

function remove_rids_not_in_rstructs(
    rstructs::Reactions,
    rids::IDs)

    collect(String,intersect(rids,keys(rstructs)))
end

function matrixify_compounds(
    rstructs::Reactions,
    cids::IDs,
    rids::IDs)

    ## Reactons in columns, Compounds in rows. Rows, columns.
    R = zeros(Int, (length(cids),length(rids))) ## Reactants
    P = zeros(Int, (length(cids),length(rids))) ## Products

    for (i,r) in enumerate(rids)
        R[:,i] = [Int(cpd in rstructs[r].left) for cpd in cids]
        P[:,i] = [Int(cpd in rstructs[r].right) for cpd in cids]
    end

    R, P
end

function matrixify_targets(
    cids::IDs,
    tids::IDs)
    
    biosystem_tids = intersect(tids,cids) ## ensures impossible to reach tids don't count

    [Int(i in biosystem_tids) for i in cids]

end
"""
Return vector of compounds with 1 in position where seeds are present and 0 where seeds are absenst.
"""
function matrixify_seeds(
    sids::IDs, 
    cids::IDs)

    [Int(i in sids) for i in cids]
end

function identify_biosystem_compounds(
    rstructs::Reactions,
    rids::IDs)

    ## only look at rids which are in my rstructs
    
    IDs(unique(Iterators.flatten([vcat(rstructs[r].left,rstructs[r].right) for r in rids])))
end

"""
If we want to take into account dg values, dispatch on this function
"""
function expandmatrices(
    R::Array{Int,2}, 
    P::Array{Int,2}, 
    x::Array{Int,1},
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)

    # initialize variables:
    RT = transpose(R)
    PT = transpose(P)

    br = vec(sum(RT, dims=2)) # sum the rows of RT. Do I need the vec call here? yes, turns it into 1 element array
    bp = vec(sum(PT, dims=2))

    # force unallowed reactions to be out of reach
    if allowed_forward != nothing
        br = br + .~allowed_forward 
    end
    if allowed_backward != nothing
        bp = bp + .~allowed_backward  
    end
    
    # find the total number of metabolites in the seed set
    k = sum(x);

    # initialize previous iteration count of metabolites in the network
    k0 = 0;

    # iteration 1 consistes of the seed set
    X = Vector{Int}[x]

    # initialize reaction accumulation matrix
    Y = Vector{Int}[];

    # while the metabolite set has not converged
    while k > k0 
        # update previous number of metabolites
        k0 = k;

        # RT*x ==> represnts the number of present metabolites in the 
        # network within each reaction; if this isequal to the total 
        # number of metabolites in each reaction, then the reaction 
        # is added to the network

        ## Forward reactions
        y = RT * x .== br
        ## Backward reactions
        yp = PT * x .== bp

        # P*y > 0 ==> represents the vector of reactions that produce 
        # metabolite i. (i in 1:m).  If this is >0, 
        # then that metabolite is producable 

        ## Forward reactions
        xnew = P * y .> 0

        ## Backward reactions
        xnewp = R * yp .> 0

        # add to previous set of metabolites (only needed to retain seed set)

        ## Add forward and backward reactions
        x = (x .| xnew .| xnewp)
        y = Array{Int}(y .| yp)

        # find new total number of metabolites in network
        k = sum(x);

        # append accumulation matricies
        push!(X, x)
        push!(Y, y)
    end
    X, Y
end

function expand(
    system::System,
    write_path::Union{String,Nothing}=nothing)
    
    expand(system.rstructs,
        system.rids,
        system.sids,
        system.tids,
        write_path,
        system.allowed_forward,
        system.allowed_backward)
end

function expand(
    rstructs::Reactions,
    rids::IDs,
    sids::IDs,
    tids::IDs=IDs(),
    write_path::Union{String,Nothing}=nothing,
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)

    rids = remove_rids_not_in_rstructs(rstructs,rids)
    cids = identify_biosystem_compounds(rstructs,rids)

    (R, P) = matrixify_compounds(rstructs,cids,rids)
    x = matrixify_seeds(sids, cids)
    t = matrixify_targets(cids,tids)

    # X, Y = Vector{Int}[], Vector{Int}[] ## same as Vector{Vector{Int}}(),Vector{Vector{Int}}()
    X, Y = expandmatrices(R, P, x, allowed_forward, allowed_backward)

    if write_path !== nothing
        simple_write_out(write_path,x,t,cids,rids,X,Y)
    end

    x, t, cids, X, Y
end

"""
Return indices of seeds within the compounds vector.
"""
function seed_indicies(sids::IDs, cids::IDs)
    # This is a generator, not an array. You can iterate over this thing exactly once
    # because it only stores the current state and what it needs to find the next state.
    (findfirst(isequal(c), cids) for c in sids)
end


function find_minimal_seed_set(
    system::System,
    write_path::Union{String,Nothing}=nothing,
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)
    
    find_minimal_seed_set(system.rstructs,
        system.rids,
        system.sids,
        system.tids,
        write_path,
        allowed_forward,
        allowed_backward)
end

"""
    find_minimal_seed_set(rstructs,rids,sids,tids,write_path)

Return: 
- seeds (binary version) (these must be provided because their ordering matters for the algorithm)
- targets (binary version)
- compound ids
- reactions ids
- compounds (binary version) present at each timestep
- reactions (binary version) availabel at each timestep
"""
function find_minimal_seed_set(
    rstructs::Reactions,
    rids::IDs,
    sids::IDs,
    tids::IDs=IDs(),
    write_path::Union{String,Nothing}=nothing,
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)

    rids = remove_rids_not_in_rstructs(rstructs,rids)
    cids = identify_biosystem_compounds(rstructs,rids)
    println(cids)
    println(sids)

    (R, P) = matrixify_compounds(rstructs,cids,rids)
    t = matrixify_targets(cids,tids)
    # X, Y = Vector{Int}[], Vector{Int}[]
    x = matrixify_seeds(sids, cids) ## This should be a vector of all 1s of length(cids)
    println(x)
    println(typeof(x))
    println(ones(Int,length(cids)))
    println(typeof(ones(Int,length(cids))))
    x != ones(Int,length(cids)) && throw(DomainError("This should be a vector of all 1s of length(cids"))
    X, Y, x = loop_and_remove_seeds(sids,cids,x,t,R,P,allowed_forward,allowed_backward)

    if write_path !== nothing
        simple_write_out(write_path,x,t,cids,rids,X,Y)
    end

    x,t,cids,rids,X,Y

end

"""
    find_minimal_seed_set()

Return system variables after finding many minimal seed sets.
"""
function find_minimal_seed_set(
    rstructs::Reactions,
    rids::IDs,
    sid_sets::Vector{IDs},
    tids::IDs=IDs(),
    write_path::Union{String,Nothing}=nothing,
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)

    rids = remove_rids_not_in_rstructs(rstructs,rids)
    cids = identify_biosystem_compounds(rstructs,rids)
    (R, P) = matrixify_compounds(rstructs,cids,rids)
    t = matrixify_targets(cids,tids)
    # X, Y = Vector{Int}[], Vector{Int}[]
    all_seed_results = Vector{}
    if write_path==nothing ## no parallel processing because of the push
       
        for (i,sids) in enumerate(sid_sets)
            x = matrixify_seeds(sids, cids) ## This should be a vector of all 1s of length(cids)
            x !== ones(Int,length(cids)) && throw(DomainError("This should be a vector of all 1s of length(cids"))
            X, Y, x = loop_and_remove_seeds(sids,cids,x,t,R,P,allowed_forward,allowed_backward)
            push!(all_seed_results,(x,t,cids,rids,X,Y))
        end

    else ## Use parallel processing

        Threads.@threads for (i,sids) in collect(enumerate(sid_sets))
            x = matrixify_seeds(sids, cids) ## This should be a vector of all 1s of length(cids)
            x !== ones(Int,length(cids)) && throw(DomainError("This should be a vector of all 1s of length(cids"))
            X, Y, x = loop_and_remove_seeds(sids,cids,x,t,R,P,allowed_forward,allowed_backward)
            simple_write_out(joinpath(write_path,"$i.json"),x,t,cids,rids,X,Y)
        end
    
    end

    all_seed_results # This will be empty if write_path!=nothing

end

function loop_and_remove_seeds(
    sids::IDs,
    cids::IDs,
    x::Vector{Int},
    t::Vector{Int},
    R::Array{Int,2}, 
    P::Array{Int,2},
    allowed_forward::Union{Vector{Bool},Nothing}=nothing,
    allowed_backward::Union{Vector{Bool},Nothing}=nothing)

    tT = transpose(t)
    sum_t = sum(t)
    ## Run 1 network expansion per seed variation
    for i in seed_indicies(sids, cids)
        x[i] = 0

        global X, Y = expandmatrices(R, P, x, allowed_forward, allowed_backward) ## global needed to access the variables defined in-loop

        (tT * X[end]) != sum_t && (x[i] = 1) ## This is a short-circuit if statement
        # if (tT * X[end]) != sum_t
        #     x[i] = 1
        # end
    end
    
    X, Y, x

end

function simple_write_out(
    path::String, 
    x::Vector{Int}, 
    t::Vector{Int},
    cids::IDs, 
    rids::IDs,
    X::Vector{Vector{Int}}, 
    Y::Vector{Vector{Int}})

    # println("Writing out single network expansion...")
    # if !ispath(path)
    #     mkpath(path)
    # end

    data = Dict()
    data["x"] = x
    data["t"] = t
    data["compounds"] = cids
    data["reactions"] = rids
    data["X"] = X
    data["Y"] = Y

    open(path,"w") do f
        JSON.print(f, data, 2) #indent=2
    end
end

## Change the below to IOs?
# rstructs = readmaster("path/to/master.json")
# rids = readids("path/to/rids.json")
# sids = readids("path/to/sids.json")
# tids = readids("path/to/tids.json")
# path = "path/to/output/netexp/results.json"

# (x, t, cids, X, Y) = expand(rstructs,rids,sids,tids)

end