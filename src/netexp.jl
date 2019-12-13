module Netexp

using Structs
using Randomize

## reaction_structs only give detailed info about reactions.
##   they are never supposed to be used as a comprehensive set of reactions
##   for the expansion process

function matrixify_compounds(
    rstructs::Reactions,
    cids::IDs,
    rids::IDs)

    ## For R and P (equivilent to R and P in segrelab's github):
    ## Reactons in columns, Compounds in rows. Rows, columns.

    R = zeros(Int, (length(cids),length(rids)))
    P = zeros(Int, (length(cids),length(rids)))

    for (i,r) in enumerate(rids)
        R[:,i] = [Int(cpd in rstructs[r].left) for cpd in cids]
        P[:,i] = [Int(cpd in rstructs[r].right) for cpd in cids]
    end

    R, P
end

function matrixify_targets(
    cids::IDs,
    tids::IDs)
    
    biosystem_tids = intersect(tids,cids)
    t = [Int(i in biosystem_tids) for i in cids]
    
    t
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

    IDs(unique(Iterators.flatten([vcat(rstructs[r].left,rstructs[r].right) for r in rids])))
end

function expandmatrices(
    R::Array{Int,2}, 
    P::Array{Int,2}, 
    x::Array{Int,1})

    # initialize variables:
    RT = transpose(R)
    PT = transpose(P)

    br = vec(sum(RT, dims=2)) # sum the rows of RT. Do I need the vec call here? yes, turns it into 1 element array
    bp = vec(sum(PT, dims=2))   
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
    write_path::{String,nothing}=nothing)
    
    expand(system.rstructs,
        system.rids,
        system.sids,
        system.tids,
        write_path)
end

function expand(
    rstructs::Reactions,
    rids::IDs,
    sids::IDs,
    tids::IDs=IDs(),
    write_path::{String,nothing}=nothing)

    cids = identify_biosystem_compounds(rstructs,rids)

    (R, P) = matrixify_compounds(rstructs,cids,rids,tids)
    x = matrixify_seeds(sids, cids)
    t = matrixify_targets(cids,tids)

    # X, Y = Vector{Int}[], Vector{Int}[] ## same as Vector{Vector{Int}}(),Vector{Vector{Int}}()
    X, Y = expandmatrices(R, P, x)

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
    write_path::{String,nothing}=nothing)
    
    find_minimal_seed_set(system.rstructs,
        system.rids,
        system.sids,
        system.tids,
        write_path)
end

function find_minimal_seed_set(
    rstructs::Reactions,
    rids::IDs,
    sids::IDs,
    tids::IDs=IDs(),
    write_path::{String,nothing}=nothing)

    cids = identify_biosystem_compounds(rstructs,rids)
    

    (R, P) = matrixify_compounds(rstructs,cids,rids,tids)
    t = matrixify_targets(cids,tids)
    # X, Y = Vector{Int}[], Vector{Int}[]
    x = matrixify_seeds(sids, cids) ## This should be a vector of all 1s of length(cids)
    x !== ones(Int,length(cids)) && throw(DomainError("This should be a vector of all 1s of length(cids"))
    X, Y, x = loop_and_remove_seeds(sids,cids,x,t,R,P)

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
    write_path::{String,nothing}=nothing)

    cids = identify_biosystem_compounds(rstructs,rids)
    (R, P) = matrixify_compounds(rstructs,cids,rids,tids)
    t = matrixify_targets(cids,tids)
    # X, Y = Vector{Int}[], Vector{Int}[]
    all_seed_results = Vector{}
    for (i,sids) in enumerate(sid_sets)
        x = matrixify_seeds(sids, cids) ## This should be a vector of all 1s of length(cids)
        x !== ones(Int,length(cids)) && throw(DomainError("This should be a vector of all 1s of length(cids"))
        X, Y, x = loop_and_remove_seeds(sids,cids,x,t,R,P)

        if write_path !== nothing
            simple_write_out(joinpath(write_path,"$i.json"),x,t,cids,rids,X,Y)
        end

        push!(all_seed_results,(x,t,cids,rids,X,Y))
    
    end

    all_seed_results # Is this going to use a massive amount of memory?

end

function loop_and_remove_seeds(
    sids::IDs,
    cids::IDs,
    x::Vector{Int},
    t::Vector{Int},
    R::Array{Int,2}, 
    P::Array{Int,2})

    tT = transpose(t)
    sum_t = sum(t)
    ## Run 1 network expansion per seed variation
    for i in seed_indicies(sids, cids)
        x[i] = 0

        X, Y = expandmatrices(R, P, x)

        (tT * X[end]) != sum_t && x[i] = 1 ## This is a short-circuit if statement
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

    println("Writing out single network expansion...")
    if !ispath(path)
        mkpath(path)
    end

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
rstructs = readmaster("path/to/master.json")
rids = readids("path/to/rids.json")
sids = readids("path/to/sids.json")
tids = readids("path/to/tids.json")
path = "path/to/output/netexp/results.json"

(x, t, cids, X, Y) = expand(rstructs,rids,sids,tids)
