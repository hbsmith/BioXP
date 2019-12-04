module Netexp

using structs

## reaction_structs only give detailed info about reactions.
##   they are never supposed to be used as a comprehensive set of reactions
##   for the expansion process

function matrixify_compounds(
    reaction_structs::Reactions,
    biosystem_compounds::IDs,
    biosystem_reactions::IDs,
    target_compounds::IDs)

    ## For R and P (equivilent to R and P in segrelab's github):
    ## Reactons in columns, Compounds in rows. Rows, columns.

    R = zeros(Int, (length(biosystem_compounds),length(biosystem_reactions)))
    P = zeros(Int, (length(biosystem_compounds),length(biosystem_reactions)))

    for (i,r) in enumerate(biosystem_reactions)
        R[:,i] = [Int(cpd in reaction_structs[r].left) for cpd in biosystem_compounds]
        P[:,i] = [Int(cpd in reaction_structs[r].right) for cpd in biosystem_compounds]
    end

    R, P
end

function matrixify_targets(
    biosystem_compounds::IDs,
    target_compounds::IDs)
    
    biosystem_target_compounds = intersect(target_compounds,biosystem_compounds)
    t = [Int(i in biosystem_target_compounds) for i in biosystem_compounds]
    
    t
end

function matrixify_seeds(
    seed_compounds::IDs, 
    biosystem_compounds::IDs)

    [Int(i in seed_compounds) for i in biosystem_compounds]
end

function netexp(
    R::Array{Int,2}, 
    P::Array{Int,2}, 
    RT::TArray{Int,2}, 
    PT::TArray{Int,2},
    br::Vector{Int}, 
    bp::Vector{Int}, 
    x::Array{Int,1})

    # initialize variables:
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
    reaction_structs::Reactions,
    biosystem_reactions::IDs,
    seed_compounds::IDs,
    target_compounds::IDs=IDs())

    biosystem_compounds = IDs(unique(Iterators.flatten([vcat(reaction_structs[r].left,reaction_structs[r].right) for r in biosystem_reactions])))

    (R, P) = matrixify_compounds(reaction_structs,biosystem_compounds,biosystem_reactions,target_compounds)
    x = matrixify_seeds(seed_compounds, biosystem_compounds)
    t = matrixify_targets(biosystem_compounds,target_compounds)

    ## Should I just move all the below into netexp itself?
    RT = transpose(R)
    PT = transpose(P)

    br = vec(sum(RT, dims=2)) # sum the rows of RT. Do I need the vec call here?
    bp = vec(sum(PT, dims=2))    

    X, Y = Vector{Int}[], Vector{Int}[]
    X, Y = netexp(R, P, RT, PT, br, bp, x)

    x, t, biosystem_compounds, X, Y
end

function simple_write_out(
    path::String, 
    x::Array{Int,1}, 
    t::Array{Int,1},
    biosystem_compounds::Array{String,1}, 
    biosystem_reactions::Array{String,1},
    X::Vector{Vector{Int}}, 
    Y::Vector{Vector{Int}})

    println("Writing out single network expansion...")
    if !ispath(path)
        mkpath(path)
    end

    data = Dict()
    data["x"] = x
    data["t"] = t
    data["compounds"] = biosystem_compounds
    data["reactions"] = biosystem_reactions
    data["X"] = X
    data["Y"] = Y

    open(path,"w") do f
        JSON.print(f, data, 2)
    end
end

reaction_structs = readmaster("path/to/master.json")
biosystem_reactions = readids("path/to/biosystem_reactions.json")
seed_compounds = readids("path/to/seed_compounds.json")
target_compounds = readids("path/to/target_compounds.json")

(x, t, biosystem_compounds, X, Y) = expand(reaction_structs,biosystem_reactions,seed_compounds,target_compounds)
simple_write_out(path, x, t, biosystem_compounds, biosystem_reactions, X, Y)
