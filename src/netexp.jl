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

# function transform_
function netexp(
    R::Array{Int,2}, 
    P::Array{Int,2}, 
    RT::TArray{Int,2}, 
    PT::TArray{Int,2},
    b::Vector{Int}, 
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
        y = RT * x .== b
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
    t = matrixify_targets(biosystem_compounds,target_compounds)
    x = matrixify_seeds(seed_compounds, biosystem_compounds)

    RT = transpose(R)
    PT = transpose(P)

    br = vec(sum(RT, dims=2)) # sum the rows of RT. Do I need the vec call here?
    bp = vec(sum(PT, dims=2))

    tT = transpose(t)
    sum_t = sum(t)

    

    X, Y = Vector{Int}[], Vector{Int}[]

    X, Y = netexp(R, P, RT, PT, br, bp, x)

    println("Writing out single network expansion...")
    if !ispath(path)
        mkpath(path)
    end
    fullpath = path*"reaction_edges_P.json"
    simple_write_out(fullpath, x, t, compounds, reactions, X, Y)

