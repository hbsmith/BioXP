## Network expansion using lists/dicts
import JSON
using LinearAlgebra

const EdgesType = Dict{String,Dict{String,Vector{String}}}
const TargetsType = Dict{String,String}

function prepare_matrices_and_targets(reaction_edges_json::String, target_json::String)
    reaction_edges = EdgesType(JSON.parsefile(reaction_edges_json))

    ## could use ["substrates""] too
    reactions = [i for i in keys(reaction_edges["products"])]

    cpd_reactants_flat = unique(Iterators.flatten(values(reaction_edges["substrates"])))
    cpd_products_flat = unique(Iterators.flatten(values(reaction_edges["products"])))
    compounds = unique(Iterators.flatten([cpd_reactants_flat,cpd_products_flat]))

    ## For R and P:
    ## Reactons in columns
    ## Compounds in rows
    ## Rows, Columns
    ## This is equivilent to R and P in segrelab's github
    R = zeros(Int, (length(compounds),length(reactions)))
    P = zeros(Int, (length(compounds),length(reactions)))

    for (i,r) in enumerate(reactions)
        R[:,i] = [Int(cpd in reaction_edges["substrates"][r]) for cpd in compounds]
        P[:,i] = [Int(cpd in reaction_edges["products"][r]) for cpd in compounds]
    end
    ################################################################################
    target_compounds = TargetsType(JSON.parsefile(target_json));
    target_compounds_org = collect(intersect(keys(target_compounds),compounds))
    t = [Int(i in target_compounds_org) for i in compounds]
    ################################################################################

    R, P, compounds, reactions, t
end

function prepare_seeds(seed_list::Vector{String}, compounds::Vector{String})
    [Int(i in seed_list) for i in compounds]
end

## Deprecated
# function prepare_seeds(seed_json::String,seed_key::String,compounds::Vector{String})	
#     seed_dict = JSON.parsefile(seed_json)	
#     [Int(i in seed_dict[seed_key]) for i in compounds]	
# end

function seed_indicies(seed_list::Vector{String}, compounds::Vector{String})
    # This is a generator, not an array. You can iterate over this thing exactly once
    # because it only stores the current state and what it needs to find the next state.
    (findfirst(isequal(c), compounds) for c in seed_list)
end

# This is just cosmetic
const TArray{T, N} = Transpose{T, Array{T, N}}

function netexp(R::Array{Int,2}, P::Array{Int,2}, RT::TArray{Int,2}, PT::TArray{Int,2},
                b::Vector{Int}, bp::Vector{Int}, x::Array{Int,1})
    
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

function simple_write_out(outpath::String, x::Array{Int,1}, t::Array{Int,1},
                          compounds::Array{String,1}, reactions::Array{String,1},
                          X::Vector{Vector{Int}}, Y::Vector{Vector{Int}})
    data = Dict()
    data["x"] = x
    data["t"] = t
    data["compounds"] = compounds
    data["reactions"] = reactions
    data["X"] = X
    data["Y"] = Y
    
    open(outpath,"w") do f
        JSON.print(f, data)
    end
end

function enumerate_minimal_seed_sets(TARGETJSON::String,EDGEDIR::String,SEEDDIR::String,OUTDIR::String)

    for FNAME in readdir(EDGEDIR) 

        FULLEDGEPATH = joinpath(EDGEDIR, FNAME) #json of all edges for organism
        FULLSEEDPATH = joinpath(SEEDDIR, FNAME) #json of all seeds for organism
        
        if isfile(FULLEDGEPATH)
            
            if last(splitext(FULLEDGEPATH)) == ".json"
                println("Finding minimal seeds for: $FNAME")
                
                R, P, compounds, reactions, t = prepare_matrices_and_targets(FULLEDGEPATH, TARGETJSON)

                RT = transpose(R)
                PT = transpose(P)

                b = vec(sum(RT, dims=2))
                bp = vec(sum(PT, dims=2))

                tT = transpose(t)
                sum_t = sum(t)
                
                all_of_the_seeds = Vector{Vector{String}}(JSON.parsefile(FULLSEEDPATH))
                for (n_seed, seed_list) in enumerate(all_of_the_seeds)
                    OUTDIRWITHORGNAME = joinpath(OUTDIR, first(splitext(FNAME)))

                    if !ispath(OUTDIRWITHORGNAME)
                        mkpath(OUTDIRWITHORGNAME)
                    end

                    # I want 1 randomizaiton per outpath
                    FULLOUTPATH = joinpath(OUTDIRWITHORGNAME, "$n_seed.json")

                    x = prepare_seeds(seed_list, compounds)

                    X, Y = Vector{Int}[], Vector{Int}[]
                    for i in seed_indicies(seed_list, compounds)
                        x[i] = 0

                        X, Y = netexp(R, P, RT, PT, b, bp, x)

                        if (tT * X[end]) != sum_t
                            x[i] = 1
                        end
                    end

                    println("Writing out randomization: $n_seed")
                    simple_write_out(FULLOUTPATH, x, t, compounds, reactions, X, Y)
                end
            end
        end
    end
end

# seed_compounds = JSON.parsefile("../seeds.json");
# x = [Int(i in seed_compounds["Enceladus_20-SAFR-032"]) for i in compounds];

# ds80_seeds = ["C00011",
# "C20298",
# "C14819",
# "C00087",
# "C00237",
# "C00058",
# "C00033",
# "C00031",
# "C00095",
# "C00124",
# "C00159",
# "C00243",
# "C00208",
# "C00282",
# "C00007",
# "C00001"]
########################################
#### CHECK MINIMAL SEED SET ######
########################################
# # seedkey = "Enceladus_20-SAFR-032"

# # SEEDJSON = "29012812801.json" ## Contains keys numbered 1-100, with values of random compounds

# const TARGETJSON = "targets/Freilich09.json"
# const EDGEDIR = "jgi/2018-09-29/ph_edge_jsons/archaea_split/a00s/"
# const SEEDDIR = "seeds/minimal_seed_randomizations/archaea/"
# const OUTDIR = "results/simple/minimal_seed_randomizations_fixed/archaea/a00s/"  #*split(SEEDJSON,".json")[1]*"/"

# if !ispath(OUTDIR)
#     mkpath(OUTDIR)
# end

# enumerate_minimal_seed_sets(TARGETJSON,EDGEDIR,SEEDDIR,OUTDIR)

########################################
#### MANY NETWORK EXPANSION RUN ######
########################################
seedkey = "Enceladus_20-SAFR-032_P"

SEEDJSON = "seeds.json"
TARGETJSON = "targets/Freilich09.json"
DATADIR = "jgi/2019-09-09/ph_edge_jsons_P/bacteria/"

# fsplit = split(DATADIR,"/")
# OUTDIR = "results/simple/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
OUTDIR = "results/simple_2019-09-09/ph_edge_jsons_P/bacteria/"

if ispath(OUTDIR)==false
    mkpath(OUTDIR)
end

## Prepare seed_list
seed_dict = JSON.parsefile(SEEDJSON,dicttype=Dict{String,Vector{String}})
seed_list = seed_dict[seedkey]

for FNAME in readdir(DATADIR)
    FULLINPATH = DATADIR*FNAME
    FULLOUTPATH = OUTDIR*FNAME
    if isfile(FULLINPATH)==true
        if split(FULLINPATH,".")[2] == "json"
            ## DO MAIN
            (R,P,compounds,reactions,t) = prepare_matrices_and_targets(FULLINPATH,TARGETJSON)

            RT = transpose(R)
            PT = transpose(P)

            b = vec(sum(RT, dims=2))
            bp = vec(sum(PT, dims=2))

            tT = transpose(t)
            sum_t = sum(t)

            x = prepare_seeds(seed_list, compounds)
            # (X,Y) = netexp(R,P,x)
            X, Y = Vector{Int}[], Vector{Int}[]
            X, Y = netexp(R, P, RT, PT, b, bp, x)
            simple_write_out(FULLOUTPATH,x,t,compounds,reactions,X,Y)
        end
    end
end

########################################
#### SINGLE NETWORK EXPANSION RUN ######
########################################
# ## Inputs
# # ds80_seeds = ["C00031","C00001"]
# seed_list = ["C00001","C00011","C00237","C00282","C00067","C00132","C06548","C00469","C00283","C00014","C00697","C01326","C01438","C01548","C06547","C11505","C20783","C01407"]#,"C00009"]
# reaction_edges_json = "kegg/2018-09-25/reaction_edges.json"
# target_json = "targets/Freilich09.json"
# # seed_json = "seeds.json"

# ## Create out path
# # fsplit = split(reaction_edges_json,"/")
# # path = "results/netexpdata_jsons/"*fsplit[end-2]*"/"*fsplit[end-1]
# # path = "results/data_glucose_test/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
# path ="results/simple_2019-09-09/kegg_edge_json/"

# if !ispath(path)
#     mkpath(path)
# end
# fullpath = path*"reaction_edges_P.json"

# ## DO MAIN
# (R,P,compounds,reactions,t) = prepare_matrices_and_targets(reaction_edges_json,target_json)

# RT = transpose(R)
# PT = transpose(P)

# b = vec(sum(RT, dims=2))
# bp = vec(sum(PT, dims=2))

# tT = transpose(t)
# sum_t = sum(t)

# x = prepare_seeds(seed_list, compounds)

# X, Y = Vector{Int}[], Vector{Int}[]
# # for i in seed_indicies(seed_list, compounds)
# #     x[i] = 0

# X, Y = netexp(R, P, RT, PT, b, bp, x)

#     # if (tT * X[end]) != sum_t
#     #     x[i] = 1
#     # end
# # end

# println("Writing out single network expansion...")
# simple_write_out(fullpath, x, t, compounds, reactions, X, Y)
