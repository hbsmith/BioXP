## Network expansion using lists/dicts
import JSON
using Profile
using LinearAlgebra

function prepare_matrices_and_targets(reaction_edges_json::String,target_json::String)

    reaction_edges = JSON.parsefile(reaction_edges_json)

    reactions = [i for i in keys(reaction_edges["products"])] ## could use ["substrates""] too

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
    target_compounds = JSON.parsefile(target_json);
    target_compounds_org = collect(intersect(keys(target_compounds),compounds))
    t = [Int(i in target_compounds_org) for i in compounds]
    ################################################################################
    (R,P,compounds,reactions,t)

end

function prepare_seeds(seed_list::Array{String,1},compounds::Array{String,1})
    [Int(i in seed_list) for i in compounds]
end

function prepare_seeds(seed_list,compounds)
    [Int(i in seed_list) for i in compounds]
end

function prepare_seeds(seed_json::String,seed_key::String,compounds::Array{String,1})
    seed_dict = JSON.parsefile(seed_json)
    [Int(i in seed_dict[seed_key]) for i in compounds]
end

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

        FULLEDGEPATH = EDGEDIR*FNAME #json of all edges for organism
        FULLSEEDPATH = SEEDDIR*FNAME #json of all seeds for organism
        
        if isfile(FULLEDGEPATH)==true
            
            if split(FULLEDGEPATH,".")[2] == "json"

                println("Finding minimal seeds for: $FNAME")
                
                (R,P,compounds,reactions,t) = prepare_matrices_and_targets(FULLEDGEPATH,TARGETJSON)

                RT = transpose(R)
                PT = transpose(P)

                b = [sum(RT[i,:]) for i in 1:size(RT)[1]]
                bp = [sum(PT[i,:]) for i in 1:size(PT)[1]]

                tT = transpose(t)
                sum_t = sum(t)
                
                for (n_seed,seed_list) in enumerate(JSON.parsefile(FULLSEEDPATH))

                    OUTDIRWITHORGNAME = OUTDIR*split(FNAME,".json")[1]*"/"

                    if ispath(OUTDIRWITHORGNAME)==false
                        mkpath(OUTDIRWITHORGNAME)
                    end

                    # I want 1 randomizaiton per outpath
                    FULLOUTPATH = OUTDIRWITHORGNAME*string(n_seed)*".json"

                    seed_list_original = deepcopy(seed_list)

                    X, Y = Vector{Int}[], Vector{Int}[]
                    for cpd in seed_list_original
                        # remove cpd from seed_list
                        deleteat!(seed_list,findfirst(isequal(cpd),seed_list))
                        
                        x = prepare_seeds(seed_list,compounds)

                        X, Y = netexp(R, P, RT, PT, b, bp, x)

                        # if all targets not produced
                        if (tT*X[end]) != sum_t
                            # put cpd back in seed_list
                            push!(seed_list, cpd)
                        end
                    end

                    x = prepare_seeds(seed_list,compounds)

                    println("Writing out randomization: $n_seed")
                    simple_write_out(FULLOUTPATH, x, t, compounds, reactions, X, Y)
                    break
                end
            end
        end
        break
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
# seedkey = "Enceladus_20-SAFR-032"

# SEEDJSON = "29012812801.json" ## Contains keys numbered 1-100, with values of random compounds
TARGETJSON = "targets/Freilich09.json"
EDGEDIR = "jgi/2018-09-29/ph_edge_jsons/archaea/"
SEEDDIR = "seeds/minimal_seed_randomizations/archaea/"
OUTDIR = "results/simple/minimal_seed_randomizations_fixed/archaea/"  #*split(SEEDJSON,".json")[1]*"/"

if ispath(OUTDIR)==false
    mkpath(OUTDIR)
end

@time enumerate_minimal_seed_sets(TARGETJSON,EDGEDIR,SEEDDIR,OUTDIR)

Profile.init(n=10^7)
@profile enumerate_minimal_seed_sets(TARGETJSON,EDGEDIR,SEEDDIR,OUTDIR)
Profile.print()

########################################
#### MANY NETWORK EXPANSION RUN ######
########################################
# seedkey = "Enceladus_20-SAFR-032"

# SEEDJSON = "seeds.json"
# TARGETJSON = "targets/Freilich09.json"
# DATADIR = "jgi/2018-09-29/kegg_edge_json/"

# # fsplit = split(DATADIR,"/")
# # OUTDIR = "results/simple/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
# OUTDIR = "results/simple/kegg_edge_json/"

# if ispath(OUTDIR)==false
#     mkpath(OUTDIR)
# end

# for FNAME in readdir(DATADIR)
#     FULLINPATH = DATADIR*FNAME
#     FULLOUTPATH = OUTDIR*FNAME
#     if isfile(FULLINPATH)==true
#         if split(FULLINPATH,".")[2] == "json"
#             ## DO MAIN
#             (R,P,compounds,reactions,t) = prepare_matrices_and_targets(FULLINPATH,TARGETJSON)
#             x = prepare_seeds(SEEDJSON,seedkey,compounds)
#             (X,Y) = netexp(R,P,x)
#             simple_write_out(FULLOUTPATH,x,t,compounds,reactions,X,Y)
#         end
#     end
# end

########################################
#### SINGLE NETWORK EXPANSION RUN ######
########################################
# ## Inputs
# ds80_seeds = ["C00031","C00001"]
# reaction_edges_json = "kegg/2018-09-25/reaction_edges.json"
# target_json = "targets/Freilich09.json"
# seed_json = "seeds.json"

# ## Create out path
# fsplit = split(reaction_edges_json,"/")
# # path = "results/netexpdata_jsons/"*fsplit[end-2]*"/"*fsplit[end-1]
# path = "results/data_glucose_test/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"

# if ispath(path)==false
#     mkpath(path)
# end
# fullpath = path*"data_glucose_test.json"

# ## DO MAIN
# (R,P,compounds,reactions,t) = prepare_matrices_and_targets(reaction_edges_json,target_json)
# x = prepare_seeds(ds80_seeds,compounds)
# (X,Y) = netexp(R,P,x)
# simple_write_out(fullpath,x,t,compounds,reactions,X,Y)
