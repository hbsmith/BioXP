## Network expansion using lists/dicts
import JSON

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

function prepare_seeds(seed_json::String,seed_key::String,compounds::Array{String,1})
    seed_dict = JSON.parsefile(seed_json)
    [Int(i in seed_dict[seed_key]) for i in compounds]
end

function netexp(R::Array{Int,2},P::Array{Int,2},x::Array{Int,1})
    
    # initialize variables:
    # find the total number of metabolites in the seed set
    k = sum(x);
    # initialize previous iteration count of metabolites in the network
    k0 = 0;
    # iteration 1 consistes of the seed set
    X = []
    push!(X,x)  ## Each row is 1 generation
    # initialize reaction accumulation matrix
    Y = [];
    
    # transpose R
    RT = transpose(R)
    PT = transpose(P)
    
    b = [sum(RT[i,:]) for i in 1:size(RT)[1]]
    bp = [sum(PT[i,:]) for i in 1:size(PT)[1]]

    # while the metabolite set has not converged
    while k > k0 
        
        # update previous number of metabolites
        k0 = k;
        
        # RT*x ==> represnts the number of present metabolites in the 
        # network within each reaction; if this isequal to the total 
        # number of metabolites in each reaction, then the reaction 
        # is added to the network
        y = Array{Int}(RT * x .== b)  ## Forward reactions
        yp = Array{Int}(PT * x .== bp) ## Backward reactions

        # P*y > 0 ==> represents the vector of reactions that produce 
        # metabolite i. (i in 1:m).  If this is >0, 
        # then that metabolite is producable 
        xnew = Array{Int}(P * y .> 0)   ## Forward reactions
        xnewp = Array{Int}(R * yp .> 0) ## Backward reactions
        
        #add to previous set of metabolites (only needed to retain seed set)
        x = x .| xnew .| xnewp ## Add forward and backward reactions
        
        # find new total number of metabolites in network
        k = sum(x);

        # append accumulation matricies
        push!(X,x)
        push!(Y,y .| yp)
    
    end
    (X,Y)
end

function simple_write_out(outpath::String,x::Array{Int,1},t::Array{Int,1},compounds::Array{String,1},reactions::Array{String,1},X::Array{Any,1},Y::Array{Any,1})
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

function prepare_to_write_netexp_results(x::Array{Int,1},t::Array{Int,1},compounds::Array{String,1},reactions::Array{String,1},X::Array{Any,1},Y::Array{Any,1})
    data = Dict("stats"=>Dict(),"generations"=>Dict())

    ## Store data about the scope, about the equilibrium
    data["stats"]["scope_compounds"] = compounds
    data["stats"]["scope_reactions"] = reactions
    data["stats"]["scope_seeds"] = [compounds[i] for i in 1:length(x) if Bool(x[i])]
    data["stats"]["scope_targets"] = [compounds[i] for i in 1:length(t) if Bool(t[i])]
    
    # [Int(i in seed_list) for i in compounds]

    ## Store data by generation
    n_generations = length(X)-1 ## Last two generations are repeated for compounds
    for gen in 1:n_generations
        
        data["generations"][gen] = Dict()
        # data["generations"][gen]["reactions_cumulative"] = []
        data["generations"][gen]["compounds_cumulative"] = []
        data["generations"][gen]["targets_cumulative"] = []
        data["generations"][gen]["reactions_new"] = []
        data["generations"][gen]["compounds_new"] = []
        data["generations"][gen]["targets_new"] = []
        

        
        ## Add cumulative reaction list for each generation
        data["generations"][gen]["reactions_cumulative"] = [reactions[i] for i in 1:length(reactions) if Bool(Y[gen][i])]
        ## Alternative way to write
        # for (i,r) in enumerate(reactions)
        #     if Y[gen][i]==1
        #         push!(data["generations"][gen]["reactions_cumulative"],reactions[i])
        #     end
        # end

        ## Add cumulative compound and target list for each generation
        # data["generations"][gen]["compounds_cumulative"] = [compounds[i] for i in 1:length(compounds) if Bool(X[gen][i])]
        for (i,c) in enumerate(compounds)
            if X[gen][i]==1
                push!(data["generations"][gen]["compounds_cumulative"],compounds[i])
                if t[i]==1
                    push!(data["generations"][gen]["targets_cumulative"],compounds[i])
                end
            end
        end

        ## Store new reactions and compounds
        if gen==1
            data["generations"][gen]["reactions_new"] = collect(Set(data["generations"][gen]["reactions_cumulative"])) 
            data["generations"][gen]["compounds_new"] = collect(Set(data["generations"][gen]["compounds_cumulative"])) 
            data["generations"][gen]["targets_new"] = collect(Set(data["generations"][gen]["targets_cumulative"])) 
        end

        if gen!=1
            data["generations"][gen]["reactions_new"] = setdiff(data["generations"][gen]["reactions_cumulative"], data["generations"][gen-1]["reactions_cumulative"])
            data["generations"][gen]["compounds_new"] = setdiff(data["generations"][gen]["compounds_cumulative"], data["generations"][gen-1]["compounds_cumulative"])
            data["generations"][gen]["targets_new"] = setdiff(data["generations"][gen]["targets_cumulative"], data["generations"][gen-1]["targets_cumulative"])
        end
    end
    
    data
end

# # pass data as a json string (how it shall be displayed in a file)
# stringdata = JSON.json(dict1)
# # write the file with the stringdata variable information
# open("write_read.json", "w") do f
#     write(f, stringdata)
# end

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

## Inputs
ds80_seeds = ["C00031","C00001"]
reaction_edges_json = "kegg/2018-09-25/reaction_edges.json"
target_json = "targets/Freilich09.json"
seed_json = "seeds.json"

## Create out path
fsplit = split(reaction_edges_json,"/")
# path = "results/netexpdata_jsons/"*fsplit[end-2]*"/"*fsplit[end-1]
path = "results/data_glucose_test/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"

if ispath(path)==false
    mkpath(path)
end
fullpath = path*"data_glucose_test.json"

## DO MAIN
(R,P,compounds,reactions,t) = prepare_matrices_and_targets(reaction_edges_json,target_json)
x = prepare_seeds(ds80_seeds,compounds)
(X,Y) = netexp(R,P,x)
simple_write_out(fullpath,x,t,compounds,reactions,X,Y)


## Not using
# data = prepare_to_write_netexp_results(x,t,compounds,reactions,X,Y)

# ## Write out
# open(path*"data_glucose_test.json","w") do f
#     JSON.print(f, data)
# end