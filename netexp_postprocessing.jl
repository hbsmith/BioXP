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

fullinpath = "results/data_glucose_test/kegg/2018-09-25data_glucose_test.json"
# fullinpath = path*"data_glucose_test.json"

D = JSON.parsefile(fullinpath)
D["x"] = x
D["t"] = t
D["compounds"] = compounds
D["reactions"] = reactions
D["X"] = X
D["Y"] = Y

newdata = formatted_netexp_results(x,t,compounds,reactions,X,Y)

fsplit = split(fullinpath,"/")
outpath = "results_formatted/data_glucose_test/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
fulloutpath = outpath*"data_glucose_test.json"

## Write out
open(fulloutpath,"w") do f
    JSON.print(f, newdata)
end