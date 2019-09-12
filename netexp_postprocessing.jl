import JSON

function formatted_netexp_results(x::Array{Int,1},t::Array{Int,1},compounds::Array{String,1},reactions::Array{String,1},X::Array{Any,1},Y::Array{Any,1})
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



function format_many(DATADIR::String,OUTDIR::String)

    for FNAME in readdir(DATADIR)
        
        FULLINPATH = DATADIR*FNAME
        FULLOUTPATH = OUTDIR*FNAME


        D = JSON.parsefile(FULLINPATH)
        x = Array{Int,1}(D["x"])
        t = Array{Int,1}(D["t"])
        compounds = Array{String,1}(D["compounds"])
        reactions = Array{String,1}(D["reactions"])
        X = Array{Any,1}(D["X"])
        Y = Array{Any,1}(D["Y"])

        newdata = formatted_netexp_results(x,t,compounds,reactions,X,Y)

        ## Write out
        open(FULLOUTPATH,"w") do f
            JSON.print(f, newdata)
        end
    
    end

end

function format_many_nested(DOMAINDIR::String,OUTDIR::String)

    for orgdir in readdir(DOMAINDIR)

        if !ispath(joinpath(OUTDIR,orgdir))
            mkpath(joinpath(OUTDIR,orgdir))
        end

        for seedname in readdir(joinpath(DOMAINDIR,orgdir))

            seedpath_in = joinpath(DOMAINDIR,orgdir,seedname)
            seedpath_out = joinpath(OUTDIR,orgdir,seedname)

            if isfile(seedpath_in) && (last(splitext(seedpath_in)) == ".json")

                D = JSON.parsefile(seedpath_in)
                x = Array{Int,1}(D["x"])
                t = Array{Int,1}(D["t"])
                compounds = Array{String,1}(D["compounds"])
                reactions = Array{String,1}(D["reactions"])
                X = Array{Any,1}(D["X"])
                Y = Array{Any,1}(D["Y"])
        
                newdata = formatted_netexp_results(x,t,compounds,reactions,X,Y)
        
                ## Write out
                open(seedpath_out,"w") do f
                    JSON.print(f, newdata)
                end
            
            end
        
        end
    
    end
    
end
            
#########################
### FORMAT MANY NESTED FILES
#########################
# for domain in ["bacteria"]
#     DATADIR = "results/simple_2019-09-09/min_seeds_partial/"*domain*"/"
#     OUTDIR = "results/formatted_2019-09-09/min_seeds_partial/"*domain*"/"

#     format_many_nested(DATADIR,OUTDIR)

# end

#########################
### FORMAT MANY FILES
#########################
# const DATADIR = "results/simple/min_seeds_partial/archaea/2506520044/"

# # fsplit = split(DATADIR,"/")
# # OUTDIR = "results/formatted/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
# const OUTDIR = "results/formatted/min_seeds_partial/archaea/2506520044/"


for domain in ["archaea","bacteria"]
    DATADIR = "results/simple_2019-09-09/minimal_seed_randomizations_together/"*domain*"/"
    OUTDIR = "results/formatted_2019-09-09/minimal_seed_randomizations_together/"*domain*"/"

    if !ispath(OUTDIR)
        mkpath(OUTDIR)
    end

    format_many(DATADIR,OUTDIR)

end

#########################
### FORMAT SINGLE FILE
#########################
# fullinpath = "results/simple_2019-09-09/kegg_edge_json/reaction_edges.json"
# # fullinpath = path*"data_glucose_test.json"

# D = JSON.parsefile(fullinpath)
# x = Array{Int,1}(D["x"])
# t = Array{Int,1}(D["t"])
# compounds = Array{String,1}(D["compounds"])
# reactions = Array{String,1}(D["reactions"])
# X = Array{Any,1}(D["X"])
# Y = Array{Any,1}(D["Y"])

# newdata = formatted_netexp_results(x,t,compounds,reactions,X,Y)

# fsplit = split(fullinpath,"/")
# # outpath = "results/formatted/"*fsplit[end-2]*"/"*fsplit[end-1]*"/"
# outpath = "results/formatted_2019-09-09/kegg_edge_json/"
# if !ispath(outpath)
#     mkpath(outpath)
# end
# fulloutpath = outpath*"reaction_edges.json"

# ## Write out
# open(fulloutpath,"w") do f
#     JSON.print(f, newdata)
# end