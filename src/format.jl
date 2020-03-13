import JSON

## As simple as...
## (run netexp)
## formatbioxpoutput(netexp_simple_results_dir)
## OR
## formatbioxpoutput(netexp_simple_results_file)

## Single file
function formatbioxpoutput(path::String,
                           write_path::Union{String,Nothing}=nothing)

    if isdir(path)
        _formatbioxpdir(path,write_path)
    elseif isfile(path)
        _formatbioxpfile(path,write_path)
    else
        error("neither a file nor a directory")
    end

end

function _formatbioxpfile(fname::String,
                         write_path::Union{String,Nothing}=nothing)

    ## write_path is guaranteed to be a file here, named "formatted/fname.ext"
    if write_path == nothing
        # path = joinpath(splitpath(abspath(path))[1:end-2]...)
        newdir = joinpath(dirname(abspath(fname)),"formatted")
        if !ispath(newdir)
            mkpath(newdir)
        end
        write_path = joinpath(newdir,basename(fname))
    ## in case write_path is a file, get the directory
    elseif !ispath(dirname(abspath(write_path)))
        mkpath(dirname(abspath(write_path)))
    end

    formattedexpansion = formatexpansion(readexpansion(fname))

    open(write_path,"w") do f
        JSON.print(f, formattedexpansion, 2) #indent=2
    end
end

function _formatbioxpdir(dir::String,
                        write_path::Union{String,Nothing}=nothing)
    
    ## write_path is guaranteed to be a dir here, named "formatted"
    if write_path == nothing
        newdir = joinpath(abspath(dir),"formatted")
        if !ispath(newdir)
            mkpath(newdir)
        end
    ## This can only possibly create a directory
    elseif !ispath(write_path)
        mkpath(write_path)
    end
    
    # "$i.json"
    for fname in filter(x->endswith(x,".json"),readdir(abspath(dir))) ## List all json files in dir
        
        formattedexpansion = formatexpansion(readexpansion(joinpath(abspath(dir),fname)))

        if write_path == nothing
            real_write_path = joinpath(newdir,basename(fname))  
        ## Only supports choosing the dirname to write when writing whole dir
        else
            real_write_path = joinpath(write_path,basename(fname))
        end
    
        open(real_write_path,"w") do f
            JSON.print(f, formattedexpansion, 2) #indent=2
        end
    end

end


function formatexpansion(nexp::Expansion)
    formatexpansion(nexp.x,
                    nexp.t,
                    nexp.compounds,
                    nexp.reactions,
                    nexp.X,
                    nexp.Y)
end

function formatexpansion(x::Vector{Int},t::Vector{Int},compounds::Array{String,1},reactions::Array{String,1},X::Vector{Vector{Int}},Y::Vector{Vector{Int}})
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