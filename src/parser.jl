import JSON

# export readmaster, readcompounds, readids

function readexpansion(f::String)
    nexp = JSON.parsefile(f)
    x = Vector{Int}(nexp["x"])
    t = Vector{Int}(nexp["t"])
    compounds = IDs(nexp["compounds"])
    reactions = IDs(nexp["reactions"])
    X = Vector{Vector{Int}}(nexp["X"])
    Y = Vector{Vector{Int}}(nexp["Y"])
    Expansion(x,t,compounds,reactions,X,Y)
end

## The point of this is to verify the structure of the 
## master file or object.
"""
    readmaster(master)

Return a `Dict{Any,Reaction}` of reactions coerced into Reaction struct.
When applied to a string, look for a file to read in. 
"""
function readmaster end

function readmaster(f::String)
    master = JSON.parsefile(f)
    readmaster(master)
end

function readmaster(master::Dict)
    reactions = Reactions()
    for (k,v) in master["reactions"]
        reactions[k] = Reaction(id=k,left=v["left"],right=v["right"],metadata=v["metadata"])
    end
    return reactions
end

"""
    readcompounds(dir)

Returns a `Dict{Symbol,Compound}` of compounds coerced into Compound struct.
"""
function readcompounds(dir::String)
    compounds = Compounds()
    for f in readdir(dir)
        if endswith(f,".json")
            compound = JSON.parsefile(joinpath(dir,f))[1]
            compound = symbolizekeys(compound)
            compounds[compound[:entry_id]] = Compound(;compound...)
        end
    end
    return compounds
end

"""
    readids(f)

Return a `Vector{String}` from a vector of reactions or compounds.
When applied to a string, look for a file to read in. 
"""
function readids end

function readids(f::String)
    ids = JSON.parsefile(f)
    readids(ids)
end

function readids(ids::Vector{<:Any})
    return convert(IDs,ids)
end

"""
    readkeyedids(f)

Return a `Dict{String,IDs}` from a dict-like with values of reactions or compounds.
When applied to a string, look for a file to read in. 
"""
function readkeyedids(f::String)
    ids = JSON.parsefile(f)
    readkeyedids(ids)
end

function readkeyedids(ids::Dict)
    return convert(Dict{String,IDs},ids)
end

"""
    symbolizekeys(d)

Converts the keys of a dict to the type `Symbol`.
"""
function symbolizekeys(d::Dict)
    newdict = Dict{Symbol,Any}()
    for (k,v) in d
        newdict[Symbol(k)] = v
    end
    return newdict
end