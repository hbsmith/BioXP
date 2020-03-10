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

#######################################

# function populatecompounds::Dict()
#     for (k,v) in Dict()
#         Symbol(k)=v
#         Mystruct(:)


# struct Point{D::Dict}
#     id::String
#     name::String
#     formula::String 
#     exact_mass::Float64
#     mol_weight::Float64
#     names::Vector{String}=Vector{String}()
#     elements::Vector{String}=Vector{String}()
#     reactions::Vector{String}=Vector{String}()   
#     comment::String=""
#     enzymes::Vector{String}=Vector{String}()   
#     rpairs::Vector{String}=Vector{String}()    
#     dblinks::Dict{String,Vector{String}}=Dict{String,Vector{String}}()     
#     pathways::Dict{String,String}=Dict{String,String}()  
#     kcf::String=""       
#     remark::String=""    
#     glycans::Vector{String}=Vector{String}()  
#     for (k,v) in D 
#         Symbol(k) = v
#     end
#     # Point{T}(x,y) where {T<:Real} = new(x,y)
# end

# struct EdgeEvidence
#     src::Int
#     dst::Int
#     evidence::Float64
#     EdgeEvidence(src, dst, evidence) = if isnan(evidence)
#         throw(DomainError(evidence, "cannot be NaN"))
#     else
#         new(src, dst, evidence)
#     end
# end

# JSON.lower(e::EdgeEvidence) = string(e)

# restore(::Type{EdgeEvidence}, j::AbstractString) = eval(Meta.parse(j))


# struct Netexp
#     reactions::a
#     targets::a
# end

# Expand(reactions::Reactions,seeds::Seeds,targets::Targets=Vector{String}()) --> produces expanded network, tracking targets (targets are none if not provided)
# Expand(reactions::Reactions,targets,minseed::MinSeed) --> produced expanded network



# #############
# const AllowedReactions = Dict{String,Dict{String,Vector{String}}}
# # cpds_left_flat =unique(Iterators.flatten([r["left"] for r in values(allowed_reactions["reactions"])]))
# const TargetCompounds = Vector{String}

# read_reactions(f::String) = AllowedReactions(JSON.parsefile(f))
# read_target(f::String) = TargetCompounds(JSON.parsefile(f))
# read_targets(f::String) = JSON.parsefile(f)





# function netexp(master_json::String,biosystem_rxns_json::String,target_compounds::String)

# end

# function netexp(master_json::String,biosystem_rxns_json::Vector,target_compounds::String)

# end

# function netexp(master_json::String,biosystem_rxns_json::Vector,target_compounds::Vector)

# end

# function netexp(master_json::String,biosystem_rxns_json::String,target_compounds::Vector)

# end

# ## Allow master_json["reactions"] dict formatting as input as well

# function netexp()
