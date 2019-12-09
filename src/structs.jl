module Structs

using Parameters

export Reaction, Reactions, IDs, TArray

@with_kw struct Reaction
    id::String
    left::Vector{String}
    right::Vector{String}
    metadata::Dict{String}=Dict{String,Any}()
end

### Trying to figure out how to automatically set fields here
## Compounds needed in order to set up the randomization on BioXP.

## I think I just need to make Compound take in a dict and then call things based on the dict,
## avoid trying to assign all the keys to values

@with_kw struct Compound
    id::String
    name::String
    formula::String 
    exact_mass::Float64
    mol_weight::Float64
    names::Vector{String}=Vector{String}()
    elements::Vector{String}=Vector{String}()
    reactions::Vector{String}=Vector{String}()   
    comment::String=""
    enzymes::Vector{String}=Vector{String}()   
    rpairs::Vector{String}=Vector{String}()    
    dblinks::Dict{String,Vector{String}}=Dict{String,Vector{String}}()     
    pathways::Dict{String,String}=Dict{String,String}()  
    kcf::String=""       
    remark::String=""    
    glycans::Vector{String}=Vector{String}()   
end

#consts basically will just check type, but wont actually keep constant
const IDs = Vector{String} #either "C00000" or "R00000" 

const Reactions = Dict{String,Reaction}

const Compounds = Dict{String,Compound}

# This is just cosmetic
# const TArray{T, N} = Transpose{T, Array{T, N}}

end