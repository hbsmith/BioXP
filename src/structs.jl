module Structs

using Parameters

export Reaction, Reactions, Compound, Compounds, IDs

@with_kw struct System 
    reaction_structs::Reactions 
    biosystem_reactions::IDs
    seed_compounds::IDs 
    target_compounds::IDs = IDs()
end

@with_kw struct Reaction
    id::String
    left::Vector{String}
    right::Vector{String}
    metadata::Dict{String}=Dict{String,Any}()
end

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
const TArray{T, N} = Transpose{T, Array{T, N}}

end