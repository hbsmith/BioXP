using Parameters, LinearAlgebra

# export System, Reaction, Reactions, Compound, Compounds, IDs

@with_kw struct Reaction
    id::String
    left::Vector{String}
    right::Vector{String}
    metadata::Dict{String}=Dict{String,Any}()
end

@with_kw struct Compound
    entry_id::String
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

@with_kw struct System 
    rstructs::Reactions 
    rids::IDs
    sids::IDs 
    tids::IDs = IDs()
end

struct Expansion
    x::Vector{Int}
    t::Vector{Int}
    compounds::IDs
    reactions::IDs
    X::Vector{Vector{Int}}
    Y::Vector{Vector{Int}}
end