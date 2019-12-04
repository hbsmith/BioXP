module Structs

using Parameters

export Reaction, Reactions, IDs, TArray

@with_kw struct Reaction
    id::String
    left::Vector{String}
    right::Vector{String}
    metadata::Dict{String}=Dict{String,Any}()
end

#consts basically will just check type, but wont actually keep constant
const IDs = Vector{String} #either "C00000" or "R00000" 

const Reactions = Dict{String,Reaction}

# This is just cosmetic
const TArray{T, N} = Transpose{T, Array{T, N}}

end