# Base.:+(a::Nothing,b::Nothing) = nothing
# Base.:-(a::Nothing,b::Nothing) = nothing

function simple_allowed_reactions(
    dgs::Vector{<:Union{Missing,Real}}, #this shouldn't contain `nothing`s though
    threshold::Real,
    allow_missing::Bool=false)

    allowed_reactions = dgs .<= threshold
    if allow_missing==true
        allowed_reactions[ismissing.(dgs)] .= true
    else
        allowed_reactions[ismissing.(dgs)] .= false
    end
    
    Vector{Bool}(allowed_reactions)
end

"""
Indicates which reactions can proceed based on the dg threshold.

Returns two Vector{Bool} of length equal to rids, one for forward reactions and one for backwards.

Note: `allow_unbalanced==false` will only allow reactions with `is_balanced==nothing`
      when `allow_nothings==true`. 
      `allow_unbalanced==true` will only allow reactions with `is_balanced==nothing`
      when `allow_nothings==true`
"""
function filter_reactions_by_dg(
    threshold::Real, 
    env_key::Any, 
    rids::IDs,
    rstructs::Reactions;
    allow_nothings::Bool=false,
    allow_unbalanced::Bool=false,
    allow_within_ci::Bool=true)

    mean_dgs = Union{Nothing,Missing,Real}[rstructs[rid].metadata["dg"][env_key]["standard_dg_prime_value"] for rid in rids]
    cis = Union{Nothing,Missing,Real}[rstructs[rid].metadata["dg"][env_key]["standard_dg_prime_ci"] for rid in rids]
    ## convert to missing because it can handle operators
    mean_dgs[mean_dgs .== nothing] .= missing
    cis[cis .== nothing] .= missing

    ## identify best case scenario dg values
    if allow_within_ci==true
        forward_dgs = Vector{Union{Missing,Real}}(mean_dgs - cis)
        backward_dgs = Vector{Union{Missing,Real}}(-mean_dgs - cis)
    else
        forward_dgs = Vector{Union{Missing,Real}}(mean_dgs)
        backward_dgs = Vector{Union{Missing,Real}}(-mean_dgs)
    end

    forward_allowed = simple_allowed_reactions(forward_dgs,threshold,allow_nothings)
    backward_allowed = simple_allowed_reactions(backward_dgs,threshold,allow_nothings)

    ## gather balance data
    balances = Union{Nothing,Missing,Bool}[rstructs[rid].metadata["dg"][env_key]["is_balanced"] for rid in rids]
    ## convert to missing because it can handle operators
    balances[balances .== nothing] .= missing

    if (allow_unbalanced == true) & (allow_nothings == true)
        ## don't use balance data
        return forward_allowed, backward_allowed

    elseif (allow_unbalanced == true) & (allow_nothings == false)
        forward_allowed[ismissing.(balances)] .= false
        backward_allowed[ismissing.(balances)] .= false
        return forward_allowed, backward_allowed

    elseif (allow_unbalanced == false) & (allow_nothings == true)
        balances[ismissing.(balances)] .= true

        forward_allowed[balances .== false] .= false
        backward_allowed[balances .== false] .= false
        return forward_allowed, backward_allowed

    elseif (allow_unbalanced == false) & (allow_nothings == false)
        balances[ismissing.(balances)] .= false

        forward_allowed[balances .== false] .= false
        backward_allowed[balances .== false] .= false
        return forward_allowed, backward_allowed
    end
end
    
    

