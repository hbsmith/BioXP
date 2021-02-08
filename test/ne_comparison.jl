## Using the reactions found in Josh's network_full.csv 

@testset "Compare to python networkexpansion data" begin
    rstructs = readmaster("data/master_with_dgs_reserialized.json")
    rids_ne = readids_from_nenetwork("data/from_ne/network_full.csv")
    sids = readkeyedids("data/seeds/seeds.json")["Goldford-SAFR-032"]

    ## Check if all of Josh's reactions are in my rstructs
    length(setdiff(rids_ne,keys(rstructs))) != 0 && println("Reactions in .py ne that aren't in rstructs. This is likely due to reactions which have been removed from the KEGG database.")
    
    ## Check differences in compound lists 
    cids_ne = readdlm("data/from_ne/network_full.csv",',',header=true)[1][:, 1]
    cids = BioXP.identify_biosystem_compounds(rstructs,collect(keys(rstructs)))
    length(setdiff(cids_ne,cids)) != 0 && println("Compounds in .py ne that aren't in rstructs. This is likely due to reactions which have been removed from the KEGG database.")

    ## Check if equations compounds on both sides of the reactions are the same for reactions which are shared
    rids_shared = intersect(rids_ne,keys(rstructs))
    rid_stoich_dict_ne = Dict([(rid,Dict([("left",[]),("right",[])])) for rid in rids_shared])
    network_full = readdlm("data/from_ne/network_full.csv",',',header=true)[1]
    for row in 1:size(network_full,1)
        cid = network_full[row,:][1]
        rid = network_full[row,:][2]
        stoich = network_full[row,:][3]
        if rid in keys(rid_stoich_dict_ne)
            stoich < 0 ? push!(rid_stoich_dict_ne[rid]["left"],(cid,stoich)) : push!(rid_stoich_dict_ne[rid]["right"],(cid,stoich)) 
        end
    end

    unequal_rids = []
    for (rid, sides) in rid_stoich_dict_ne
        if Set([i[1] for i in sides["left"]]) != Set(rstructs[rid].left) && Set([i[1] for i in sides["left"]]) != Set(rstructs[rid].right)
            push!(unequal_rids,rid)
        elseif Set([i[1] for i in sides["right"]]) != Set(rstructs[rid].right) && Set([i[1] for i in sides["right"]]) != Set(rstructs[rid].left)
            push!(unequal_rids,rid)
        end
    end
    length(unequal_rids) != 0 && println("There are different reaction definitions in .py ne and rstructs. This is likely due to reactions and compounds which have been removed from the KEGG database.")

    ## Check if removing all reactions containing compounds missing from rstructs fixes issue
    ## Maybe not a good use of my time... it seems like it won't fix the issue because of reactions
    ## like R07506, that have more compounds in rstructs than they do in Josh's data.

    ## Check if filtering out elementally unbalanced reactions gives the same list of rids 

    
    ## ADD ABILITY in my code to get rid of stoichiometrically unbalanced reactions
    ## Throw away any reactions containing compounds without formula (e.g. R00002)
    ## Keep R groups, e- if balanced

    ## Check if filtering out stoichiometrically unbalanced reactions gives same list of rids


    ## If all of this works...
    ## Write out my reactions in the same format as his, to be read into his py code

    
    
    
    ## Start with 9388 unique reactions in network_full.csv
    ## 9074 reactions that are elementally balanced via `reactions_consistent.csv`
    ## 6880 reactions that are stoichiometrically balanced via `reactions_balanced.csv`

    # rids_consistent = Set(readdlm("test/data/from_ne/reactions_consistent.csv",header=true)[1])
    # rids = intersect(rids,rids_consistent)
end