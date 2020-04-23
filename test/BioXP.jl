# using Random
### This doesn't error!
# @testset "does the expand crash?" begin

#     rstructs = readmaster("data/master.json")
#     rids = readids("data/taxon_reactions/1234567892.json")
#     sids = readkeyedids("data/seeds/seeds.json")["Enceladus_20-SAFR-032"]
#     tids = readkeyedids("data/seeds/seeds.json")["targets_Freilich09"]

#     S = System(rstructs,
#                rids,
#                sids,
#                tids)

#     expand(S,"results/simple_test.json")

# end

### This doesn't error!
@testset "minimal seed sets" begin

    rstructs = readmaster("data/master.json")
    rids = readids(["R09735","R08827"])
    sids = readids(["C19848","C06232","C18237","C00020","C00001","C00355","C08538","C08543","C00001"])
    tids = readids(["C06232"])
    # sids = readkeyedids("data/seeds/seeds.json")["Enceladus_20-SAFR-032"]
    # tids = readkeyedids("data/seeds/seeds.json")["targets_Freilich09"]

    S = System(rstructs=rstructs,
               rids=rids,
               sids=sids,
               tids=tids)

    ## Below tests just make sure that nothing errors
    @test isa(find_minimal_seed_set(S,"results/simple_minimalseed_test.json"), Any)
    @test isa(find_minimal_seed_set(S), Any)

    allowed_forward = [true,false]
    allowed_backward = [false,true]

    @test isa(find_minimal_seed_set(S, nothing, allowed_forward, allowed_backward), Any)
    
    # myrng = MersenneTwister(1234)
    # [rand(myrng,[0,1]) for i in 1:length(rids)]

end

@testset "test BioXP 'private' functions" begin 

    rstructs = readmaster("data/master-subset-3rxns.json")
    rids = readids(["R09735","R08827"])

    @testset "identify_biosystem_compounds" begin 
        
        expected_cids = ["C19848","C06232","C18237","C00020","C00001","C00355","C08538","C08543","C00001"]

        @test Set(expected_cids) == Set(BioXP.identify_biosystem_compounds(rstructs,rids))
        @test 8 == length(BioXP.identify_biosystem_compounds(rstructs,rids))
    end

    @testset "matrixify_compounds" begin 

        cids = BioXP.identify_biosystem_compounds(rstructs,rids)
        expected_cid_order = ["C19848", "C06232", "C18237", "C00020", "C00001", "C00355", "C08538", "C08543"]
        # leftids = ["C19848","C06232"] ["C00355","C08538"]
        # rightids = ["C18237","C00020","C00001"] ["C08543","C00001"]
        expected_R = [1 0;1 0;0 0;0 0;0 0;0 1;0 1;0 0]
        expected_P = [0 0;0 0;1 0;1 0;1 1;0 0;0 0;0 1]

        R,P = BioXP.matrixify_compounds(rstructs,cids,rids)
        @test expected_R == R
        @test expected_P == P    

    end

    @testset "matrixify_seeds" begin 

        sids = ["C19848","C06232","C00001"]
        cids = BioXP.identify_biosystem_compounds(rstructs,rids)
        expected_cid_order = ["C19848", "C06232", "C18237", "C00020", "C00001", "C00355", "C08538", "C08543"]

        expected_x = [1,1,0,0,1,0,0,0]
        @test expected_x == BioXP.matrixify_seeds(sids, cids)

    end

    @testset "matrixify_targets" begin

        tids = ["C18237","C00020","C00001"]
        cids = BioXP.identify_biosystem_compounds(rstructs,rids)
        expected_cid_order = ["C19848", "C06232", "C18237", "C00020", "C00001", "C00355", "C08538", "C08543"]

        expected_t = [0,0,1,1,1,0,0,0]
        @test expected_t == BioXP.matrixify_targets(cids,tids)

    end

    ### test BioXP simple expansion" 

    @testset "expandmatrices" begin

        sids = ["C19848","C06232","C00001"]
        cids = BioXP.identify_biosystem_compounds(rstructs,rids)
        R,P = BioXP.matrixify_compounds(rstructs,cids,rids)
        # leftids = ["C19848","C06232"] ["C00355","C08538"]
        # rightids = ["C18237","C00020","C00001"] ["C08543","C00001"]
        # expected_cid_order = ["C19848", "C06232", "C18237", "C00020", "C00001", "C00355", "C08538", "C08543"]
        x = BioXP.matrixify_seeds(sids, cids)

        ## 1st X column should be same as seeds
        ## 2nd X column should also include new compounds 
        ##     which can be made from first column
        expected_X = [[1,1,0,0,1,0,0,0], # seeds
                      [1,1,1,1,1,0,0,0], # seeds+new compounds
                      [1,1,1,1,1,0,0,0]] # can't make anymore new compounds
        
        expected_Y = [[1,0],
                      [1,0]]

        X,Y = BioXP.expandmatrices(R,P,x)
        @test expected_X == X
        @test expected_Y == Y

    end

    @testset "expand" begin

    end

end

## Next steps:
## - check old results in `test/data/submission_results` against results from running new code
##    - the formatting should be identical. but perhaps there's no guarantee that the ordering of the compounds will be the same.
##      So I will have to figure out some way to test. 
@testset "test expanding original master prokaryotes" begin
    newformat_dir = "../test/data/submission_results/archaea_jgi_newformat/"
    archaea_simpleresults_P_dir = "../test/data/submission_results/archaea_simpleresults_P/"
    seeds_and_targets = readkeyedids("../test/data/seeds/seeds.json")

    sids = seeds_and_targets["Enceladus_20-SAFR-032_P"]
    tids = seeds_and_targets["targets_Freilich09"]
    rstructs = readmaster("../test/data/master_from_redges-og-submission.json")

    for fname in readdir(newformat_dir)
        rids = readids(newformat_dir*fname)
        x, t, cids, X, Y = expand(rstructs,rids,sids,tids) ## new results
        
        e = JSON.parsefile(archaea_simpleresults_P_dir*fname) ## expected results

        @test length(x) == length(e["x"])
        @test length(t) == length(e["t"])
        @test length(cids) == length(e["compounds"])
        @test length(X) == length(e["X"])
        @test length(Y) == length(e["Y"])

        # Compare calculated against expected results
        for gen in 1:length(e["X"])
            
            expected_cpds = Set(zip(e["compounds"],e["X"][gen]))
            calcualted_cpds = Set(zip(cids,X[gen]))
            
            @test expected_cpds == calcualted_cpds
            
            if gen != length(e["X"])
                expected_rxns = Set(zip(e["reactions"],e["Y"][gen]))
                calculated_rxns = Set(zip(rids,Y[gen]))
                
                @test expected_rxns == calculated_rxns
            end
        end
    end
end

## This test takes about 60 seconds
# @testset "test expanding original master KEGG w/P" begin
#     expected_results_fname = "../test/data/submission_results/kegg_edge_json_P/reaction_edges_P.json"
#     seeds_and_targets = readkeyedids("../test/data/seeds/seeds.json")
#     sids = seeds_and_targets["Enceladus_20-SAFR-032_P"]
#     tids = seeds_and_targets["targets_Freilich09"]
#     rstructs = readmaster("../test/data/master_from_redges-og-submission.json")
    
#     rids = readids("../test/data/taxon_reactions/kegg.json")
#     x, t, cids, X, Y = expand(rstructs,rids,sids,tids) ## new results
    
#     e = JSON.parsefile(expected_results_fname) ## expected results

#     @test length(x) == length(e["x"])
#     @test length(t) == length(e["t"])
#     @test length(cids) == length(e["compounds"])
#     @test length(X) == length(e["X"])
#     @test length(Y) == length(e["Y"])

#     # Compare calculated against expected results
#     for gen in 1:length(e["X"])
        
#         expected_cpds = Set(zip(e["compounds"],e["X"][gen]))
#         calcualted_cpds = Set(zip(cids,X[gen]))
        
#         @test expected_cpds == calcualted_cpds
        
#         if gen != length(e["X"])
#             expected_rxns = Set(zip(e["reactions"],e["Y"][gen]))
#             calculated_rxns = Set(zip(rids,Y[gen]))
            
#             @test expected_rxns == calculated_rxns
#         end
#     end
# end