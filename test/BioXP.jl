

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
