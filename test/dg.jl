## dg.jl



@testset "test simple_allowed_reactions allow_missing" begin

    @testset "some dgs" begin
        dgs = [101.1,15,-15,missing]
        threshold = -10.2

        @test [false,false,true,true]  == BioXP.simple_allowed_reactions(dgs,threshold,true)
        @test [false,false,true,false] == BioXP.simple_allowed_reactions(dgs,threshold,false)
    end

    @testset "value = threshold" begin
        dgs = [missing,missing,0,8,900]
        threshold = 8

        @test [true,true,true,true,false]   == BioXP.simple_allowed_reactions(dgs,threshold,true)
        @test [false,false,true,true,false] == BioXP.simple_allowed_reactions(dgs,threshold,false)
    end

end

@testset "test filter_reactions_by_dg" begin 

    threshold = -75
    env_key = "11pH_300mM"
    rids = readids(["R02745","R05080","R00778","R08348","R06706"])
    rstructs = readmaster("data/master_with_dgs_reserialized.json")

    @testset "unbalanced only" begin
        allow_nothings=false
        allow_unbalanced=true
        allow_within_ci=false

        forward_allowed=[false,false,true,true,false]
        backward_allowed=[false,false,false,false,true]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)

    end

    @testset "nothings and unbalanced" begin
        allow_nothings=true
        allow_unbalanced=true
        allow_within_ci=false

        forward_allowed=[true,true,true,true,false]
        backward_allowed=[true,true,false,false,true]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)
    end

    @testset "nothings only" begin
        allow_nothings=true
        allow_unbalanced=false
        allow_within_ci=false

        forward_allowed=[false,false,true,true,false]
        backward_allowed=[false,false,false,false,false]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)
    end

    @testset "ci only" begin 
        allow_nothings=false
        allow_unbalanced=false
        allow_within_ci=true

        forward_allowed=[false,false,true,true,false]
        backward_allowed=[false,false,true,false,false]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)

    end

    @testset "nothings and ci" begin #this doesn't change anything either
        allow_nothings=true
        allow_unbalanced=false
        allow_within_ci=true

        forward_allowed=[false,false,true,true,false]
        backward_allowed=[false,false,true,false,false]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)
    end

    @testset "unbalanced and ci" begin #ci doesn't open any new reactions here
        allow_nothings=false
        allow_unbalanced=true
        allow_within_ci=true

        forward_allowed=[false,false,true,true,false]
        backward_allowed=[false,false,true,false,true]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)
    end

    @testset "all allowed" begin #ci doesn't open any new reactions here
        allow_nothings=true
        allow_unbalanced=true
        allow_within_ci=true

        forward_allowed=[true,true,true,true,false]
        backward_allowed=[true,true,true,false,true]
        @test (forward_allowed, backward_allowed) == BioXP.filter_reactions_by_dg(
                                                    threshold, 
                                                    env_key, 
                                                    rids,
                                                    rstructs,
                                                    allow_nothings=allow_nothings,
                                                    allow_unbalanced=allow_unbalanced,
                                                    allow_within_ci=allow_within_ci)
    end

    ## need to test where ci allows new reaction--change threshold


end
        