## dg.jl

# @testset "test filter_reactions_by_dg" begin 

#     rstructs = readmaster("data/master_with_dgs_reserialized.json")
#     rids = readids(["R02745","R05080","R00778","R08348","R06706"])
#     env_key = "11pH_300mM"

#     @testset ""

@testset "test simple_allowed_reactions" begin

    dgs = [101.1,15,-15,missing]
    threshold = -10.2

    @testset "allow_missing=true" begin
        expected_output = [false,false,true,true]
        @test expected_output == BioXP.simple_allowed_reactions(dgs,threshold,true)
    end
    # @testset "allow_missing=false"

    # @testset "identify_biosystem_compounds" begin 
        
    #     expected_cids = ["C19848","C06232","C18237","C00020","C00001","C00355","C08538","C08543","C00001"]

    #     @test Set(expected_cids) == Set(BioXP.identify_biosystem_compounds(rstructs,rids))
    #     @test 8 == length(BioXP.identify_biosystem_compounds(rstructs,rids))
    # end
end