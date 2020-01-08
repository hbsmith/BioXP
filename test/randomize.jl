@testset "sort_biosystem_compounds" begin

    compounds = readcompounds("data/compound")

    @testset "test sortkey" begin