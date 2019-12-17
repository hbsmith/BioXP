using Test
using BioXP

@testset "read master from json" begin
    reactions = readmaster("../ecg/test/userdata/master.json")
    @test typeof(reactions) == Reactions
end