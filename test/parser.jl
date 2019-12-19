# using Test
# using BioXP
using JSON

@testset "read master from file" begin
    # reactions = readmaster("../../ecg/test/userdata/kegg/master.json")
    reactions = readmaster("data/master-subset-3rxns.json")
    expected_reactions = ["R10514", "R09735", "R08827"]
    @test typeof(reactions) == Reactions
    @test Set(keys(reactions)) == Set(expected_reactions)

    @testset "correct reaction data" begin
        expected_left = [
            ["C00024","C20638"],
            ["C19848","C06232"],
            ["C00355","C08538"]]
        expected_right = [
            ["C00010","C20672"],
            ["C18237","C00020","C00001"],
            ["C08543","C00001"]]
        for (i,r) in enumerate(expected_reactions)
            @test reactions[r].left == expected_left[i]
            @test reactions[r].right == expected_right[i]
        end
    end
end

@testset "read master from open json" begin
    master = JSON.parsefile("data/master-subset-3rxns.json")
    @test typeof(master) <: Dict

    reactions = readmaster(master)
    expected_reactions = ["R10514", "R09735", "R08827"]
    @test typeof(reactions) == Reactions
    @test Set(keys(reactions)) == Set(expected_reactions)

    @testset "correct reaction data" begin
        expected_left = [
            ["C00024","C20638"],
            ["C19848","C06232"],
            ["C00355","C08538"]]
        expected_right = [
            ["C00010","C20672"],
            ["C18237","C00020","C00001"],
            ["C08543","C00001"]]
        for (i,r) in enumerate(expected_reactions)
            @test reactions[r].left == expected_left[i]
            @test reactions[r].right == expected_right[i]
        end
    end
end

@testset "read compounds from dir" begin
    