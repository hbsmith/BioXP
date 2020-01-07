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
    compounds = readcompounds("data/compound")

    @testset "spot check compound properties C00379" begin
        c = compounds["C00379"]
        @test c.entry_id == "C00379"
        @test c.name == "Xylitol"
        @test c.formula == "C5H12O5"
        @test c.exact_mass == 152.0685
        @test c.mol_weight == 152.1458
        @test Set(c.reactions) == Set([
            "R01431",
            "R01896",
            "R01904",
            "R02136",
            "R05831",
            "R07152",
            "R09477"])
        
        expected_pathways = Dict("map00040"=> "Pentose and glucuronate interconversions",
        "map01100"=> "Metabolic pathways",
        "map02010"=> "ABC transporters")
        @test Set(keys(expected_pathways)) == Set(keys(c.pathways))
        for k in keys(expected_pathways)
            @test c.pathways[k] == expected_pathways[k]
        end

        expected_enzymes = [
            "1.1.1.9",
            "1.1.1.10",
            "1.1.1.14",
            "1.1.1.15",
            "1.1.1.21",
            "1.1.1.307",
            "1.1.3.41",
            "2.7.1.122"]
        @test Set(keys(expected_enzymes)) == Set(keys(c.enzymes))
    end

    @testset "spot check compound properties C00380" begin
        c = compounds["C00380"]
        @test c.entry_id == "C00380"
        @test c.name == "Cytosine"
        @test c.formula == "C4H5N3O"
        @test c.exact_mass == 111.0433
        @test c.mol_weight == 111.102
        @test Set(c.reactions) == Set([
            "R00510",
            "R00974",
            "R02137",
            "R02296",
            "R10837"])
        
        expected_pathways = Dict("map00240"=>  "Pyrimidine metabolism",
        "map01100"=> "Metabolic pathways")
        @test Set(keys(expected_pathways)) == Set(keys(c.pathways))
        for k in keys(expected_pathways)
            @test c.pathways[k] == expected_pathways[k]
        end

        expected_enzymes = [
            "2.4.2.2",
            "2.4.2.57",
            "3.2.2.8",
            "3.2.2.10",
            "3.5.4.1"]
        @test Set(keys(expected_enzymes)) == Set(keys(c.enzymes))
    end

end
@testset "read reaction IDs" begin
    @testset "from dir" begin
        reactions = readids("data/taxon_reactions/1234567892.json")
        @test Set(reactions) == Set([
            "R01773",
            "R01775",
            "R02003",
            "R05556"])

    end

    @testset "from vector" begin
        reactions = readids([
            "R01773",
            "R01775",
            "R02003",
            "R05556"])
        
        @test Set(reactions) == Set([
            "R01773",
            "R01775",
            "R02003",
            "R05556"])

    end
end

@testset "read seed and target IDs" begin 
    @testset "from dir" begin
        seeds = readkeyedids("data/seeds/seeds.json")
        
        @test Set(seeds["Goldford-SAFR-032"]) == Set(["C00001","C00011","C00288","C00283","C00014","C00697","C00058","C00033"])
        @test Set(seeds["Enceladus_20-SAFR-032"]) == Set(["C00001","C00011","C00237","C00282","C00067","C00132","C06548","C00469","C00283","C00014","C00697","C01326","C01438","C01548","C06547","C11505","C20783","C01407"])
        @test Set(seeds["Enceladus_20-SAFR-032_P"]) == Set(["C00001","C00011","C00237","C00282","C00067","C00132","C06548","C00469","C00283","C00014","C00697","C01326","C01438","C01548","C06547","C11505","C20783","C01407","C00009"])
        @test Set(seeds["targets_Freilich09"]) == Set(["C00001","C00002","C00003","C00004","C00005","C00006","C00008","C00015","C00016","C00020","C00024","C00025","C00035","C00037","C00041","C00043","C00044","C00047","C00049","C00054","C00055","C00062","C00063","C00064","C00065","C00073","C00075","C00078","C00079","C00082","C00097","C00105","C00112","C00116","C00123","C00131","C00135","C00144","C00148","C00152","C00183","C00188","C00234","C00239","C00249","C00255","C00286","C00350","C00360","C00362","C00364","C00399","C00407","C00458","C00459","C00641","C00748","C01050","C05764","C05890","C05894","C05899","C05980","C06040","C15672","C16221"])
    end

    @testset "from dict" begin
        seeds = readkeyedids(Dict("Goldford-SAFR-032"=>["C00001","C00011","C00288","C00283","C00014","C00697","C00058","C00033"],
                 "Enceladus_20-SAFR-032"=>["C00001","C00011","C00237","C00282","C00067","C00132","C06548","C00469","C00283","C00014","C00697","C01326","C01438","C01548","C06547","C11505","C20783","C01407"]))

        @test Set(seeds["Goldford-SAFR-032"]) == Set(["C00001","C00011","C00288","C00283","C00014","C00697","C00058","C00033"])
        @test Set(seeds["Enceladus_20-SAFR-032"]) == Set(["C00001","C00011","C00237","C00282","C00067","C00132","C06548","C00469","C00283","C00014","C00697","C01326","C01438","C01548","C06547","C11505","C20783","C01407"])

    end

    @testset "from vector" begin
        seeds = readids(["C00001","C00011","C00288","C00283","C00014","C00697","C00058","C00033"])

        @test Set(seeds) == Set(["C00001","C00011","C00288","C00283","C00014","C00697","C00058","C00033"])
    end
end 

@testset "symbolize keys" begin 
    testdicts = [Dict("key1"=>1,"key2"=>2), Dict("key1"=>"1","key2"=>"2")]

    for d in testdicts
        s =  BioXP.symbolizekeys(d)
        for k in keys(s)
            @test typeof(k) == Symbol
        end
    end

end 

