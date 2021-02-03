import Random
using Random

@testset "sort_biosystem_compounds" begin

    compounds = readcompounds("data/compound")

    @testset "sortkey: exact_mass" begin

        biosystem_compounds = readids(["C00068","C00069","C00379","C00380","C00381"])
        sortkey = :exact_mass
        zero_mass_behavior = "end"

        expected_exact_mass = [("C00380",111.0433),
                               ("C00379",152.0685),
                               ("C00068",425.045),
                               ("C00381",0.0),
                               ("C00069",0.0)]

        Random.seed!(1234);
        @test expected_exact_mass == BioXP.sort_biosystem_compounds(compounds,
                                                    biosystem_compounds,
                                                    sortkey,
                                                    zero_mass_behavior)
        
        for i in 1:1000
            @test expected_exact_mass[1:3] == BioXP.sort_biosystem_compounds(compounds,
                                                    biosystem_compounds,
                                                    sortkey,
                                                    zero_mass_behavior)[1:3]

        end
        
    end

    @testset "sortkey: mol_weight" begin

        biosystem_compounds = readids(["C00068","C00069","C00379","C00380","C00381"])
        sortkey = :mol_weight
        zero_mass_behavior = "end"
        
        expected_mol_weight = [("C00380",111.102),
                               ("C00379",152.1458),
                               ("C00068",425.3144),
                               ("C00381",0.0),
                               ("C00069",0.0)]

        Random.seed!(1234);
        @test expected_mol_weight == BioXP.sort_biosystem_compounds(compounds,
                                                    biosystem_compounds,
                                                    sortkey,
                                                    zero_mass_behavior)
                                                    
        
        for i in 1:1000
            @test expected_mol_weight[1:3] == BioXP.sort_biosystem_compounds(compounds,
                                                    biosystem_compounds,
                                                    sortkey,
                                                    zero_mass_behavior)[1:3]

        end
    end

    @testset "zero_mass_behavior: random" begin

        biosystem_compounds = readids(["C00068","C00069","C00379","C00380","C00381"])
        sortkey = :exact_mass
        zero_mass_behavior = "random"

        expected_exact_mass = [("C00069",0.0),
                            ("C00380",111.0433),
                            ("C00379",152.0685),
                            ("C00381",0.0),
                            ("C00068",425.045)]
        
        expected_order = [("C00380",111.0433),
                          ("C00379",152.0685),
                          ("C00068",425.045)]

        Random.seed!(1234);
        # @test expected_exact_mass == BioXP.sort_biosystem_compounds(compounds,
        #                                             biosystem_compounds,
        #                                             sortkey,
        #                                             zero_mass_behavior)
        
        for i in 1:1000
            sorted_cpd_masses = BioXP.sort_biosystem_compounds(compounds,
                                            biosystem_compounds,
                                            sortkey,
                                            zero_mass_behavior)
            
            ## Test that the compounds with masses are always in order of inc mass
            pos = 0
            for cm in expected_order
                newpos = findfirst(isequal(cm),sorted_cpd_masses)
                @test newpos>pos
                pos = newpos
            end

        end
                            

    end

    @testset "flip probability" begin

        beta = 20

        @test .5 == BioXP.getflipprobability(("C00069",0.0),("C00066",0.0),beta)
        @test .5 == BioXP.getflipprobability(("C00069",102.0),("C00066",0.0),beta)
        @test .5 == BioXP.getflipprobability(("C00069",0.000),("C00066",6663.0),beta)
        @test .12857 ≈ BioXP.getflipprobability(("C00379",152.0685),("C00380",111.0433),beta) atol=0.00001
        @test 1.0 == BioXP.getflipprobability(("C00380",111.0433),("C00379",152.0685),beta)

    end

    @testset "randomizecompounds" begin
        biosystem_compounds = readids(["C00068","C00069","C00379","C00380","C00381"])
        n_runs = 3

        # println(randomizecompounds(biosystem_compounds,compounds,n_runs))
        runs = randomizecompounds(biosystem_compounds,compounds,n_runs)
        for i in runs
            @test length(biosystem_compounds) == length(i)
        end
    end

end