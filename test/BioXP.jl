

@testset "expand crash?" begin

    rstructs = readmaster("data/master.json")
    rids = readids("data/taxon_reactions/1234567892.json")
    sids = readkeyedids("data/seeds/seeds.json")["Enceladus_20-SAFR-032"]
    tids = readkeyedids("data/seeds/seeds.json")["targets_Freilich09"]

    S = System(rstructs,
               rids,
               sids,
               tids)

    expand(S,"results/simple_test.json")

end