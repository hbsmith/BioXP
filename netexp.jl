## Network expansion using lists/dicts
import JSON

reaction_edges = JSON.parsefile("../kegg/2018-09-25/reaction_edges.json");
reactions = [i for i in keys(reaction_edges["products"])] ## could use ["substrates""] too

cpd_reactants_flat = unique(Iterators.flatten(values(reaction_edges["substrates"])))
cpd_products_flat = unique(Iterators.flatten(values(reaction_edges["products"])))
compounds = unique(Iterators.flatten([cpd_reactants_flat,cpd_products_flat]));


## For R and P:
## Reactons in columns
## Compounds in rows
## Rows, Columns
## This is equivilent to R and P in segrelab's github
R = zeros(Int, (length(compounds),length(reactions)))
P = zeros(Int, (length(compounds),length(reactions)));

for (i,r) in enumerate(reactions)
    R[:,i] = [Int(cpd in reaction_edges["substrates"][r]) for cpd in compounds]
    P[:,i] = [Int(cpd in reaction_edges["products"][r]) for cpd in compounds]
end

# seed_compounds = JSON.parsefile("../seeds.json");
# x = [Int(i in seed_compounds["Enceladus_20-SAFR-032"]) for i in compounds];

# ds80_seeds = ["C00011",
# "C20298",
# "C14819",
# "C00087",
# "C00237",
# "C00058",
# "C00033",
# "C00031",
# "C00095",
# "C00124",
# "C00159",
# "C00243",
# "C00208",
# "C00282",
# "C00007",
# "C00001"]

ds80_seeds = ["C00031","C00001"]

x = [Int(i in ds80_seeds) for i in compounds];


target_compounds = JSON.parsefile("../targets/Freilich09.json");
target_compounds_org = collect(intersect(keys(target_compounds),compounds));
t = [Int(i in target_compounds_org) for i in compounds];

function netexp(R::Array{Int,2},P::Array{Int,2},x::Array{Int,1})
    
    # initialize variables:
    # find the total number of metabolites in the seed set
    k = sum(x);
    # initialize previous iteration count of metabolites in the network
    k0 = 0;
    # iteration 1 consistes of the seed set
    X = []
    push!(X,x)  ## Each row is 1 generation
    # initialize reaction accumulation matrix
    Y = [];
    
    # transpose R
    RT = transpose(R)
    PT = transpose(P)
    
    b = [sum(RT[i,:]) for i in 1:size(RT)[1]]
    bp = [sum(PT[i,:]) for i in 1:size(PT)[1]]

    # while the metabolite set has not converged
    while k > k0 
        
        # update previous number of metabolites
        k0 = k;
        
        # RT*x ==> represnts the number of present metabolites in the 
        # network within each reaction; if this isequal to the total 
        # number of metabolites in each reaction, then the reaction 
        # is added to the network
        y = Array{Int}(RT * x .== b)  ## Forward reactions
        yp = Array{Int}(PT * x .== bp) ## Backward reactions

        # P*y > 0 ==> represents the vector of reactions that produce 
        # metabolite i. (i in 1:m).  If this is >0, 
        # then that metabolite is producable 
        xnew = Array{Int}(P * y .> 0)   ## Forward reactions
        xnewp = Array{Int}(R * yp .> 0) ## Backward reactions
        
        #add to previous set of metabolites (only needed to retain seed set)
        x = x .| xnew .| xnewp ## Add forward and backward reactions
        
        # find new total number of metabolites in network
        k = sum(x);

        # append accumulation matricies
        push!(X,x)
        push!(Y,y .| yp)
    
    end
    (X,Y)
end
    

(X,Y) = netexp(R,P,x)

size(X)

size(Y)

X

for i in 1:length(X)
    println(sum(X[i]))
end

for i in 1:length(Y)
    println(sum(Y[i]))
end

# unique([Y[1][i]==1 ? reactions[i] : "SKIP" for (i,r) in enumerate(reactions)])


length(reactions)

length(compounds)

reactions_by_generation = Dict()
for j in 1:size(Y)[1]
    reactions_by_generation[j] = []
    for (i,r) in enumerate(reactions)
        if Y[j][i]==1
            push!(reactions_by_generation[j],reactions[i])
        end
    end
end

compounds_by_generation = Dict()
for j in 1:size(X)[1]
    compounds_by_generation[j] = []
    for (i,c) in enumerate(compounds)
        if X[j][i]==1
            push!(compounds_by_generation[j],compounds[i])
        end
    end
end

reactions_by_generation

compounds_by_generation
