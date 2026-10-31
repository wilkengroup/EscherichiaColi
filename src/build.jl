
function build_model(simplify_biomass = false)

    model = convert(CM.Model, A.load(joinpath(pkgdir(@__MODULE__),"data", "iML1515.json")))
    if simplify_biomass
        simplify_biomass!(model)
    end

    # load data from complex DB, and uniprot
    proteome_data = begin
        x = CSV.File(joinpath(pkgdir(@__MODULE__), "data", "e_coli_uniprot.tsv"), delim = '\t')
        Dict{String,Tuple{Float64,String}}(
            gp => (m, coalesce(k, "")) for (gp, m, k) in zip(x.gene_product, x.mass, x.kind)
        )
    end

    complex_data = begin
        x = CSV.File(joinpath(pkgdir(@__MODULE__), "data", "e_coli_complex.tsv"), delim = '\t')
        res = Dict{String,Dict{String,Float64}}()
        for c in x.complex
            res[c] = Dict{String,Float64}()
        end
        for (c, gp, s) in zip(x.complex, x.gene_product, x.stoichiometry)
            res[c][gp] = s
        end
        res
    end

    reaction_isozymes = get_reaction_isozymes(model, proteome_data, complex_data)

    ecoli = Model()

    for (rid, rxn) in model.reactions
        ecoli.reactions[rid] = Reaction(
            name = rxn.name,
            lower_bound = rxn.lower_bound,
            upper_bound = rxn.upper_bound,
            stoichiometry = rxn.stoichiometry,
            objective_coefficient = rxn.objective_coefficient,
            gene_association = get(reaction_isozymes, rid, nothing),
            annotations = rxn.annotations,
            notes = rxn.notes,
        )
    end

    for (mid, met) in model.metabolites
        ecoli.metabolites[mid] = Metabolite(
            name = met.name,
            compartment = met.compartment,
            formula = met.formula,
            charge = met.charge,
            annotations = met.annotations,
            notes = met.notes,
        )
    end

    for (gid, gn) in model.genes
        ecoli.genes[gid] = Gene(
            name = gn.name,
            annotations = gn.annotations,
            notes = gn.notes,
        )
    end

    ecoli
end
