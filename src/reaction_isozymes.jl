
"""
Checks if reaction has a gene reaction rule assigned to it.
"""
function has_reaction_grr(model, rid)
    grr = A.reaction_gene_association_dnf(model, rid)
    if isnothing(grr) || isempty(grr) || isnothing(first(grr)) || isempty(first(grr))
        return false
    end
    return true
end

"""
Assemble reaction isozymes from Uniprot proteome data and ComplexPortal data. In
a set of complexes, this also gets rid of low confidence complexes if any
complex in the set can be found in the ComplexPortal.
"""
function get_reaction_isozymes(model, proteome_data, complex_data)
    # protein stoich map, infer from uniprot
    mer_map = Dict(
        "Homomonomer" => 1,
        "Monomer" => 1,
        "Homodimer" => 2,
        "Homotrimer" => 3,
        "Homotetramer" => 4,
        "Homopentamer" => 5,
        "Homohexamer" => 6,
        "Homoheptamer" => 7,
        "Homooctamer" => 8,
        "Homodecamer" => 10,
        "Homododecamer" => 12,
    )

    # infer protein stoichiometry from uniprot annotations
    reaction_isozymes = Dict{String,Dict{String,Isozyme}}()
    multi_component_enzymes = [] # for use later
    for rid in A.reactions(model)
        if has_reaction_grr(model, rid)
            grrs = A.reaction_gene_association_dnf(model, rid)
            for (i, grr) in enumerate(grrs)
                isozyme_id = "isozyme_" * string(i)
                d = get!(reaction_isozymes, rid, Dict{String,Isozyme}())

                if length(grr) == 1 # only assign homomers
                    gid = first(grr)
                    if haskey(proteome_data, gid) # has uniprot data
                        idx = proteome_data[gid][2]
                        ss_counts = [get(mer_map, idx, 1.0)]
                    else # no data
                        ss_counts = [1.0]
                    end
                else # assume complexes have uni-stoichiometry, fix in next step
                    ss_counts = fill(1.0, length(grr))
                    push!(multi_component_enzymes, (rid, isozyme_id))
                end

                d[isozyme_id] = Isozyme(
                    gene_product_stoichiometry = Dict(grr .=> ss_counts),
                    kcat_forward = nothing,
                    kcat_reverse = nothing
                )
            end
        end
    end

    # fix complex stoichiometry using ComplexPortal data
    fixed_multi_component_enzymes = []
    not_fixed_multi_component_enzymes = []
    for (rid, isozyme_id) in multi_component_enzymes
        grr = collect(keys(reaction_isozymes[rid][isozyme_id].gene_product_stoichiometry))

        stoichs = []
        for v in values(complex_data)
            if length(intersect(collect(keys(v)), grr)) == length(grr) &&
               length(intersect(collect(keys(v)), grr)) == length(v)
                push!(stoichs, v)
            end
        end

        if length(stoichs) == 1
            push!(fixed_multi_component_enzymes, (rid, isozyme_id))
            d = first(stoichs)
            for (k, v) in d
                if v == 0
                    d[k] = 1.0
                end
            end
            reaction_isozymes[rid][isozyme_id].gene_product_stoichiometry = d
        else
            push!(not_fixed_multi_component_enzymes, (rid, isozyme_id))
        end
    end

    # use only the entries where the stoichiometry is known, if any are known
    for (rid, isozyme_id) in not_fixed_multi_component_enzymes
        if rid in first.(fixed_multi_component_enzymes)
            delete!(reaction_isozymes[rid], isozyme_id)
        end
    end

    return reaction_isozymes
end
