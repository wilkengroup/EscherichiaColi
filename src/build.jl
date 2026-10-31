
# load data from complex DB, and uniprot
proteome_data = begin
    x = CSV.File("data/ecoli/e_coli_uniprot.tsv", delim = '\t')
    Dict{String,Tuple{Float64,String}}(
        gp => (m, coalesce(k, "")) for (gp, m, k) in zip(x.gene_product, x.mass, x.kind)
    )
end

complex_data = begin
    x = CSV.File("data/ecoli/e_coli_complex.tsv", delim = '\t')
    res = Dict{String,Dict{String,Float64}}()
    for c in x.complex
        res[c] = Dict{String,Float64}()
    end
    for (c, gp, s) in zip(x.complex, x.gene_product, x.stoichiometry)
        res[c][gp] = s
    end
    res
end

reaction_isozymes = get_reaction_isozymes!(model, proteome_data, complex_data)
