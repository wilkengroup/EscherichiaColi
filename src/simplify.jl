
"""
Return the chemical element of `x`.
"""
to_element(x::String) = begin
    sym = length(x) > 1 ? Symbol(uppercase(first(x)) * x[2:end]) : Symbol(uppercase(first(x)))
    elements[sym]
end

"""
Return molar mass of `mid`. Return `NaN` if formula does not exist.
"""
get_molar_mass(model, mid) = begin
    rs = A.metabolite_formula(model, mid)
    isnothing(rs) && return NaN # if metabolite has no molar mass
    sum(v * to_element(k).atomic_mass for (k, v) in rs).val
end

"""
Remove low fraction biomass constituents. Deleting these from the biomass
function might help reduce numerical trouble with various analyses, like
gradient descent etc.
"""
function simplify_biomass!(model)

    delete!(model.reactions, "BIOMASS_Ec_iML1515_WT_75p37M")

    # list of biomass components to delete:
    bmids = [
        "10fthf_c" # -0.000223,
        "2fe2s_c" # -2.6e-05,
        "2ohph_c" # -0.000223,
        "4fe4s_c" # -0.00026,
        "amet_c" # -0.000223,
        "btn_c" # -2e-06,
        "ca2_c" # -0.005205,
        "cl_c" # -0.005205,
        "coa_c" # -0.000576,
        "cobalt2_c" # -2.5e-05,
        "cu2_c" # -0.000709,
        "fad_c" # -0.000223,
        "fe2_c" # -0.006715,
        "fe3_c" # -0.007808,
        "k_c" # -0.195193,
        "mg2_c" # -0.008675,
        "mlthf_c" # -0.000223,
        "mn2_c" # -0.000691,
        "mobd_c" # -7e-06,
        "nad_c" # -0.001831,
        "nadp_c" # -0.000447,
        "ni2_c" # -0.000323,
        "pheme_c" # -0.000223,
        "pydx5p_c" # -0.000223,
        "ribflv_c" # -0.000223,
        "sheme_c" # -0.000223,
        "so4_c" # -0.004338,
        "succoa_c" # -9.8e-05,
        "thf_c" # -0.000223,
        "thmpp_c" # -0.000223,
        "udcpdp_c" # -5.5e-05,
        "zn2_c" # -0.000341
    ]
    for mid in bmids
        delete!(model.reactions["BIOMASS_Ec_iML1515_core_75p37M"].stoichiometry, mid)
    end


    d = A.reaction_stoichiometry(model, "BIOMASS_Ec_iML1515_core_75p37M")
    atp_baseline = d["adp_c"]
    atprxn = ["adp_c", "atp_c", "h_c", "h2o_c", "pi_c"]
    for mid in keys(d)
        if mid in atprxn
            d[mid] += (atp_baseline * (mid in ["atp_c", "h2o_c"] ? +1 : -1))
        end
    end

    x = sum(v * get_molar_mass(model, k) for (k, v) in d) / -1000.0
    reweight = 1/x

    for mid in keys(d)
        if mid in atprxn
            d[mid] *= reweight 
            d[mid] -= (mid in ["atp_c", "h2o_c"] ? +1 : -1)*atpm_baseline
        else
            d[mid] *= reweight
        end
    end

end
