
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

end
