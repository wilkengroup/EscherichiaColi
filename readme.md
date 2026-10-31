# A genome-scale metabolic model of _Escherichia coli_

This model is fully based on the iML1515 model, which can be found here:
https://bigg.bio/models/iML1515.

This repo merely serves as a convenient location to reproducibly build the gene
reaction stoichiometry information that is collated from multiple sources,
including the ComplexPortal and Uniprot.

The primary use-case for this repo is to build the model in a Julia-language
data format, that makes further advanced analyses possible. 


## Build instructions

If you would like to build the model from scratch, then you must first install
`EscherichiaColi` like any (unregistered) Julia package.
```julia
] add https://github.com/wilkengroup/EscherichiaColi.git
```

Thereafter, you only need to call a single function to build the entire model.
```julia
using EscherichiaColi

model = EscherichiaColi.build_model()
```

## Notes

1. The `build_model` function has an optional keyword argument
   `simplify_biomass` that is set to `false` by default. When this argument is
   `true`, then low abundance biomass constituents are removed from the biomass
   function, and it is proportionately re-weighted.
2. The `build_model` function has an optional keyword argument `rename_biomass`
   that is set to `false` by default. When this is set to `true`, the biomass
   reaction is renamed to `BIOMASS`.
