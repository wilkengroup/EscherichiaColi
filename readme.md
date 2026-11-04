# A genome-scale metabolic model of _Escherichia coli_

This repo allows you to build an _E. coli_ model that is either based on
[iML1515](https://bigg.bio/models/iML1515) or
[iCH360](https://github.com/marco-corrao/iCH360).

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

1. The `build_model` function has an optional keyword argument `simplify` that
   is set to `false` by default. When this argument is `true`, then the iCH360
   based model is loaded, otherwise the iML1515 based model is loaded. 
2. The `build_model` function has an optional keyword argument `rename_biomass`
   that is set to `false` by default. When this is set to `true`, the biomass
   reaction is renamed to `BIOMASS`.
