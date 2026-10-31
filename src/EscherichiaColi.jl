module EscherichiaColi

using CSV, JSON
using DocStringExtensions
using AbstractFBCModels ,JSONFBCModels
using PeriodicTable # to correctly simplify biomass

import AbstractFBCModels as A
import .A.CanonicalModel as CM

include("model.jl")
include("simplify.jl")
include("reaction_isozymes.jl")
include("build.jl")

end
