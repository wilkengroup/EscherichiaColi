module EscherichiaColi

using CSV, JSON
import AbstractFBCModels as A
using DocStringExtensions

include("model.jl")
include("reaction_isozymes.jl")
include("build.jl")

end
