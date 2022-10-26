module MachineLearningPotentials

    using Molly
    using UnPack

    include("interface.jl")
    include("trajectories.jl")
    include("descriptors.jl")
    include("training.jl")

    export MLPotentialProblem, solve
    export lj_homogenous, lj_heterogenous
    export transform
    export learn 
end 