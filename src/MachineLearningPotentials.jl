module MachineLearningPotentials

    using Molly
    using UnPack
    using LinearAlgebra
    using Optim
    

    include("utils.jl")
    
    include("basis.jl")
    export HatFunctionBasis,PiecewiseLinearBasis,getparams,setparams!
    
    include("interface.jl")
    export MLPotentialProblem, solve

    include("datagen.jl")
    export random_data,random_derivative
    export lj_homogenous, lj_heterogenous

    include("models.jl")
    export Maggioni,optimize,parameters,loss
end 