using Plots
using LinearAlgebra

# Domain  and data
x = 0.0:0.01:4π
y = sin.(x) .+ cos.(x)

# Visualize the data
plot(x,y,title="Test case for interpolation basis functions.")


# Radial Basis Basis - Squared Exponential basis
abstract type GenericBasis end

mutable struct SquaredExponentialBasis <: GenericBasis
    σ::Float64
    basis::Vector{Function}
    
    function SquaredExponentialBasis(D,σ)
        kernel(x,y) = exp.(-norm((x.-y)./σ^2))
        basis = []
        for i=1:length(x)
            local temp(x) = kernel(x,D[i])
            push!(basis,temp)
        end 
        new(σ,basis)
    end 
end 

function w_matrix(base::SquaredExponentialBasis,x)
    k = length(base.basis)
    @show k
    basis = base.basis
    W = zeros(length(x),k)
    for i=1:length(k)
        for j=1:length(x)
            W[j,i] = basis[i](x[j])
        end 
    end 
    W
end 

# Test
σ = 2
seb = SquaredExponentialBasis(x,σ)
Kw = w_matrix(seb,x)

# Matern Basis
mutable struct MaternBasis <: GenericBasis

end

# Piecewise Constant Basis
mutable struct PiecewiseConstantBasis <: GenericBasis


end 

# PiecewiseLinearBasis
mutable struct PiecewiseLinearBasis <: GenericBasis


end

# PiecewiseQuadraticBasis
mutable struct PiecewiseQuadraticBasis <: GenericBasis


end 