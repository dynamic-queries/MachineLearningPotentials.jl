using Optim
using MachineLearningPotentials

#--------------------------------------------------------------#
# Script starts here.
#--------------------------------------------------------------#

# Generate trajectory and velocity data.
trajectories = random_data()
derivatives = random_derivative(trajectories)

# Perform optimization 
k = 20
basis = PolynomialBasis(k)
model = Maggioni(trajectories,derivatives,basis)
# result = optimize(model) 

# Non-linear Optimization

params = parameters(model)
lossfunction = loss(model)
result = optimize(lossfunction,params,LSTSQ())
#--------------------------------------------------------------#