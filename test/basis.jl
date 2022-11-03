using MachineLearningPotentials
using Plots

# Define limits
xmin = 0
xmax = 1 
nintervals = 10

# Constructor
basis = PiecewiseLinearBasis([xmin,xmax],nintervals)

# Setparameters
slopes = rand(nintervals)
intercepts = rand(nintervals)
setparams!([slopes,intercepts],basis)

# Getparameters 
getparams(basis)

# Evaluate the basis
nevalpoint = 20
xeval = LinRange(xmin,xmax,nevalpoint)
yeval = basis(xeval)
plot(xeval,yeval)