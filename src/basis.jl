#--------------------------------------------------------------#
# Basis functions
#--------------------------------------------------------------#

abstract type AbstractBasis end 

# Hat function
mutable struct HatFunctionBasis <: AbstractBasis
    #inputs
    xmin
    xmax 
    nintervals
    xrange
    h 
    # params
    coefficients
    
    function HatFunctionBasis(limits,nintervals)
        xmin,xmax = limits
        xran = LinRange(xmin,xmax,nintervals)
        coefficients = zero(xran)
        h = xran[2]-xran[1]
        new(xmin,xmax,nintervals,xran,h,coefficients)
    end 
end 

function (basis::HatFunctionBasis)(position::Float64)
    nintervals = basis.nintervals
    coeff = basis.coefficients
    xrange = basis.xrange
    h = basis.h
    k = 1/h
    val = 0
    # 1st interval
    if(position>=xrange[1])&(position<xrange[2])
        val += k*coeff[1]*(xrange[2]-position)
    end 

    # Middle intervals
    for i=2:nintervals-1
        if(position>=xrange[i-1])&(position<xrange[i+1])
            val += k*coeff[i]*(position-xrange[i-1])
            val += k*coeff[i]*(xrange[i+1]-position)
        end 
    end   

    # Last interval
    if(position>=xrange[end-1])&(position<xrange[end])
        val += k*coeff[end]*(position-xrange[end-1])
    end     
    val
end 

function (basis::HatFunctionBasis)(positions::Union{Vector,StepRangeLen,LinRange})
    values = zero(positions)
    for (i,pos) in enumerate(positions)
        values[i] = basis(pos)
    end 
    values
end 

function setparams!(params,basis::HatFunctionBasis)
    basis.coefficients = params
    nothing
end

function getparams(basis::HatFunctionBasis)
    return basis.coefficients
end 

# Piecewise Linear function
mutable struct PiecewiseLinearBasis <: AbstractBasis
    # tunable params
    xmin
    xmax
    nintervals
    h
    # parameters
    slopes
    intercepts
    xrange

    function PiecewiseLinearBasis(limits,nintervals)
        xmin,xmax = limits
        h = (xmax - xmin) / nintervals
        sl = zeros(nintervals)
        in = zeros(nintervals)
        xran = LinRange(xmin,xmax,nintervals+1)
        new(xmin,xmax,nintervals,h,sl,in,xran)  
    end
end 

function (basis::PiecewiseLinearBasis)(position::Float64)
    nintervals = basis.nintervals
    slopes = basis.slopes
    intercepts = basis.intercepts
    xrange = basis.xrange
    val = 0
    for i=1:nintervals
        if (position >= xrange[i]) & (position < xrange[i+1])
            val += slopes[i]*position + intercepts[i]
        end 
    end 
    val
end 

function (basis::PiecewiseLinearBasis)(positions::Union{Vector,StepRangeLen,LinRange})
    values = zero(positions)
    for (i,position) in enumerate(positions)
        values[i] = basis(position)
    end 
    values
end 

function setparams!(params,basis::PiecewiseLinearBasis)
    basis.slopes,basis.intercepts = params
    nothing
end

function getparams(basis::PiecewiseLinearBasis)
    [basis.slopes,basis.intercepts]
end

using Plots


# Define limits
xmin = 0
xmax = 1 
nintervals = 10

# Constructor
basis = HatFunctionBasis([xmin,xmax],nintervals)

# Setparameters
coefficients = rand(nintervals)
setparams!(coefficients,basis)

# Getparameters 
getparams(basis)

# Evaluate the basis
nevalpoint = 20
xeval = LinRange(xmin,xmax,nevalpoint)
yeval = basis(xeval)
plot(xeval,yeval)