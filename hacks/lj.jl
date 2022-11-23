using OrdinaryDiffEq
using Plots
using LinearAlgebra
using ForwardDiff
using StatsBase


# LJ potential
function LJPotential(r,ε,σ)
    4*ε*((σ/r)^12 - (σ/r)^6)
end 

function ∇LJPotential(r,ε,σ)
    24*ε*(((σ)^6/(r)^7)- 2*((σ)^12/(r)^13))
end 

# Setup initial conditions
# We assume a unit cube, with a resolution of 1e-6
Natoms = 100
ndims = 3
domain = 0.0:1e-6:1.0
X = sample(domain,(Natoms,ndims),replace=false)
scatter3d(X[:,1],X[:,2],X[:,3],markersize=1,title="Initial positions of Ar atoms in a unit cube")
savefig("figures/initial_positions.png")

# Compute displacements
function displacements(X::Matrix)
    Natoms = size(X,1)
    disp = zeros(Natoms*(Natoms-1),3)
    k = 1
    for i=1:Natoms
        for j=1:Natoms
            if i!=j
                disp[k,:] .= X[i,:] .- X[j,:]
                k += 1 
            end 
        end 
    end
    disp 
end

disp = displacements(X) 

# Compute distances
function distances(X::Matrix)
    n = size(X,1)
    distances = zeros(n)
    for i=1:n
        distances[i] = norm(X[i,:])
    end
    distances 
end 
dist = distances(disp)
scatter(dist[:,1],markersize=0.2,title="Distances between the molecules in the domain.",ylabel="Distances in A",xlabel="Atom pairs")
savefig("figures/distances.png")


# Compute the LJ potential
R = sort(dist) 
# Parameters from molmod Xe I data
ϵ = 32.30 # 1/(kᵦ K)
σ = 2.82 # 1/A
U = LJPotential.(R,ϵ,σ)
∇U = ∇LJPotential.(R,ϵ,σ)

# Visualize the potential and its derivative.
plot(R,U,title="LJ Potential for Ne.",xlabel="rᵢⱼ",ylabel="Vₗⱼ")
savefig("figures/LJ_neon.png")
plot(R,∇U,xlabel="rᵢⱼ",ylabel="∇Vₗⱼ",label="Analytical",title="Gradient of LJ Potential for neon")
savefig("figures/∇LJ_neon.png")