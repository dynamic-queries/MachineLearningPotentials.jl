using Molly
using LinearAlgebra
using Plots

function munge(VoS::Vector)
    t = length(VoS)
    N = length(VoS[1])
    p = length(VoS[1][1].data)
    trajectories = Array{Float64,3}(undef,t,N,p)
    for ts = 1:t
        for ns = 1:N
            for ps = 1:p
                trajectories[ts,ns,ps] = (VoS[ts][ns].data)[ps].val
            end 
        end 
    end 
    trajectories
end 

function LJ(natoms)
    # Parameters
    atom_mass = 10.0u"u"
    atoms = [Atom(mass=atom_mass,σ=0.3u"nm",ϵ=0.2u"kJ * mol^-1") for i=1:natoms]
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")
    coords = place_atoms(natoms, boundary; min_dist=0.3u"nm")
    temp = 100.0u"K"
    velocities = [velocity(atom_mass, temp) for i in 1:natoms]
    pairwise_inters = (LennardJones(),)

    # Define a system
    sys = System(
        atoms=atoms,
        pairwise_inters=pairwise_inters,
        coords=coords,
        velocities=velocities,
        boundary=boundary,
        loggers=(
            temp=TemperatureLogger(1),
            coords=CoordinateLogger(1),
        ),
    )

    # Define a simulator
    simulator = VelocityVerlet(
        dt=0.002u"ps",
        coupling=AndersenThermostat(temp, 1.0u"ps"),
    )

    # Perform simulation
    simulate!(sys, simulator, 1_000)

    # Get coordinates
    trajs = values(sys.loggers.coords)
    trajectories = munge(trajs)
    t = 0:0.001:1.0
    return t,trajectories
end 


function compute_velocities(trajectories,time)
    # First derivative
    vel = trajectories[2:end,:,:] .- trajectories[1:end-1,:,:]
    vel = vel[1:end,:,:] ./ (time[2:end] .- time[1:end-1])
    vel
end 

function compute_distances(trajectories)
    t,N,_ = size(trajectories)
    distances = Array{Float64,3}(undef,t,N*N,1)
    rel_coords = Array{Float64,3}(undef,t,N*N,3)
    for ts=1:t
        for ns=1:N
            for ons=1:N
                k = (ns-1)*N + ons
                temp = trajectories[ts,ns,:] - trajectories[ts,ons,:]
                rel_coords[ts,k,:] .= temp
                distances[ts,k,1] = norm(temp)
            end 
        end 
    end 
    rel_coords,distances
end 

struct RadialBasisFunction
    distances
    nmodes
    rdomain
    lengthscale
    
    function RadialBasisFunction(dist,nmodes,ε)
        rmin = minimum(dist)
        rmax = maximum(dist)
        rdomain = rmin:(rmax-rmin)/(nmodes-1):rmax
        new(dist,nmodes,rdomain,ε)
    end 
end 

function eval(rbf::RadialBasisFunction,r::Union{Vector,StepRangeLen})
    k = length(rbf.rdomain)
    rdomain = rbf.rdomain
    ε = rbf.lengthscale
    Kernel = zeros(length(r),k)
    for i=1:length(r)
        for j=1:k
            Kernel[i,j] =  exp(-((r[i]-rdomain[j])^2)/(ε)^2)
        end 
    end 
    Kernel ./= k
end 


function eval(rbf::RadialBasisFunction, r::Array{Float64,3})
    t,N2,_ = size(r)
    rv = vec(r)
    kernel = eval(rbf,rv)
    reshape(kernel,t,N2,rbf.nmodes)
end 

function ∇Potential(A::Union{Vector,Matrix},basis::RadialBasisFunction,r::Union{Array,Vector,StepRangeLen})
    nmodes = basis.nmodes
    ψ = eval(basis,r)
    reshape(ψ,(:,nmodes))*A
end 


function learn(natoms,time,rtrajectories,velocities,trajectories,basis)
    # Convert trajectories from coordinate space to distance space
    rel_coords,distance = compute_distances(trajectories)

    # Setup a radial basis function for this problem
    evaluated_hypo_space = eval(basis,distance)

    # Setup and solve regression
    nmodes = basis.nmodes
    evaluated_rel_coords = permutedims(rel_coords,(3,1,2))
    Σ = reshape(evaluated_rel_coords,(size(trajectories,3),:)) * reshape(evaluated_hypo_space,(:,nmodes))
    s = size(distance,1)*size(distance,2)
    v = reshape(permutedims(velocities,(3,1,2)),(3,:))
    A = Σ \ v
    Ar = sum(A;dims=2) / size(A,2)
    return Ar 
end


function reconstruct(A,basis,evalpoints,)

end 

# Learn the potential for 20 atoms with 1000 timesteps and 1 initial condition
natoms = 20
time,rtrajectories = LJ(natoms)
velocities = compute_velocities(rtrajectories,time)
trajectories = rtrajectories[1:end-1,:,:]
nmodes = 100
ε = 5
basis = RadialBasisFunction(distance,nmodes,ε)
Ar = learn(natoms,time,rtrajectories,velocities,trajectories,basis)

# Reconstruct the dynamics of the system
