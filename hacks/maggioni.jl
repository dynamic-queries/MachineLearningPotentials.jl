using Molly
using Molly:VelocityVerlet
using LinearAlgebra
using Plots
using SparseArrays

# Use Molly's LJ to generate trajectories
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

# Compute velocities
function velocities(X,t)
    nsteps = size(X,1)
    shape = collect(size(X))
    shape[1] = shape[1]-1
    vel = zeros(shape...)
    dt = t[2:end] - t[1:end-1]
    for step = 1:nsteps-1
        vel[step,:,:] .= (X[step+1,:,:] .- X[step,:,:])/(dt[step])
    end     
    permutedims(vel,(3,2,1))
end 

# Change of basis matrix
function transition_matrix(r)
    ϕ(r1,r2) = exp.(-((r1-r2)^2)/σ)
    N = length(r)
    K = spzeros(N,N)
    for i=1:N
        for j=1:N
            K[i,j] = ϕ(r[i],r[j])
        end
    end 
    K
end

# Compute the displacement and distance of the atoms
function displacement_and_distance(X)
    _,N,L = size(X)
    Rij = spzeros(3*N,N*N)
    rij = []
    for i=1:N
        hs = (i-1)*N + 1
        he = hs + N - 1
        vs = (i-1)*3 + 1
        ve = vs + 2
        Rij[vs:ve,hs:he] = displacement(X[:,i,1],X[:,:,1])
        append!(rij,distance(Rij[vs:ve,hs:he]))
    end 
    Rij = (1/N)*Rij
    Rij,rij
end

function distance(X)
    n = size(X,2)
    Z = spzeros(n)
    for i=1:n
        Z[i] = norm(X[:,i])
    end 
    Z
end 

function displacement(X,Y)
    ndims = length(X)
    N = size(Y,2)
    M = spzeros(ndims,N)
    for i=1:N 
        M[:,i] = Y[:,i] - X
    end 
    M
end 


function displacement_and_distance(X,t)
    _,N,L = size(X)
    @show N,L
    Rij = spzeros(3*N*L,N*N*L)
    rij = []
    for k=1:L
        for i=1:N
            hs = (k-1)*N*N + (i-1)*N + 1
            he = hs + N - 1
            vs = (k-1)*3*N + (i-1)*3 + 1
            ve = vs + 2
            Rij[vs:ve,hs:he] = displacement(X[:,i,k],X[:,:,k])
            append!(rij,distance(Rij[vs:ve,hs:he]))
        end 
    end 
    Rij = (1/(L*N))*Rij
    Rij,rij
end 

function construct_rbf_bases(R)
    n = length(R)
    ϕ(r1,r2) = exp.(-((r1 .- r2)^2)/σ)
    bases = []
    for i=1:n
        push!(bases,(r,R)->ϕ(r,R[i]))
    end
    bases 
end 

function setup_regression(X,t,nb)
    R,r = displacement_and_distance(X,t) 
    T = transition_matrix(r) 
    T = T + 1e-8 * randn(size(T))
    T = T[:,1:nb]
    A = Array(R*T)
    b = reshape(V,(prod(size(V)),))
    x = A\b
    bases = construct_rbf_bases(r[1:nb])

    x,r[1:nb],bases
end 

function inferred_ϕ(bases,a)
    ff = function(r,R)
        potential = 0
        for (i,basis) in enumerate(bases)
            potential += basis(r,R)*a[i]
        end 
        potential
    end 
    ff
end 

function compute_displacement(X)
    _,N = size(X)
    R = zeros(3,N*N)
    k = 1
    for i=1:N
        for j=1:N
            R[:,k] = X[:,i] - X[:,j]
            k += 1
        end 
    end 
    R
end 

function compute_distance(R)
    N = size(R,2)
    r = zeros(N)
    for i=1:N
        r[i] = norm(R[:,i]) 
    end 
    r
end 

# Molly struct for pairwise interaction
struct LJ_reconstructed <: Molly.PairwiseInteraction
    phi::Function
    N::Int
    
    function LJ_reconstructed()

    end
end 

function Molly.force(inter::LJ_reconstructed,dr,coord_i,coord_j,atom_i,atom_j,boundary)
    
end 


# Inferring
begin
    natoms = 7
    traw,Xraw = LJ(natoms)
    σ = 1e-3
    ntimesteps = 10   
    nbasis = 20

    t,X = traw[1:ntimesteps],Xraw[1:ntimesteps,:,:]
    V = velocities(X,t)
    X = permutedims(X[1:end-1,:,:],(3,2,1))
    a,rb,bases = setup_regression(X,t,nbasis)
    phi = inferred_ϕ(bases,a)
end

# encompassing function ϕ
r = 0.0:0.01:5.0
Φ = [phi(ri,rb) for ri in r]
plot(r,Φ,title="ϕ",xlabel="r",label="ϕ")
savefig("figures/ϕ.png")