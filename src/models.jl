abstract type RegressionModels end

mutable struct Maggioni <: RegressionModels
    positions
    velocities
    basis
    distances

    function Maggioni(positions,velocities,basis)
        new(positions,velocities,basis,compute_distances(positions))
    end 

    function Maggioni(positions,basis)
        new(positions,compute_velocities(positions),basis,compute_distances(positions))
    end 
end 

function compute_forces!(ic,p,t,model,forces)
    @unpack distances,positions,basis = model
    s = size(distances) # ic,ts,n,n
    s1 = size(forces) # ic,n,paramters_dims,ts   
    num_particles = s[4]
    for i=1:num_particles
        forces += basis(distances[ic,t,p,i]) * (positions[ic,p,:,t] - positions[ic,i,:,t])
    end 
    forces /= num_particles
end 

function forcefield(model::Maggioni)
    @unpack positions,distances,basis,velocities = model
    ics,n,p,tsteps = size(positions)
    forces = zero(positions)
    for ic in 1:ics
        for tstep in 1:tsteps
            for particle in 1:particles
                compute_forces!(ic,particle,tstep,model,forces)
            end 
        end 
    end 
    return forces
end     

function Optim.optimize(model::Maggioni)
    
end

function parameters(model::Maggioni)
    getparams(model.basis)
end 

function loss(model::Maggioni)
    localloss = function(params)
        setparams!(params,model)
        norm(model.velocities - forcefield(model)[:,:,:,1:end-1])
    end 
    return localloss
end