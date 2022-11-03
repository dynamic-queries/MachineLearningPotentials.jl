function compute_pairwise_distance(x1,x2)
    return norm(x1-x2)
end 

function compute_distances(X::Array)
    ics,n,p,tsteps = size(X)
    working_array = zeros(n,n)
    distances  = zeros(ics,tsteps,n,n)
    for ic in 1:ics
        for t in 1:tsteps
            for i=1:n
                for j=1:n
                    working_array[i,j] = compute_pairwise_distance(X[ic,i,:,t],X[ic,j,:,t])
                end 
            end
            distances[ic,t,:,:] .= working_array 
        end 
    end
    distances
end 

function compute_velocities(positions::Array)

end 