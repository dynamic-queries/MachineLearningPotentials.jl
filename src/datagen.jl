function random_data()
    # Parameters for the simulation
    ics = 10    
    trajs = 50
    num_particles = 100
    p = 3
    # generate data
    data = rand(ics,num_particles,p,trajs)
    data
end 

function random_derivative(data)
    s = collect(size(data))
    s[3] -= 1
    rand(s...)
end     

function Base.Matrix(obj::Vector)
    n = length(obj)
    p = length(obj[1])
    X = zeros(n,p)
    for i=1:n
        for j=1:p
            X[i,j] = obj[i][j].val
        end 
    end 
    X
end 

function lj_homogenous()
    n_atoms = 100
    boundary = CubicBoundary(2.0u"nm", 2.0u"nm", 2.0u"nm")
    temp = 298.0u"K"
    atom_mass = 10.0u"u"

    atoms = [Atom(mass=atom_mass, σ=0.3u"nm", ϵ=0.2u"kJ * mol^-1") for i in 1:n_atoms]
    coords = place_atoms(n_atoms, boundary; min_dist=0.3u"nm")
    velocities = [velocity(atom_mass, temp) for i in 1:n_atoms]
    pairwise_inters = (LennardJones(),)
    simulator = VelocityVerlet(
        dt=0.002u"ps",
        coupling=AndersenThermostat(temp, 1.0u"ps"),
    )

    sys = System(
        atoms=atoms,
        pairwise_inters=pairwise_inters,
        coords=coords,
        velocities=velocities,
        boundary=boundary,
        loggers=(temp=TemperatureLogger(100),),
    )

    simulate!(sys, simulator, 10_000)
    Matrix(sys.coords)
end 