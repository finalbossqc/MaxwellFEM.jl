"""
Time-explicit Finite Element Method solver for Maxwell's equations in linear media
Solves the time-dependent Maxwell equations:
∂E/∂t = (1/ε)(∇ × H - J)
∂H/∂t = -(1/μ)(∇ × E)
where E is electric field, H is magnetic field, ε is permittivity, μ is permeability, J is current density
"""

"""
Create a rectangular mesh with triangular elements
"""
function create_rectangular_mesh(nx::Int, ny::Int, Lx::Float64, Ly::Float64, material::Material)
    nodes = Node[]
    elements = Element[]
    
    # Create nodes
    node_id = 1
    for j in 1:ny+1
        for i in 1:nx+1
            x = (i-1) * Lx / nx
            y = (j-1) * Ly / ny
            push!(nodes, Node(x, y, node_id))
            node_id += 1
        end
    end
    
    # Create triangular elements
    element_id = 1
    for j in 1:ny
        for i in 1:nx
            # Lower triangle
            n1 = (j-1)*(nx+1) + i
            n2 = (j-1)*(nx+1) + i + 1
            n3 = j*(nx+1) + i
            push!(elements, Element([n1, n2, n3], material, element_id))
            element_id += 1
            
            # Upper triangle
            n1 = (j-1)*(nx+1) + i + 1
            n2 = j*(nx+1) + i + 1
            n3 = j*(nx+1) + i
            push!(elements, Element([n1, n2, n3], material, element_id))
            element_id += 1
        end
    end
    
    # Identify boundary nodes
    boundary_nodes = Int[]
    for node in nodes
        if node.x ≈ 0.0 || node.x ≈ Lx || node.y ≈ 0.0 || node.y ≈ Ly
            push!(boundary_nodes, node.id)
        end
    end
    
    return Mesh(nodes, elements, boundary_nodes)
end

"""
Compute element mass matrix for triangular element
"""
function compute_element_mass_matrix(element::Element, nodes::Vector{Node}, material_property::Float64)
    # Get node coordinates
    x = [nodes[element.nodes[i]].x for i in 1:3]
    y = [nodes[element.nodes[i]].y for i in 1:3]
    
    # Compute element area
    area = 0.5 * abs((x[2]-x[1])*(y[3]-y[1]) - (x[3]-x[1])*(y[2]-y[1]))
    
    # Mass matrix for linear triangular elements
    M_e = zeros(3, 3)
    for i in 1:3
        for j in 1:3
            if i == j
                M_e[i,j] = material_property * area / 6.0
            else
                M_e[i,j] = material_property * area / 12.0
            end
        end
    end
    
    return M_e
end

"""
Compute element curl matrix for triangular element
"""
function compute_element_curl_matrix(element::Element, nodes::Vector{Node})
    # Get node coordinates
    x = [nodes[element.nodes[i]].x for i in 1:3]
    y = [nodes[element.nodes[i]].y for i in 1:3]
    
    # Compute element area
    area = 0.5 * abs((x[2]-x[1])*(y[3]-y[1]) - (x[3]-x[1])*(y[2]-y[1]))
    
    # Shape function derivatives
    b = [y[2]-y[3], y[3]-y[1], y[1]-y[2]] / (2*area)
    c = [x[3]-x[2], x[1]-x[3], x[2]-x[1]] / (2*area)
    
    # Curl matrix (for 2D TM mode: curl of scalar field)
    C_ex = zeros(3, 3)  # ∂Ez/∂y term
    C_ey = zeros(3, 3)  # -∂Ez/∂x term
    
    for i in 1:3
        for j in 1:3
            C_ex[i,j] = b[j] * area / 3.0  # ∂N_j/∂y integrated over element
            C_ey[i,j] = -c[j] * area / 3.0  # -∂N_j/∂x integrated over element
        end
    end
    
    return C_ex, C_ey
end

"""
Assemble global matrices
"""
function assemble_matrices(mesh::Mesh)
    n_nodes = length(mesh.nodes)
    n_elements = length(mesh.elements)
    
    # Initialize sparse matrices
    I_eps = Int[]
    J_eps = Int[]
    V_eps = Float64[]
    
    I_mu = Int[]
    J_mu = Int[]
    V_mu = Float64[]
    
    I_cx = Int[]
    J_cx = Int[]
    V_cx = Float64[]
    
    I_cy = Int[]
    J_cy = Int[]
    V_cy = Float64[]
    
    # Assemble element contributions
    for element in mesh.elements
        # Mass matrices
        M_eps_e = compute_element_mass_matrix(element, mesh.nodes, element.material.epsilon)
        M_mu_e = compute_element_mass_matrix(element, mesh.nodes, element.material.mu)
        
        # Curl matrices
        C_ex_e, C_ey_e = compute_element_curl_matrix(element, mesh.nodes)
        
        # Add to global matrices
        for i in 1:3
            for j in 1:3
                global_i = element.nodes[i]
                global_j = element.nodes[j]
                
                # Mass matrices
                push!(I_eps, global_i)
                push!(J_eps, global_j)
                push!(V_eps, M_eps_e[i,j])
                
                push!(I_mu, global_i)
                push!(J_mu, global_j)
                push!(V_mu, M_mu_e[i,j])
                
                # Curl matrices
                push!(I_cx, global_i)
                push!(J_cx, global_j)
                push!(V_cx, C_ex_e[i,j])
                
                push!(I_cy, global_i)
                push!(J_cy, global_j)
                push!(V_cy, C_ey_e[i,j])
            end
        end
    end
    
    # Create sparse matrices
    M_eps = sparse(I_eps, J_eps, V_eps, n_nodes, n_nodes)
    M_mu = sparse(I_mu, J_mu, V_mu, n_nodes, n_nodes)
    C_x = sparse(I_cx, J_cx, V_cx, n_nodes, n_nodes)
    C_y = sparse(I_cy, J_cy, V_cy, n_nodes, n_nodes)
    
    return M_eps, M_mu, C_x, C_y
end

"""
Initialize Maxwell solver
"""
function initialize_maxwell_solver(mesh::Mesh, dt::Float64)
    M_eps, M_mu, C_x, C_y = assemble_matrices(mesh)
    
    n_nodes = length(mesh.nodes)
    E_z = zeros(n_nodes)
    H_x = zeros(n_nodes)
    H_y = zeros(n_nodes)
    
    # Combine curl matrices for convenience
    C = [C_x C_y]
    
    return MaxwellSolver(mesh, M_eps, M_mu, C, dt, 0.0, E_z, H_x, H_y)
end

"""
Apply boundary conditions (Perfect Electric Conductor - PEC)
"""
function apply_boundary_conditions!(solver::MaxwellSolver)
    # Set electric field to zero on boundary nodes (PEC boundary)
    for node_id in solver.mesh.boundary_nodes
        solver.E_z[node_id] = 0.0
    end
end

"""
Add current source (Gaussian pulse)
"""
function add_current_source!(solver::MaxwellSolver, source_x::Float64, source_y::Float64, 
                           amplitude::Float64, t0::Float64, sigma::Float64)
    # Find closest node to source location
    min_dist = Inf
    source_node = 1
    for (i, node) in enumerate(solver.mesh.nodes)
        dist = sqrt((node.x - source_x)^2 + (node.y - source_y)^2)
        if dist < min_dist
            min_dist = dist
            source_node = i
        end
    end
    
    # Gaussian pulse current source
    J_z = amplitude * exp(-((solver.t - t0)/sigma)^2)
    
    # Add source term to the node
    solver.E_z[source_node] += solver.dt * J_z / solver.mesh.elements[1].material.epsilon
end

"""
Time step using explicit scheme
"""
function time_step!(solver::MaxwellSolver)
    n_nodes = length(solver.mesh.nodes)
    
    # Update magnetic field: H^{n+1/2} = H^{n-1/2} - dt/μ * curl(E^n)
    # For 2D TM mode: Hx = -∂Ez/∂y, Hy = ∂Ez/∂x
    
    # Compute curl of E_z
    curl_E_x = zeros(n_nodes)  # -∂Ez/∂y
    curl_E_y = zeros(n_nodes)  # ∂Ez/∂x
    
    for element in solver.mesh.elements
        C_ex_e, C_ey_e = compute_element_curl_matrix(element, solver.mesh.nodes)
        
        for i in 1:3
            global_i = element.nodes[i]
            for j in 1:3
                global_j = element.nodes[j]
                curl_E_x[global_i] += C_ex_e[i,j] * solver.E_z[global_j]
                curl_E_y[global_i] += C_ey_e[i,j] * solver.E_z[global_j]
            end
        end
    end
    
    # Update H field
    M_mu_inv = inv(Matrix(solver.M_mu))  # In practice, use iterative solver
    dH_x = -solver.dt * M_mu_inv * curl_E_x
    dH_y = -solver.dt * M_mu_inv * curl_E_y
    
    solver.H_x .+= dH_x
    solver.H_y .+= dH_y
    
    # Update electric field: E^{n+1} = E^n + dt/ε * (curl(H^{n+1/2}) - J^{n+1/2})
    # For 2D TM mode: curl(H) = ∂Hy/∂x - ∂Hx/∂y
    
    curl_H_z = zeros(n_nodes)
    for element in solver.mesh.elements
        C_ex_e, C_ey_e = compute_element_curl_matrix(element, solver.mesh.nodes)
        
        for i in 1:3
            global_i = element.nodes[i]
            for j in 1:3
                global_j = element.nodes[j]
                # curl(H)_z = ∂Hy/∂x - ∂Hx/∂y
                curl_H_z[global_i] += -C_ey_e[i,j] * solver.H_y[global_j] - C_ex_e[i,j] * solver.H_x[global_j]
            end
        end
    end
    
    # Update E field
    M_eps_inv = inv(Matrix(solver.M_eps))  # In practice, use iterative solver
    dE_z = solver.dt * M_eps_inv * curl_H_z
    
    solver.E_z .+= dE_z
    
    # Apply boundary conditions
    apply_boundary_conditions!(solver)
    
    # Update time
    solver.t += solver.dt
end

"""
Run Maxwell simulation
"""
function run_maxwell_simulation(nx::Int, ny::Int, Lx::Float64, Ly::Float64, 
                               T_final::Float64, dt::Float64; 
                               source_x::Float64=Lx/2, source_y::Float64=Ly/2,
                               save_interval::Int=10)
    # Create material (vacuum)
    material = Material(1.0, 1.0, 0.0)  # ε₀, μ₀, σ=0
    
    # Create mesh
    mesh = create_rectangular_mesh(nx, ny, Lx, Ly, material)
    
    # Initialize solver
    solver = initialize_maxwell_solver(mesh, dt)
    
    # Time stepping
    n_steps = Int(T_final / dt)
    
    # Storage for visualization
    E_history = []
    times = []
    
    println("Starting Maxwell simulation...")
    println("Mesh: $(nx)×$(ny) elements, $(length(mesh.nodes)) nodes")
    println("Time steps: $n_steps, dt = $dt")
    
    for step in 1:n_steps
        # Add current source (Gaussian pulse)
        if step * dt < 3.0  # Pulse duration
            add_current_source!(solver, source_x, source_y, 1.0, 1.0, 0.5)
        end
        
        # Time step
        time_step!(solver)
        
        # Save data for visualization
        if step % save_interval == 0
            push!(E_history, copy(solver.E_z))
            push!(times, solver.t)
            
            # Print progress
            if step % (n_steps ÷ 10) == 0
                println("Progress: $(100*step÷n_steps)%, t = $(round(solver.t, digits=3))")
            end
        end
    end
    
    return solver, E_history, times
end

"""
Visualize results
"""
function visualize_results(solver::MaxwellSolver, E_history, times)
    # Create coordinate arrays for plotting
    x_coords = [node.x for node in solver.mesh.nodes]
    y_coords = [node.y for node in solver.mesh.nodes]
    
    # Plot final electric field
    scatter(x_coords, y_coords, c=solver.E_z, ms=2, 
           title="Electric Field Ez at t=$(round(solver.t, digits=3))",
           xlabel="x", ylabel="y", colorbar_title="Ez")
    
    # Create animation frames
    println("Creating animation frames...")
    anim = @animate for i in 1:length(E_history)
        scatter(x_coords, y_coords, c=E_history[i], ms=2, 
               title="Electric Field Ez at t=$(round(times[i], digits=3))",
               xlabel="x", ylabel="y", colorbar_title="Ez",
               clims=(-maximum(abs.(E_history[end])), maximum(abs.(E_history[end]))))
    end
    
    return anim
end

# Example usage
function main()
    # Parameters
    nx, ny = 40, 40           # Mesh resolution
    Lx, Ly = 2.0, 2.0        # Domain size
    T_final = 5.0            # Final time
    dt = 0.001               # Time step (must satisfy CFL condition)
    
    # Run simulation
    solver, E_history, times = run_maxwell_simulation(nx, ny, Lx, Ly, T_final, dt)
    
    # Visualize results
    anim = visualize_results(solver, E_history, times)
    
    # Save animation
    gif(anim, "maxwell_simulation.gif", fps=10)
    println("Animation saved as maxwell_simulation.gif")
    
    return solver, E_history, times
end

# Run the simulation
# Uncomment the line below to run the simulation
# main()