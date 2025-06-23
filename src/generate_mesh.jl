function sample_domain_2D(N::Int64, Lx::Float64, Ly::Float64; rng=Random.default_rng())
    xs = rand(rng, N) .* Lx
    ys = rand(rng, N) .* Ly
    points = Meshes.Point.(xs, ys)
    return Vector{Meshes.Point}(points)
end

function sample_boundary_2D(N::Int64, Lx::Float64, Ly::Float64; rng=Random.default_rng())
    boundary_points = [Meshes.Point(point.coords.x.val, 0) for point in sample_domain_1D(Int(sqrt(N)), Lx)]
    append!(boundary_points, [Meshes.Point(point.coords.x.val, Ly) for point in sample_domain_1D(Int(sqrt(N)), Lx)])
    append!(boundary_points, [Meshes.Point(0, point.coords.x.val) for point in sample_domain_1D(Int(sqrt(N)), Ly)])
    append!(boundary_points, [Meshes.Point(Lx, point.coords.x.val) for point in sample_domain_1D(Int(sqrt(N)), Ly)])
    return Vector{Meshes.Point}(boundary_points)
end

function sample_domain_1D(N::Int64, L::Float64; rng=Random.default_rng())
    xs = rand(rng, N) .* L
    points = Meshes.Point.(xs)
    return Vector{Meshes.Point}(points)
end

function generate_2D_triangular(N::Int64, Lx::Float64, Ly::Float64, material::Material; method=DelaunayTesselation(), rng=Random.default_rng())
    points = sample_domain_2D(N, Lx, Ly; rng=rng)
    boundary_points = sample_boundary_2D(N, Lx, Ly; rng=rng)
    append!(points, boundary_points)
    mesh = tesselate(points, method)
    return Mesh{2,3}(mesh.vertices, mesh.topology.connec, boundary_points)
end