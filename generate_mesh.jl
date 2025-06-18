function generate_random_points_2D(N::Int64, Lx::Float64, Ly::Float64; rng=Random.default_rng())
    xs = rand(rng, N) .* Lx
    ys = rand(rng, N) .* Ly

    println(xs)

    points = Meshes.Point.(xs, ys)
    
    return PointSet(points)
end

function generate_random_mesh_2D(N::Int, Lx::Float64, Ly::Float64, material::Material; method=DelaunayTesselation(), rng=Random.default_rng())
    points = generate_random_points_2D(N, Lx, Ly; rng=rng)
    mesh = tesselate(points, method)
    return mesh
end