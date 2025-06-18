function generate_random_points(N::Int64, Lx::Float64, Ly::Float64; rng=Random.default_rng())
    xs = rand(rng, N) .* Lx
    ys = rand(rng, N) .* Ly

    points = [Point2(x, y) for (x, y) in zip(xs, ys)]

    return points
end

function generate_random_mesh_2D(N::Int, Lx::Float64, Ly::Float64, material::Material; method=DelaunayTesselation(), seed=Random.default_rng())
    points = generate_random_points(N, Lx, Ly; seed=seed)
    mesh = tesselate(points, method)
    return mesh
end