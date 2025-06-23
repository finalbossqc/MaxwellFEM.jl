using Pkg
Pkg.activate(".")

include("../src/init.jl")
include("../src/generate_mesh.jl")

function test1()
    N = 100; Lx = 10.0; Ly = 10.0; material = Material(1.0, 1.0, 0.0);
    mesh = generate_2D_triangular(N, Lx, Ly, material)
    return mesh
end

test1()