using LinearAlgebra
using SparseArrays
using Plots
using WriteVTK
using Meshes
using Random
using CoordRefSystems
using Makie

struct Material
    epsilon::Float64    # Relative permittivity
    mu::Float64        # Relative permeability
    sigma::Float64     # Conductivity
end

struct Node
    x::Float64
    y::Float64
    id::Int
end

struct Element
    nodes::Vector{Int}  # Node indices for triangular element
    material::Material
    id::Int
end

struct Mesh
    nodes::Vector{Node}
    elements::Vector{Element}
    boundary_nodes::Vector{Int}
end

struct MaxwellSolver
    mesh::Mesh
    M_eps::SparseMatrixCSC{Float64, Int}  # Mass matrix for electric field
    M_mu::SparseMatrixCSC{Float64, Int}   # Mass matrix for magnetic field
    C::SparseMatrixCSC{Float64, Int}      # Curl matrix
    dt::Float64
    t::Float64
    E_z::Vector{Float64}                  # Electric field (z-component)
    H_x::Vector{Float64}                  # Magnetic field (x-component)
    H_y::Vector{Float64}                  # Magnetic field (y-component)
end