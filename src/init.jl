using LinearAlgebra
using SparseArrays
using Plots
using WriteVTK
using Meshes
using Random
using CoordRefSystems
using Makie
using GLMakie

struct Material
    epsilon::Float64    # Relative permittivity
    mu::Float64        # Relative permeability
    sigma::Float64     # Conductivity
end

struct Node{n}
    point::Meshes.Point{𝔼{n}}
end

struct Element{n}
    nodes::Vector{Node{n}}
end

function nodes(nodelist::Vector, connec::Meshes.Connectivity)
    res = Vector{Node}()
    for i in indices(connec)
        push!(res, nodelist[i])
    end
    return res
end

struct Mesh{n}
    domain_nodes::Vector{Node{n}}
    domain_elements::Vector{Element{n}}
    boundary_nodes::Vector{Node{n}}
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