CrystalNode = Union{Int64, Vector{Int64}}

struct CrystalGraph
    graph
    nodes
    incidence
end

@doc raw"""
    CrystalBase(B, num_elem, num_ops, e, f, eps, phi)

A crystal for quantum groups ``U_q(\mathfrak{g})`` of a Kac-Moody algebra ``\mathfrak g``.
Contains the necessary data to compute the corresponding crystal graph.
"""
struct CrystalBase
    B::Vector{CrystalNode} # must have unique values for each
    num_elem
    num_ops
    e # 
    f #
    eps # 
    phi #
    graph::Union{CrystalGraph, Missing} # TODO???
end

mutable struct WeylGroup
    cox_grp::CoxGrp
    simple_refls::Vector{CoxElt}
    elements
end

mutable struct RootSystem
    type::Symbol
    rank::Int64
    # calculated on demand
    vector_rep_crystal # ::CrystalBase
    gap_lalg
    gap_rs
    gap_weyl_group
    positive_roots
    cartan_matrix
    cartan_matrix_inv
    functional_to_root_matrix
    root_to_functional_matrix
    functional_to_fundamental_weight_matrix
    fundamental_weight_to_functional_matrix
    weyl_group
    character_ring
    uea
    uea_ordering
end

mutable struct Weight
    rs::RootSystem
    as_root
    as_fundamental_weight
    as_functional
end


mutable struct HighestWeightModule
    rs::RootSystem
    weight::Weight
    character
    dimension
    crystal
    crystal_graph
    ideal
    polytope
end

struct TensorProductWrapper
    rs::RootSystem
    maximal_vectors::Vector{Vector{Int64}}
    first_component_length
end

WeylGroupElement = Vector{T} where {T <: Integer}

mutable struct DemazureModule
    rs::RootSystem
    weight::Weight
    weyl_group_element
    character
    dimension
    crystal
    crystal_graph
end

struct WeightLatticeAlgebra
    lpr::LaurentPolynomialRing
end