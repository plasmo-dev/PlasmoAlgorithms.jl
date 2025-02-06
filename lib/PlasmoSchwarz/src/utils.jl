"""
    _find_boundary_edges(graph::OptiGraph, subgraphs::Vector{OptiGraph}) -> Vector{Vector{OptiEdge}}

Identify the boundary edges for a set of subgraphs within a given hypergraph. 
Boundary edges are those incident to nodes at the interface between the subgraphs 
and the rest of the graph.

Args:
  - `graph::OptiGraph`: The main graph containing the subgraphs and edges.
  - `subgraphs::Vector{OptiGraph}`: A vector of subgraphs for which boundary edges 
    need to be identified.

Returns:
    A vector of vectors, where each inner vector contains the boundary edges 
    associated with a specific subgraph.
"""
function _find_boundary_edges(graph::OptiGraph, subgraphs::Vector{OptiGraph})
    projection = hyper_projection(graph)
    boundary_edge_list = Vector{Vector{OptiEdge}}()
    for subgraph in subgraphs
        subgraph_nodes = all_nodes(subgraph)
        boundary_edges = Plasmo.incident_edges(projection, subgraph_nodes)
        push!(boundary_edge_list, boundary_edges)
    end
    return boundary_edge_list
end

"""
    _extract_constraints(subgraph_boundary_edges::Vector{Vector{OptiEdge}}) -> Vector{Vector{ConstraintRef}}

Retrieve the constraints associated with boundary edges for a set of subgraphs. 
For each subgraph, the constraints tied to its boundary edges are collected.

Args:
  - `subgraph_boundary_edges::Vector{Vector{OptiEdge}}`: A vector of vectors, 
    where each inner vector contains the boundary edges of a specific subgraph.

Returns:
    A vector of vectors, where each inner vector contains the constraints 
    associated with the boundary edges of a specific subgraph.
"""
function _extract_constraints(subgraph_boundary_edges::Vector{Vector{OptiEdge}})
    subgraph_boundary_constraints = Vector{Vector{ConstraintRef}}()
    for edge_vector in subgraph_boundary_edges
        edge_constraints = Vector{ConstraintRef}()
        for edge in edge_vector
            append!(edge_constraints, all_constraints(edge))
        end
        push!(subgraph_boundary_constraints, edge_constraints)
    end
    return subgraph_boundary_constraints
end

"""
    _is_hierarchical(graph::OptiGraph)

Check whether the given graph is hierarchical (i.e. contains nodes and subgraphs.)
"""
function _is_hierarchical(graph::OptiGraph)
    return !isempty(graph.subgraphs) && !isempty(graph.optinodes)
end

# TODO: utilities for helpful checks

# check that overlap is at least 1
function _check_overlap() end

# check whether partitions are contiguous
function _check_contiguous_partitions() end

# modify non_contiguous partitions to make them contiguous
function _fix_non_contiguous_partitions() end
