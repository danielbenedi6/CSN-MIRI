#pragma once

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/graph_traits.hpp"

typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;

typedef boost::adjacency_list<
    boost::listS,           // OutEdgeList: use listS for a sparse graph
    boost::vecS,            // VertexList: use vecS for vertex storage
    boost::undirectedS,     // Undirected graph
    boost::no_property,     // Vertex properties
    boost::no_property      // Edge properties. If we need weighet, it could be EdgeWeightProperty
> Graph;

