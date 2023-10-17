#include "random_graph.h"
#include <random>
#include <boost/random.hpp>
#include <boost/graph/random.hpp>

Graph createRandomBinomialGraph(int n, int m){
    Graph g(0);

    // Create a random number generator and seed it.
	// Seed with a real random value, if available
	// Use static variables in order to avoid initialization
	// every time, because it is very expensive
    static std::random_device r;
    static boost::mt19937 gen(r());

    boost::generate_random_graph(g, n, m, gen, false);
    
    return g;
}

Graph createSwitchingModel(const Graph& graph, int Q){
    Graph g(graph);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> edgeDist(0, boost::num_edges(graph) - 1);

    for(int it = 0; it < Q*boost::num_edges(graph); it++){
        // Choose two random edges
        int edge1 = edgeDist(gen);
        int edge2 = edgeDist(gen);
        
        Graph::edge_iterator edge_it1, edge_it2, edge_end;
        std::tie(edge_it1,edge_end) = boost::edges(g);
        std::tie(edge_it2,edge_end) = boost::edges(g);

        for(int i = 0; i < edge1; i++)
            edge_it1++;
        for(int i = 0; i < edge2; i++)
            edge_it2++;
        
        int u = edge_it1->m_source, v = edge_it1->m_target;
        int s = edge_it2->m_source, t = edge_it2->m_target;

        if( u == t || s == t) { // Discard to avoid self-loops
            continue;
        }

        if(boost::edge(u,t,g).second || boost::edge(s,v,g).second) { // Discard to avoid creating double edge
            continue;
        }
        
        // Swap the selected vertices with their neighbors
        boost::remove_edge(u, v, g);
        boost::remove_edge(s, t, g);
        boost::add_edge(u, t, g);
        boost::add_edge(s, v, g);
    }

    return g;
}