#include "random_graph.h"
#include <random>
#include <boost/random.hpp>
#include <boost/graph/random.hpp>

#include <vector>
#include <set>
#include <utility>

#include <chrono>

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
	//typedef boost::edge_list<Graph::edge> GraphEdge;
	
	std::vector<std::pair<int,int>> edges;
	std::set<std::pair<int,int>> contained;
    
	Graph::edge_iterator edge_it, edge_end;
    std::tie(edge_it,edge_end) = boost::edges(graph);
	for(edge_it; edge_it != edge_end; edge_it++){
		edges.push_back({edge_it->m_source, edge_it->m_target});
		contained.insert({edge_it->m_source, edge_it->m_target});
	}

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> edgeDist(0, boost::num_edges(graph) - 1);
    std::uniform_int_distribution<int> bitDist(0, 1);

	//std::chrono::steady_clock::time_point begin, end;
	//double time_edge = 0., time_quality = 0., time_swap = 0.;

    for(int it = 0; it < Q*boost::num_edges(graph); it++){
		//begin = std::chrono::steady_clock::now();
        // Choose two random edges
        int edge1 = edgeDist(gen);
        int edge2 = edgeDist(gen);
        
		if(edge1 == edge2)
			continue;

		
		int u = edges[edge1].first, v = edges[edge1].second,
			s = edges[edge2].first, t = edges[edge2].second;
		
		bool kind_swap = bitDist(gen);
		
		if(kind_swap) {
			// Swap so u and s contains the smaller index
			if(u > t){
				std::tie(t,u) = std::tie(u,t);
			}
			if(s > v) {
				std::tie(v,s) = std::tie(s,v);
			}
		} else {
			// Swap so u and s contains the smaller index
			if(u > s){
				std::tie(s,u) = std::tie(u,s);
			}
			if(v > t) {
				std::tie(t,v) = std::tie(v,t);
			}
		}

		//end = std::chrono::steady_clock::now();
		//time_edge += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

		//begin = std::chrono::steady_clock::now();
		if(kind_swap) {
			// Discard to avoid self-loops
    	    if( u == t || s == t) {
	            continue;
	        }

			// Discard to avoid creating double edge
        	if(contained.count(std::pair<int,int>{u,t}) == 1 || contained.count(std::pair<int,int>{s,v}) == 1) { 
    	        continue;
	        }
		} else {
			// Discard to avoid self-loops
    	    if( u == s || v == t) { 
	            continue;
	        }

			// Discard to avoid creating double edge
        	if(contained.count(std::pair<int,int>{u,s}) == 1 || contained.count(std::pair<int,int>{v,t}) == 1) { 
    	        continue;
	        }
		}
		
		//end = std::chrono::steady_clock::now();
		//time_quality += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
        
		//begin = std::chrono::steady_clock::now();
        // Swap the selected vertices with their neighbors
		contained.erase(std::pair<int,int>{u,v});
		contained.erase(std::pair<int,int>{s,t});
		if(kind_swap) {
			contained.insert(std::pair<int,int>{u,t});
			contained.insert(std::pair<int,int>{s,v});
	        edges[edge1].first = u; edges[edge1].second = t;
    	    edges[edge2].first = s; edges[edge2].second = v;
		} else {
			contained.insert(std::pair<int,int>{u,s});
			contained.insert(std::pair<int,int>{v,t});
	        edges[edge1].first = u; edges[edge1].second = s;
    	    edges[edge2].first = v; edges[edge2].second = t;
		}
		//end = std::chrono::steady_clock::now();
		//time_swap += std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();

		//if(it%1000 == 0) {
		//	std::cout << time_edge/it << " " << time_quality/it << " " << time_swap/it << std::endl;
		//}
    }

	Graph g(boost::num_vertices(graph));
	for(auto edge : edges)
		boost::add_edge(edge.first,edge.second,g);
    return g;
}
