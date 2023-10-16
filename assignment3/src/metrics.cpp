#include "metrics.h"
#include <iostream>
#include <boost/graph/johnson_all_pairs_shortest.hpp>

double exactClosenessCentrality(const Graph& g){
	std::vector<std::vector<int>> D(boost::num_vertices(g), std::vector<int>(boost::num_vertices(g)));
	
	auto weight_map = boost::make_constant_property<Graph::edge_descriptor>(1);

	bool success = boost::johnson_all_pairs_shortest_paths(
					g,
					D, 
					boost::weight_map(weight_map)
				   );
	if(!success) {
		std::cerr << "This sould  not happenn!" << std::endl;
		return 0.0;
	}
	
	const int N = boost::num_vertices(g);
    double c = 0;
	#pragma omp parallel for reduction(+ : c)
	for (int i = 0; i < N; ++i){
        for (int j = i+1; j < N; ++j){
            c += 1.0 / (double)D[i][j];
        }
    }

    return 2.0 * c / (double)(N*(N-1));
}

double montecarloClosenessCentrality(const Graph& g, int T, double delta) {
    

	return  0.0;
}

std::pair<double,double> boundedClosenessCentrality(const Graph& g) {
    return std::make_pair<double,double>(0.0,0.0);
}