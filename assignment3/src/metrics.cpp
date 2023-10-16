#include "metrics.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <random>
#include <unordered_map>
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
    const int Nprime = std::max(2,(int)(delta*boost::num_vertices(g)));
    double C = 0.;

    // Create a random number generator and seed it.
    // Seed with a real random value, if available
    // Use static variables in order to avoid initialization
    // every time, because it is very expensive
    static std::random_device r;
    static std::mt19937 gen(r());

    std::vector<int> vertices(boost::num_vertices(g));
    boost::integer_range<int> range(0, boost::num_vertices(g));
    std::copy(range.begin(), range.end(), vertices.begin());

    for(int rep = 0; rep < T; rep++) {
        double Cprime = 0.;
        std::shuffle(vertices.begin(), vertices.end(), gen);

        // Initialize should_save so the first Nprime vertices are saved
        std::vector<bool> should_save(boost::num_vertices(g), false);
        for(int j = 0; j < Nprime; j++)
            should_save[vertices[j]] = true;
        for(int i = 0; i < Nprime; i++) {
            // From 0 to i should not be saved because we can multiply by two
            // So we can simply set the i-th to false
            should_save[vertices[i]] = false;

            std::queue<std::pair<int,int>> Q;
            Q.push({vertices[i],0});
            std::unordered_map<int,int> distance;
            std::vector<bool> visited(boost::num_vertices(g), false);
            while(distance.size() != Nprime-1 && !Q.empty()) {
                int vertex, time;
                std::tie(vertex, time) = Q.front();
                Q.pop();

                auto iters = boost::adjacent_vertices(vertex, g);
                for(auto it = iters.first; it != iters.second; it++) {
                    if(!visited[*it]) {
                        Q.push({*it, time+1});
                        if(should_save[*it]) {
                            distance.insert({*it, time+1});
                        }
                        visited[*it] = true;
                    }
                }
            }
            
            for(auto vertex_time : distance) {
                Cprime += 2.0 / (double)vertex_time.second;
            }
        }

        Cprime /= (double)(Nprime*(Nprime-1));

        C += Cprime;
    }

	return  C/(double)T;
}

std::pair<double,double> boundedClosenessCentrality(const Graph& g) {
    return std::make_pair<double,double>(0.0,0.0);
}