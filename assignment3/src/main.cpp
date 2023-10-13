#include <iostream>
#include <vector>
#include <assert.h>
#include <iomanip>
#include "defines.h"
#include "input.h"
#include <boost/graph/random.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/random.hpp>
#include <boost/graph/closeness_centrality.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graphviz.hpp>

// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> Graph;
const double alpha = 0.05;

double compute_nh(){
    double C = 0;
    int N = 3; // Number of vertices
    int M = 3; // Number of edges
    Graph g(0);

    // Create a random number generator and seed it.
    unsigned int random_seed = 42;
    boost::mt19937 gen(random_seed);

    std::vector<int> component(N);
    boost::generate_random_graph(g, N, M, gen, false);

    // Sanity check: Print the adjacency matrix.
    // boost::write_graphviz(std::cout, g); 
    
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (boost::edge(i, j, g).second) {
                std::cout << "1 ";
            } else {
                std::cout << "0 ";
            }
        }
        std::cout << std::endl;
   }


    int num_components = boost::connected_components(g, &component[0]);
    assert(num_components == 1);
    std::vector<double> Ci(N);
    
    
    boost::closeness_centrality(g, 
                                boost::make_iterator_property_map(Ci.begin(), 
                                    boost::get(boost::vertex_index, g))
                                );
    /*
    double Ci;
    double C = 0;
    for (int i=0; i <= N ; i++){
        Ci = 1/(N-1);
        C = C+Ci;
    }
    */

    // Iterate through the vertices and print their closeness centrality.
    Graph::vertex_iterator vi, vi_end;
    for (boost::tie(vi, vi_end) = boost::vertices(g); vi != vi_end; ++vi){ 
        C = C + Ci[*vi];
        std::cout << "Vertex " << *vi << " - Closeness Centrality: " << Ci[*vi] << std::endl;
    }
    return C/(double)N;
}

double estimate_pvalue(double x, int T){
    std::cout << "In estimate p-value funciton" << std::endl;
    // x: Closeness centrality
    // T: number of repetitions. 
    
    // assert (1/(double)T < alpha); // Statistically significant
    int f = 0;
    for (int t=0; t <= T; t++){
        // produce a random network following the null hypothesis
        // Calculate x_NH on that network
        double x_nh = compute_nh();
        if (x_nh >= x){
           f ++;
        }
    }
    return f/(double)T; 
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        std::cerr << "Wrong number of parameters" << std::endl;
        std::cerr << "Usage: " << argv[0] << " input_file [input_file ...]" << std::endl;
        return 1;
    }

    for(int i = 1; i < argc; i++){
        Graph g;
        if(!read_file(std::string(argv[i]), g)) {
            return 1;
        }

        // Graph is ready to use :D

        // First print the language summary. 
        // Maybe write to file?
        int N = boost::num_vertices(g),
            E = boost::num_edges(g);
        std::cout << argv[i] << ";" << std::endl <<
                     "Number of vertices (N): " << N << ";" << std::endl <<
                     "Number of edges (E): " << E << ";" << std::setprecision(7) << std::endl <<
                     "Mean degree (k): " << 2.*double(E)/double(N) << ";" << std::setprecision(7) << std::endl <<
                     "Density of edfes (delta): " << 2.*double(E)/double(N*(N-1)) << std::endl;

    }

    // ****
    // Estimation of p-value via MC method (fraom class slides, 9 of 29)
    // ****
    int T=5;
    double p_val = estimate_pvalue(2, T);

    return 0;
}
