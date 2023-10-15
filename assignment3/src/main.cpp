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

void printAdjacencyMatrix(const Graph& g){
    int N = boost::num_vertices(g);
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
}

double computeClosenessCentrality(const Graph& g){
    std::vector<double> centrality(boost::num_vertices(g));
    double c = 0;
    for (int i = 0; i < centrality.size(); ++i){
        double ci = 0;
        for (int j = 0; j < centrality.size(); ++j){
            int di = boost::degree(j, g);
            ci = ci + di;
        }
        c = ci * 1/(centrality.size()-1); 
    }

    return c * 1/centrality.size();
}

Graph createRandomGraph(int n, int m){
    int N = 3; // Number of vertices
    int M = 3; // Number of edges
    Graph g(0);

    // Create a random number generator and seed it.
    unsigned int random_seed = 42;
    boost::mt19937 gen(random_seed);

    boost::generate_random_graph(g, N, M, gen, false);
    
    return g;
}


double estimate_pvalue_binomial(double x, int T, int N, int M){
    // x: Closeness centrality
    // T: number of repetitions. 
    
    int f = 0;
    for (int t = 0; t < T; t++){
        // produce a random network following the null hypothesis
        Graph g = createRandomGraph(N, M);
        // Calculate x_NH on that network
        double x_nh = computeClosenessCentrality(g);
        if (x_nh >= x){
           f ++;
        }
    }
    return f/(double)T; 
}

double estimate_pvalue_degree_sequence(std::vector<int> deg_sequence, int T){
    int f = 0;

    // Use handshaking lemma to find number of vertices
    int E = 0.5*std::accumulate(deg_sequence.begin(), deg_sequence.end(), 0);

    // We must ensure connectivity to recreate a graph that could be a linguistic network
    Graph g(deg_sequence.size());
    bool connected = true;
    int i = 0;
    /***
    while(connected and i < deg_sequence.size()){
        for (int j = 0; j < deg_sequence[i]; ++j){
            // Connect the highest degree d to the next d vertices
            boost::add_edge(g, i, 0);
        }
        i++; 
    }
    */
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
        std::vector<int> deg_sequence(N);
        for(int i = 0; i < deg_sequence.size(); ++i){
            deg_sequence[i] = boost::degree(i, g);
        }
        std::cout << argv[i] << ";" << std::endl <<
                     "Number of vertices (N): " << N << ";" << std::endl <<
                     "Number of edges (E): " << E << ";" << std::setprecision(7) << std::endl <<
                     "Mean degree (k): " << 2.*double(E)/double(N) << ";" << std::setprecision(7) << std::endl <<
                     "Density of edfes (delta): " << 2.*double(E)/double(N*(N-1)) << std::endl;

    }

    // ****
    // Estimation of p-value via MC method (fraom class slides, 9 of 29)
    // ****
    int N = 3;
    int E = 3;
    std::vector<int> deg_sequence = {2, 2, 1}; 

    int T = 1;
    // assert (1/(double)T < alpha); // Statistically significant
    
    // double p_val = estimate_pvalue_binomial(2, T, N, E);
    double p_val = estimate_pvalue_degree_sequence(deg_sequence, T);

    return 0;
}
