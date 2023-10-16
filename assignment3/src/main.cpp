#include <iostream>
#include <assert.h>
#include <iomanip>
#include <numeric>
#include "defines.h"
#include "input.h"
#include "random_graph.h"
#include "metrics.h"

/**
 * @brief Prints the Adjacency Matrix of the graph
 * 
 * @param g Input graph
 */
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
    std::cout << std::endl;
}

/**
 * @brief Estimates the p-value of the closeness centrality x
 *        using T repetitions of a binomial graph with N vertices
 *        and M edges.
 * 
 * @param x         Closeness Centrality of the hypothesis
 * @param T         Number of repetitions of the montecarlo estimation
 * @param N         Number of vertices in the graph
 * @param M         Number of edges in the graph
 * @return double   p-value of the null-hypothesis x_NH >= x
 */
double estimate_pvalue_binomial(double x, int T, int N, int M){
    
    int f = 0;
    #pragma omp parallel for reduction( + : f)
    for (int t = 0; t < T; t++){
        // produce a random network following the null hypothesis
        Graph g = createRandomBinomialGraph(N, M);
        // Calculate x_NH on that network
        double x_nh = montecarloClosenessCentrality(g, 1, 0.1);
		//std::cout << x_nh << " " << std::flush;
        
        f += x_nh >= x;
    }
	//std::cout << std::endl;
    return (double)f/(double)T; 
}

/**
 * @brief Estimates the p-value of the closeness centrality x
 *        using T repetitions of a graph using the switching
 *        model given the degree sequence
 * 
 * @param deg_sequence  Degree sequence
 * @param T             Number of repetitions 
 * @return double       p-value of the null-hypothesis x_NH >= x
 */
double estimate_pvalue_degree_sequence(std::vector<int> deg_sequence, int T){
    int f = 0;

    // Use handshaking lemma to find number of vertices
    int E = std::accumulate(deg_sequence.begin(), deg_sequence.end(), 0)/2;

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
                     "Density of edges (delta): " << 2.*double(E)/double(N*(N-1)) << std::endl;
		
        double C = montecarloClosenessCentrality(g, 3, 0.1);
		std::cout << "Montecarlo Closeness Centrality: " << C << std::endl;
		//double C = exactClosenessCentrality(g);
		//std::cout << "Closeness Centrality: " << C << std::endl;

		double p_val_bin = estimate_pvalue_binomial(C, 10, N, E);
		std::cout << "p-value (binomial): " << p_val_bin << std::endl;
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
