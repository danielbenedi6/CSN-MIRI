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
