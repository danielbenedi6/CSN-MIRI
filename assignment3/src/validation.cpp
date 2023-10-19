#include <iostream>
#include <fstream>
#include "defines.h"
#include "input.h"
#include "random_graph.h"
#include "metrics.h"

const int REP = 2;
const std::string language = "Catalan";
const std::string filename = "./data/"+language+"_syntactic_dependency_network.txt";

int main(int argc, char *argv[]) {
    std::ofstream table("./validation.csv");

    if(!table.is_open()) {
        std::cerr << "Could not open files to write validation tables." << std::endl;
        return 1;
    }

    table << "Type;Model;Closeness" << std::endl;;
	Graph g;
	if(!read_file(filename, g) ) {
		return 1;
	}
	double exactC = exactClosenessCentrality(g);
	table << "Exact;"+language+";" << exactC << std::endl;

	std::vector<double> C(REP, 0.);
	#pragma omp parallel for
	for(int i = 0; i < REP; i++) {
        double c = montecarloClosenessCentrality(g, 1, 0.1);
		C[i] = c;
	}
	for(int i = 0; i < REP; i++) {
		table << "Montecarlo;"+language+";" << C[i] << std::endl;
	}
    
	const int N = boost::num_vertices(g);
    const int E = boost::num_edges(g);

    Graph g_bin = createRandomBinomialGraph(N, E);
	exactC = exactClosenessCentrality(g_bin);
	table << "Exact;Binomial;" << exactC << std::endl;

	#pragma omp parallel for
	for(int i = 0; i < REP; i++) {
        double c = montecarloClosenessCentrality(g_bin, 1, 0.1);
		C[i] = c;
	}
	for(int i = 0; i < REP; i++) {
		table << "Montecarlo;Binomial;" << C[i] << std::endl;
	}
    
    Graph g_swi = createSwitchingModel(g, 1 + (int)std::log10(E));
	exactC = exactClosenessCentrality(g_swi);
	table << "Exact;Switching;" << exactC << std::endl;

	#pragma omp parallel for
	for(int i = 0; i < REP; i++) {
        double c = montecarloClosenessCentrality(g_swi, 1, 0.1);
		C[i] = c;
	}
	for(int i = 0; i < REP; i++) {
		table << "Montecarlo;Switching;" << C[i] << std::endl;
	}

    table.close();

    return 0;
}
