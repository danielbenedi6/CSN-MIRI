#include <iostream>
#include <assert.h>
#include <iomanip>
#include <numeric>
#include <fstream>
#include <omp.h>
#include "defines.h"
#include "input.h"
#include "random_graph.h"
#include "metrics.h"
#include <random>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/visitors.hpp>

const std::string header_table_1 = "\\begin{table}[!htp]\n\\centering\n\\resizebox{\\columnwidth}{!}{\n\\begin{tabular}{lllll}\nLanguage & N & E & $\\langle k \\rangle$ & $\\delta$ \\\\ \\hline\n";
const std::string header_table_2 = "\\begin{table}[!htp]\n\\centering\n\\resizebox{\\columnwidth}{!}{\n\\begin{tabular}{llll}\nLanguage & Metric & $p$-value (binomial) & $p$-value (switching) \\\\ \\hline\n";

const std::string end_table_1 = "\\end{tabular}}\n\\label{tab:summary}\n\\caption{Summary of the properties of the degree sequences.}\n\\end{table}";
const std::string end_table_2 = "\\end{tabular}}\n\\label{tab:hypothesis}\n\\caption{Estimation of the $p$-values of the hypothesis $\\mathcal{C}_{NH} \\geq \\mathcal{C}$}\n\\end{table}";

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
double estimate_pvalue_binomial(double x, int T, int N, int M, std::string filename){
    std::vector<double> X_NH(T, 0.);
    
    int f = 0;
    #pragma omp parallel for reduction( + : f)
    for (int t = 0; t < T; t++){
        // produce a random network following the null hypothesis
        Graph g = createRandomBinomialGraph(N, M);
        // Calculate x_NH on that network
        double x_nh = montecarloClosenessCentrality(g, 1, 0.1);
		//std::cout << x_nh << " " << std::flush;

        X_NH[t] = x_nh;

        f += x_nh >= x;
    }

    std::ofstream out(filename + ".csv");
    for(int t = 0; t < T; t++)
        out << X_NH[t] << std::endl;
    out.close();

	//std::cout << std::endl;
    return (double)f/(double)T; 
}

/**
 * @brief Estimates the p-value of the closeness centrality x
 *        using T repetitions of a graph using the switching
 *        model given the degree sequence
 * 
 * @param graph     Input graph for degree sequence
 * @param x         Closeness Centrality of the hypothesis
 * @param Q         Number of repetitions of the switching model
 * @param T         Number of repetitions of the montecarlo estimation
 * @return double 
 */
double estimate_pvalue_switching(const Graph& graph, double x, int Q, int T, std::string filename){
    std::vector<double> X_NH(T, 0.);

    int f = 0;
    #pragma omp parallel for reduction( + : f)
    for (int t = 0; t < T; t++){
        // produce a random network following the null hypothesis
        Graph g = createSwitchingModel(graph, Q);
        double x_nh = montecarloClosenessCentrality(g, 1, 0.1);
		//double x_nh = exactClosenessCentrality(g);
		//std::cout << "Hypothesis value" << std::endl;

        X_NH[t] = x_nh;

        f += x_nh >= x;
    }

    std::ofstream out(filename + ".csv");
    for(int t = 0; t < T; t++)
        out << X_NH[t] << std::endl;
    out.close();

    return (double)f/(double)T; 
}

int main(int argc, char *argv[]) {
    if(argc < 2) {
        std::cerr << "Wrong number of parameters" << std::endl;
        std::cerr << "Usage: " << argv[0] << " input_file [input_file ...]" << std::endl;
        return 1;
    }

    std::ofstream table1("./tables/summary.tex");
    std::ofstream table2("./tables/hypothesis.tex");

    if(!table1.is_open() || !table2.is_open()) {
        std::cerr << "Could not open files to write tables. Check that folder \"tables\" exists." << std::endl;
        return 1;
    }

    table1 << header_table_1;
    table2 << header_table_2;

    for(int i = 1; i < argc; i++){
        Graph g;
        std::string file = std::string(argv[i]);
        if(!read_file(file, g)) {
            return 1;
        }
        const int N = boost::num_vertices(g);
        const int E = boost::num_edges(g);

        std::string language = "";
        size_t pos = 0;
        while ((pos = file.find("/")) != std::string::npos) {
            language = file.substr(0, pos);
            file.erase(0, pos + 1);
        }
        language = file;
        language = language.substr(0, language.find("_"));

        // Graph is ready to use :D
        #ifdef DEBUG
        std::cout << language << std::endl <<
                     "Number of vertices (N): " << N << std::endl <<
                     "Number of edges (E): " << E << std::endl <<
                     "Mean degree (k): " << std::setprecision(7) << 2.*double(E)/double(N) << std::endl <<
                     "Density of edges (delta): " << std::setprecision(7) << 2.*double(E)/double(N*(N-1)) << std::endl;
        #endif

        table1 << language << " & "
               << N << " & "
               << E << " & "
               << std::setprecision(5) << 2.*double(E)/double(N) << " & "
               << std::setprecision(5) << 2.*double(E)/double(N*(N-1));
        if(i < argc-1) 
            table1 << " \\\\ ";
        table1 << std::endl;
		
        double C = montecarloClosenessCentrality(g, 5, 0.1);
		double p_val_bin = estimate_pvalue_binomial(C, 100, N, E, language + "_binomial");
		double p_val_sw = estimate_pvalue_switching(g, C, 1+(int)std::log10(E), 100, language + "_switching");

        #ifdef DEBUG
		double exactC = exactClosenessCentrality(g);

		std::cout << "Montecarlo Closeness Centrality: " << C << std::endl;
		std::cout << "Closeness Centrality: " << exactC << std::endl;
		std::cout << "p-value (binomial): " << p_val_bin << std::endl;
		std::cout << "p-value (switching): " << p_val_sw << std::endl;
        #endif

        table2 << language << " & "
               << std::setprecision(5) << C << " & "
               << std::setprecision(5) << p_val_bin << " & "
               << std::setprecision(5) << p_val_sw;
        if(i < argc-1) 
            table2 << " \\\\ ";
        table2 << std::endl;
    }

    table1 << end_table_1;
    table2 << end_table_2;

    table1.close();
    table2.close();

    return 0;
}
