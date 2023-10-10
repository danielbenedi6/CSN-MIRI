#include <iostream>
#include <iomanip>
#include "defines.h"
#include "input.h"

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
        std::cout << argv[i] << ";" << 
                     N << ";" <<
                     E << ";" << std::setprecision(7) <<
                     2.*double(E)/double(N) << ";" << std::setprecision(7) <<
                     2.*double(E)/double(N*(N-1)) << std::endl;

    }

    return 0;
}