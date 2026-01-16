#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_set>    
#include <cstdlib>
#include <cstring>
#include <fstream>

#include "methods.h"
//#include "methods_HKPR.h"

int main(int argc, char *argv[])
{

    SFMTinitialize();

    string filename;
    //string algo = "RSFPPR";
    string algo;

    int num_query_nodes = 20;
    double l1_error = 1e-12;
    double eps = 0.1;
    int k = 20;
    //double reps = 1e-12;
    double reps;
    double lambda=0.05;
    int N0=500;
    //double lambda1=19.4437;
    //double lambda1=57.1334;

    int i=1;
    char *endptr;
    while(i<argc) {
        if(!strcmp(argv[i], "-d")) {
            filename = argv[++i];
        } else if(!strcmp(argv[i], "-algo")) {
            algo = argv[++i];
        } else if(!strcmp(argv[i], "-n")) {
            num_query_nodes = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-e")) {
            eps = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-k")) {
            k = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-re")) {
            reps = strtod(argv[++i], &endptr);
        } else if(!strcmp(argv[i], "-N0")) {
            N0 = strtod(argv[++i], &endptr);
        }
        i++;
    }

    graph g;
    // g.read_graph_weighted(filename);
    g.read_graph(filename);
    int landmark = 0;
    for (int i=0;i<g.n;i++)
    {
        if (g.degree[i] > g.degree[landmark])
        {
            landmark=i;
        }
    }
    cout<<"Landmark node: "<<landmark<<endl;

    if(algo == "GEN_QUERY") {
        ifstream infile("data/" + filename + ".query");
        if(infile.good()) {
            cerr << "query file already generated!" << endl;
            return 0;
        }
        ofstream outf("data/" + filename + ".query");
        generateQueryNode(2*num_query_nodes, g.n, outf);
        outf.close();
    } else if(algo == "groundtruth") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s, t;
            infile >> s;
            infile >> t;
            cout<< "source node: "<<s<< " sink nodes: "<<t<<endl;

            groundtruth(g, s, t, 1e-9, filename);

        }
        infile.close();
    } else if(algo == "lanczos") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s, t;
            infile >> s;
            infile >> t;
            cout<< "source node: "<<s<< " sink nodes: "<<t<<endl;
            lanczos(g, s, t, k, filename);
        }
        infile.close();
    } else if(algo == "lanczos_push") {
        ifstream infile("data/" + filename + ".query");
        for(int i=0; i<num_query_nodes; i++) {
            int s, t;
            infile >> s;
            infile >> t;
            cout<< "source node: "<<s<< " sink nodes: "<<t<<endl;
            lanczos_Local(g, s, t, eps, k, filename);
        }
        infile.close();
    } 
    
    else if(algo == "plot_results") {
        plot_results(g, num_query_nodes, filename);
    } 
    return 0;
}
