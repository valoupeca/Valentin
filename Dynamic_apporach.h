//
// Created by lamure on 27/06/18.
//

#ifndef VALENTIN_DYNAMIC_APPORACH_H
#define VALENTIN_DYNAMIC_APPORACH_H


#include <string>

#include <fstream>

#include <utility>

#include <iostream>

#include <cstring>

#include <vector>

#include <boost/thread.hpp>

#include <boost/thread/thread.hpp>

#include <boost/thread/pthread/condition_variable.hpp>

#include <boost/graph/graphml.hpp>

#include <boost/config.hpp>

#include <boost/iterator/iterator_facade.hpp>

#include <boost/graph/graph_traits.hpp>

#include <boost/graph/properties.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/named_function_params.hpp>

#include <boost/graph/graph_utility.hpp>

#include <boost/property_map/property_map_iterator.hpp>

#include <boost/property_map/property_map.hpp>

#include <boost/graph/lookup_edge.hpp>

#include <boost/graph/random.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/graph/copy.hpp>


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "armadillo"

#include "ot.h"
#include "emd.h"

using namespace std;

using namespace arma;

using namespace boost;



struct vertex_info {

    std::string name;

    std::string label;

    string country;

    string rir;

    long prefixnum;

    long prefixall;

    long ftime;

    long ctime;

    float size;

    string asn;

    int r,g,b;

    float x,y;

    long count;

    int id;

};

/*  edges -> correspond aux arcs*/

struct edge_info {

    double weight;

    long count;

    long ftime;

    long ctime;

};

typedef adjacency_list <boost::listS, boost::vecS, boost::undirectedS,vertex_info,

        edge_info >  Graph_t;

typedef adjacency_list <boost::listS, boost::vecS, boost::directedS,vertex_info,

        edge_info >  Graph_t1;



typedef graph_traits < Graph_t >::vertex_iterator VertexIterator;

typedef graph_traits < Graph_t >::edge_iterator EdgeIterator;

typedef graph_traits < Graph_t >::adjacency_iterator AdjacencyIterator;

typedef graph_traits < Graph_t >::vertex_descriptor vertex_descriptor;

typedef graph_traits < Graph_t >::edge_descriptor edge_descriptor;

typedef property_map < Graph_t, vertex_index_t >::type IndexMap;

typedef boost::iterator_property_map < double*, IndexMap, double, double& > DistanceMap;

/*property map*/

/*in & out edges*/

typedef graph_traits<Graph_t>::in_edge_iterator inEdge;

typedef graph_traits<Graph_t>::out_edge_iterator outEdge;

class My_visitor : boost::default_bfs_visitor{
protected:
    std::set<int> dests;
public:
    My_visitor( const std::set<int> destset)
            : dests(destset) {};

    void initialize_vertex(const vertex_descriptor &s, const Graph_t &g) const {}
    void discover_vertex(const vertex_descriptor &s, const Graph_t &g) const {}
    void examine_vertex(const vertex_descriptor &s, const Graph_t &g) const {}
    void examine_edge(const edge_descriptor &e, const Graph_t &g) const {}
    void edge_relaxed(const edge_descriptor &e, const Graph_t &g) const {}
    void edge_not_relaxed(const edge_descriptor &e, const Graph_t &g) const {}
    void finish_vertex(const vertex_descriptor &s, const Graph_t &g)  {
        dests.erase(s);
        if (dests.empty())
            throw(2);
    }
};


mat Dynamic_calcul(vector<std::string> landmark);

void calcul_emd(Graph_t &graph_t, vector<int> &subsetOfLandmarks, set<int> &setOfLandmarksNeighbors, set<int> &setofChangedAS,
                set<int> &setOfChangedNeighbors, mat &matrix, mat &curvMat, int indexTab);
void calcul_ot(Graph_t &graph_t, vector<int> &subsetOfLandmarks, set<int> &setOfLandmarksNeighbors, set<int> &setofChangedAS,
               set<int> &setOfChangedNeighbors, mat &matrix, mat &curvMat, int indexTab);
void calcul_chemin(Graph_t graph_t, std::vector<int> sources, std::set<int> dests, mat &mat, int indexTab  );


void readGraphMLFile(Graph_t& graphToBuild, string gmlFileToRead);

#endif //VALENTIN_DYNAMIC_APPORACH_H
