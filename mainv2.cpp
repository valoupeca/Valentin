//

// Created by lamure on 05/06/18.

//

/* git_code.cpp

                *

                *  Created on: 17 mai 2018

*      Author: lamure

*/

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



using namespace std;

using namespace arma;

using namespace boost;

//On définit les caractéristiques de nos noeuds, arc pour notre graphe

// vertex = correspond aux noeuds

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

/*def graph*/

typedef adjacency_list <boost::listS, boost::vecS, boost::undirectedS,vertex_info,

        edge_info >  Graph_t;

typedef adjacency_list <boost::listS, boost::vecS, boost::directedS,vertex_info,

        edge_info >  Graph_t1;

/* Declare dynamic properties for graphml format. */

boost::dynamic_properties dp;

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

/*class visitor pour djisktra*/
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

/*  Declare function prototype */

void readGraphMLFile(Graph_t& graphToBuild, string gmlFileToRead);

/*

 * ===  FUNCTION  ======================================================================

 *         Name:  readGraphMLFile

 *  Description:  Read a graphml file and update the graph designG.

 * =====================================================================================

 */

void readGraphMLFile ( Graph_t1& designG, string fileName )

{

    dp.property("id", get(&vertex_info::id, designG));

    dp.property("name", get(&vertex_info::name, designG));

    dp.property("label", get(&vertex_info::label, designG));

    dp.property("asn", get(&vertex_info::asn, designG));

    dp.property("rir", get(&vertex_info::rir, designG));

    dp.property("country", get(&vertex_info::country, designG));

    dp.property("prefixnum", get(&vertex_info::prefixnum, designG));

    dp.property("prefixall", get(&vertex_info::prefixall, designG));

    dp.property("ftime", get(&vertex_info::ftime, designG));

    dp.property("ctime", get(&vertex_info::ctime, designG));

    dp.property("size", get(&vertex_info::size, designG));

    dp.property("r", get(&vertex_info::r, designG));

    dp.property("g", get(&vertex_info::g, designG));

    dp.property("b", get(&vertex_info::b, designG));

    dp.property("x", get(&vertex_info::x, designG));

    dp.property("y", get(&vertex_info::y, designG));

    dp.property("count", get(&vertex_info::count, designG));

    dp.property("weight", get(&edge_info::weight, designG));

    dp.property("count", get(&edge_info::count, designG));

    dp.property("ftime", get(&edge_info::ftime, designG));

    dp.property("ctime", get(&edge_info::ctime, designG));

    ifstream gmlStream;

    gmlStream.open(fileName.c_str(), ifstream::in);

    boost::read_graphml(gmlStream, designG, dp);

    gmlStream.close();

}               /* -----  end of method ExpressionGraph::readGraphMLFile  ----- */

/*

 * ===  FUNCTION  ======================================================================

 *         Name:  main

 *  Description:

 * =====================================================================================

 */

void calcul_chemin(Graph_t graph_t, std::vector<int> sources, std::set<int> dests, mat &mat, int indexTab  ) {

    IndexMap index = get(vertex_index, graph_t);

    graph_traits<Graph_t>::vertex_iterator i, iend;

    int cnt=0;

    std::vector<vertex_descriptor> parent(num_vertices(graph_t));

//    std::vector<double> distance(num_vertices(graph_t));

    std::vector<double> distance(num_vertices(graph_t));

    My_visitor vis(dests);

    for (auto l : sources) {

        parent.clear();

        distance.clear();

        vertex_descriptor source = vertex(l, graph_t);

        auto weights = get(&edge_info::weight,graph_t);

        try {
            dijkstra_shortest_paths(graph_t, source,
                    //    boost::no_named_parameters().weight_map(weights));
                    //.weight_map(get(&edge_info::weight, graph_t)));
                                    boost::weight_map(boost::get(&edge_info::weight,graph_t))
                                            .distance_map(boost::make_iterator_property_map(distance.begin(), boost::get(boost::vertex_index,graph_t)))
                                            .predecessor_map(boost::make_iterator_property_map(parent.begin(), boost::get(boost::vertex_index,graph_t)))
                                            .visitor(vis)
            );
        }
        catch (int exception) {
            DistanceMap distanceMap(&distance[0], index);
            int cnt2 =0;
            for (auto d: dests){
                int v= vertex(d, graph_t);
                mat(cnt2,indexTab+cnt)=distanceMap[v];
                cnt2++;
            }
        }
        cnt++;


//        dijkstra_shortest_paths(graph_t, source, predecessor_map(&parent[0]).distance_map(&distance[0]));



        /*     typedef std::vector<Graph_t::edge_descriptor> PathType;

             PathType path;

             cout << "node : " << landmark << endl;

             vertex_descriptor v = landmark;

             cout << "name vertex landmark " << get(vertex_bundle, graph_t)[v].name << endl;

             for (vertex_descriptor u = parent[v];

                  u != v; // Keep tracking the path until we get to the source

                  v = u, u = parent[v]) // Set the current vertex to the current predecessor,     and the predecessor to one level up

             {

                 std::pair<Graph_t::edge_descriptor, bool> edgePair = boost::edge(u, v, graph_t);

                 Graph_t::edge_descriptor edge = edgePair.first;

                 path.push_back(edge);

             }

             // Write shortest path

             std::cout << "Shortest path from " << noeud_principal << " to " << v << "  "

                       << get(vertex_bundle, graph_t)[noeud_principal].name;

             cout << " et noeud d'arrivée " << get(vertex_bundle, graph_t)[landmark].name << std::endl;

             float totalDistance = 0;

             for (PathType::reverse_iterator pathIterator = path.rbegin(); pathIterator != path.rend(); ++pathIterator) {

                 std::cout << get(vertex_bundle, graph_t)[boost::source(*pathIterator, graph_t)].name << " -> " <<

                           get(vertex_bundle, graph_t)[boost::target(*pathIterator, graph_t)].name

                           << " = " << get(boost::edge_bundle, graph_t)[*pathIterator].weight << std::endl;

                 totalDistance = totalDistance + get(boost::edge_bundle, graph_t)[*pathIterator].weight;

             }

             mat[1][index_tab] = distanceMap[landmark];

             std::cout << "Distance: " << distanceMap[landmark] << std::endl;

          */

    }

}

int main(int argc, char** argv)

{

    cout << "Armadillo version: " << arma_version::as_string() << endl;


    string file = "//home/lamure/Documents/developpement/test.graphml";

    Graph_t1 graph_t1;

    ifstream is(file.c_str());

    cout<<"DEBUT de LECTURE "<<endl;

    readGraphMLFile(graph_t1, file);

    cout<<"FICHIER LU! "<<endl;

    Graph_t graph_t;

    copy_graph(graph_t1, graph_t);

    int N = num_vertices(graph_t);

    /* v2 , prendre noeud arrivée et de départ, chercher les plus court voisins et développé la matrice OT*/

    vector<std::string> landmark(20);

    thread *threads[landmark.size()];

    landmark[0] ="24482";

    landmark[1] ="131284";

    landmark[2] ="1455";

    landmark[3] ="531";

    landmark[4] ="1501";

    landmark[5] ="52519";

    landmark[6] ="263276";

    landmark[7] ="22773";

    landmark[8] ="12400";

    landmark[9] ="42926";

    landmark[10] ="266084";

    landmark[11] ="11663";

    landmark[12] ="393223";

    landmark[13] ="9549";

    landmark[14] ="35816";

    landmark[15] ="52376";

    landmark[16] ="266117";

    landmark[17] ="56660";

    landmark[18]="52752";

    int timeWindowLow = 1515117900;

    int timeWindowSize = 60;

    int numOfThreads = thread::hardware_concurrency();
    cout<< "Num of possible thread :"<<numOfThreads<<endl;

    std::map<std::string, int> asMap;

    IndexMap index = get(vertex_index,graph_t);

    Graph_t::vertex_descriptor v;

    std::pair<VertexIterator, VertexIterator> vp;

    std::set<int> setofChangedAS;

    int maxTime=0;

    for (vp = vertices(graph_t); vp.first != vp.second; ++vp.first) {

        v= *vp.first;

        asMap[graph_t[v].asn] = index[v];

        if ( (graph_t[v].ctime >= timeWindowLow) && (graph_t[v].ctime < timeWindowLow+timeWindowSize)) {

            setofChangedAS.insert(index[v]);

        }

    }

    cout<<"Num of changed AS:"<<setofChangedAS.size()<<endl;

    std::set<int> setOfLandmarks;

    for(int c = 0; c < landmark.size(); c++){

        setOfLandmarks.insert(asMap[landmark[c]]);

    }

    cout<<"Num of landmarks :"<<setOfLandmarks.size()<<endl;

    std::set<int> setOfLandmarksNeighbors;

    std::pair<AdjacencyIterator, AdjacencyIterator> neighbors;

    for(auto l : setOfLandmarks) {

        neighbors= boost::adjacent_vertices(vertex(l,graph_t), graph_t);

        setOfLandmarksNeighbors.insert(l);

        for(; neighbors.first != neighbors.second; ++neighbors.first) {

            setOfLandmarksNeighbors.insert(index[*neighbors.first]);

        }

    }

    cout<<"Size of setOfLandmarkNeighbors :"<<setOfLandmarksNeighbors.size()<<endl;

    std::set<int> setOfChangedNeighbors;

    for(auto l : setofChangedAS) {

        setOfChangedNeighbors.insert(l);

        neighbors= boost::adjacent_vertices(vertex(l,graph_t), graph_t);

        for(; neighbors.first != neighbors.second; ++neighbors.first) {

            setOfChangedNeighbors.insert(index[*neighbors.first]);

        }

    }

    cout<<"Size of setOfChangedNeighbors :"<<setOfChangedNeighbors.size()<<endl;

    // A la fin de cette partie nous avons les données suivante :

    // setOfLandmarks = ensemble des index des Landmarks

    // setOfLandmarksNeighbors = ensemble des index dex voisins des landmarks

    // setOfChangedAS = ensemble des index des AS qui ont change dans la fenetre

    // setofChangedNeighbors = ensemble des index des voisins des AS qui ont change dans la fenetre

    mat dsp(setOfChangedNeighbors.size(), setOfLandmarksNeighbors.size());

    //boost::numeric::ublas::matrix<double> mat(setOfChangedNeighbors.size(), setOfLandmarksNeighbors.size());

    int numOfElementsPerThread = setOfLandmarksNeighbors.size()/numOfThreads;

    std::vector<int> landmarks;
    int count = 1;
    int threadNum =0;
    for (auto l :setOfLandmarksNeighbors){
        landmarks.push_back(l);
        if (count % numOfElementsPerThread ==0 ){
            threads[threadNum] = new thread(calcul_chemin, std::ref(graph_t),  landmarks , setOfChangedNeighbors, std::ref(dsp) , threadNum*numOfElementsPerThread);
            threadNum ++;
            cout<<"Thread Num:"<<threadNum<<endl;
            landmarks.clear();
        }
        count ++;
    }

    for(int c = 0; c < threadNum; c++)
    {
        threads[c]->join();
        delete threads[c];
    }
    for (int i=0; i<20; i++){

        for (int j=0; j<20; j++){

            cout << dsp(i,j)<<" ";

        }

        cout <<endl;

    }


    cout << "sum(sum(mzt)) : " << sum(sum(dsp)) << endl;
    std::cout << std::endl;

}