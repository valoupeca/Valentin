//
// Created by lamure on 27/06/18.
//

#include "Dynamic_apporach.h"


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



    }

}

void calcul_ot(Graph_t &graph_t, vector<int> &subsetOfLandmarks, set<int> &setOfLandmarksNeighbors, set<int> &setofChangedAS,
               set<int> &setOfChangedNeighbors, mat &matrix, mat &curvMat, int indexTab){
    int dist;
    int i;
    mat costMat;
    OTResult results;
    int I=0,J=0;
    std::pair<AdjacencyIterator, AdjacencyIterator> neighbors;
    for(auto l : subsetOfLandmarks) {
        i=0;
        I=0;
        uvec indexOfLandmarkNeighbors(degree(vertex(l,graph_t), graph_t)+1);
        vec sourceDist(degree(vertex(l,graph_t), graph_t)+1);
        sourceDist(i)=0.5;
        indexOfLandmarkNeighbors(i++) = distance(setOfLandmarksNeighbors.begin(), setOfLandmarksNeighbors.find(l));
        neighbors= boost::adjacent_vertices(vertex(l,graph_t), graph_t);
        for(; neighbors.first != neighbors.second; ++neighbors.first) {
            sourceDist(i)=0.5/(sourceDist.n_rows-1);
            indexOfLandmarkNeighbors(i++) =  distance(setOfLandmarksNeighbors.begin(), setOfLandmarksNeighbors.find(*neighbors.first));
        }
        for (auto k : setofChangedAS) {
            i=0;
            uvec indexOfChangedASNeighbors(degree(vertex(k,graph_t), graph_t)+1);
            rowvec destDist(degree(vertex(k,graph_t), graph_t)+1);
            destDist(i)=0.5;
            indexOfChangedASNeighbors(i++) = distance(setOfChangedNeighbors.begin(), setOfChangedNeighbors.find(k));
            neighbors= boost::adjacent_vertices(vertex(k,graph_t), graph_t);
            for(; neighbors.first != neighbors.second; ++neighbors.first) {
                destDist(i)=0.5/(destDist.n_cols-1);
                indexOfChangedASNeighbors(i++) =  distance(setOfChangedNeighbors.begin(), setOfChangedNeighbors.find(*neighbors.first));
            }
            costMat=matrix(indexOfChangedASNeighbors,indexOfLandmarkNeighbors);
            results = sinkhorn_knopp(destDist, sourceDist, costMat, 8e-2, 1.9 , 1e-6, 1e4);
            curvMat(I,J+indexTab)=1-results.optCost/costMat(0,0);
//            cout<<"Opt. Cost from "<<I<< " to "<<J<<" is "<<curvMat(I,J)<<endl;
            I++;
        }
        J++;
    }
}


void calcul_emd(Graph_t &graph_t, vector<int> &subsetOfLandmarks, set<int> &setOfLandmarksNeighbors, set<int> &setofChangedAS,
                set<int> &setOfChangedNeighbors, mat &matrix, mat &curvMat, int indexTab)
{


    int dist;
    int i;
    mat costMat;
    double results;
    int I=0,J=0;


    // double* G;

    double* alpha(0);
    double* beta(0);
    double* cost(0);
    std::pair<AdjacencyIterator, AdjacencyIterator> neighbors;
    for(auto l : subsetOfLandmarks) {
        i=0;
        I=0;
        uvec indexOfLandmarkNeighbors(degree(vertex(l,graph_t), graph_t)+1);
        vec sourceDist(degree(vertex(l,graph_t), graph_t)+1);
        sourceDist(i)=0.5;
        indexOfLandmarkNeighbors(i++) = distance(setOfLandmarksNeighbors.begin(), setOfLandmarksNeighbors.find(l));
        neighbors= boost::adjacent_vertices(vertex(l,graph_t), graph_t);
        for(; neighbors.first != neighbors.second; ++neighbors.first) {
            sourceDist(i)=0.5/(sourceDist.n_rows-1);
            indexOfLandmarkNeighbors(i++) =  distance(setOfLandmarksNeighbors.begin(), setOfLandmarksNeighbors.find(*neighbors.first));
        }
        for (auto k : setofChangedAS) {
            i=0;
            uvec indexOfChangedASNeighbors(degree(vertex(k,graph_t), graph_t)+1);
            rowvec destDist(degree(vertex(k,graph_t), graph_t)+1);
            destDist(i)=0.5;
            indexOfChangedASNeighbors(i++) = distance(setOfChangedNeighbors.begin(), setOfChangedNeighbors.find(k));
            neighbors= boost::adjacent_vertices(vertex(k,graph_t), graph_t);
            for(; neighbors.first != neighbors.second; ++neighbors.first) {
                destDist(i)=0.5/(destDist.n_cols-1);
                indexOfChangedASNeighbors(i++) =  distance(setOfChangedNeighbors.begin(), setOfChangedNeighbors.find(*neighbors.first));
            }
            costMat=matrix(indexOfChangedASNeighbors,indexOfLandmarkNeighbors);


            //curvMat(I,J)=results;

            //    curvMat(I,J+indexTab)=1-results/costMat(0,0);
//            cout<<"Opt. Cost from "<<I<< " to "<<J<<" is "<<curvMat(I,J)<<endl;



            alpha =  0;
            beta = 0;
            cost= 0;

            int maxIter = 150;

            int n1 = (costMat.n_rows);
            int n2 = (costMat.n_cols);


            //cout << "cosT MAT " << costMat.memptr() << endl;

            // curvMat.print(std::cout);
            results =  EMD_wrap(n1,n2, destDist.memptr(), sourceDist.memptr() ,costMat.memptr() , curvMat.memptr() , alpha, beta, maxIter);

            //   cout << " cost emd : " << results << endl;
            curvMat(I,J+indexTab)=1-results/costMat[0,0];
            cout<<I<<","<<J<<","<<curvMat(I,J+indexTab)<<endl;
//          cout<<"Opt. Cost from "<<I<< " to "<<J<<" is "<<curvMat(I,J)<<endl;
            I++;
        }
        J++;
    }

}
void readGraphMLFile ( Graph_t1& designG, string fileName )

{
    boost::dynamic_properties dp;

    dp.property("count", get(&edge_info::count, designG));

    dp.property("ftime", get(&edge_info::ftime, designG));

    dp.property("ctime", get(&edge_info::ctime, designG));

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

    ifstream gmlStream;

    gmlStream.open(fileName.c_str(), ifstream::in);

    boost::read_graphml(gmlStream, designG, dp);

    gmlStream.close();

}               /* -----  end of method ExpressionGraph::readGraphMLFile  ----- */



mat Dynamic_calcul(vector<std::string> landmark)
{



        cout << "Armadillo version: " << arma_version::as_string() << endl;


        string file = "/home/lamure/Documents/developpement/test.graphml";

        Graph_t1 graph_t1;

        ifstream is(file.c_str());

        cout<<"DEBUT de LECTURE "<<endl;

        readGraphMLFile(graph_t1, file);

        cout<<"FICHIER LU! "<<endl;

        Graph_t graph_t;

        copy_graph(graph_t1, graph_t);
        graph_t1.clear();

        int N = num_vertices(graph_t);

        /* v2 , prendre noeud arrivée et de départ, chercher les plus court voisins et développé la matrice OT*/



        thread *threads[landmark.size()];

        int timeWindowLow = 1515117900;

        int timeWindowSize = 60;

        int numOfThreads = thread::hardware_concurrency();
        cout<< "Num of possible thread :"<<numOfThreads<<endl;
        numOfThreads = 7;

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


        // Now let's calculate the curvature matrix
        mat curvMat1=mat(setofChangedAS.size(),setOfLandmarks.size());
        cout<<"fin de calcul de la matrice des distances"<<endl;
        landmarks.clear();
        threadNum = 0;
        numOfElementsPerThread = setOfLandmarks.size()/numOfThreads+1;
        count = 1;
        threadNum =0;
        for (auto l :setOfLandmarks){
            landmarks.push_back(l);
            if (count % numOfElementsPerThread ==0 ){
                threads[threadNum] = new thread(calcul_ot, std::ref(graph_t),  landmarks , setOfLandmarksNeighbors,
                                                setofChangedAS, setOfChangedNeighbors, std::ref(dsp), std::ref(curvMat1),
                                                threadNum*numOfElementsPerThread);
                cout<<"Thread Num:"<<threadNum<<endl;
                threadNum ++;
                landmarks.clear();
            }
            count ++;
        }
        for(int c = 0; c < threadNum; c++)
        {
            threads[c]->join();
            delete threads[c];
        }


        // Now let's calculate the curvature matrix with optimal transport

        mat curvMat=mat(setofChangedAS.size(),setOfLandmarks.size());

        cout<<"fin de calcul de la matrice des distances"<<endl;


        cout << " calcul EMD " << endl;

        landmarks.clear();
        numOfThreads=1;
        numOfElementsPerThread = setOfLandmarks.size()/numOfThreads;
        count = 1;
        threadNum =0;
        for (auto l :setOfLandmarks){
            landmarks.push_back(l);
            if (count % numOfElementsPerThread ==0 ){
                threads[threadNum] = new thread(calcul_emd, std::ref(graph_t),  landmarks , setOfLandmarksNeighbors,
                                                setofChangedAS, setOfChangedNeighbors, std::ref(dsp), std::ref(curvMat),
                                                threadNum*numOfElementsPerThread);
                cout<<"Thread Num:"<<threadNum<<endl;
                threadNum ++;
                landmarks.clear();
            }
            count ++;


        }
        for(int c = 0; c < threadNum; c++)
        {
            threads[c]->join();
            delete threads[c];
        }

        /* vector<int> setoflands(setOfLandmarks.begin(),setOfLandmarks.end());
         int result = 0;
         calcul_emd((graph_t),  setoflands , setOfLandmarksNeighbors,
                 setofChangedAS, setOfChangedNeighbors, (dsp), (curvMat),
                 150);*/

        cout << "fin du calcul  = "  << endl;





        mat difference=mat(setofChangedAS.size(),setOfLandmarks.size());

        difference = curvMat - curvMat1;

        difference.print(std::cout);

        /* difference des matrices */



    return difference;

}