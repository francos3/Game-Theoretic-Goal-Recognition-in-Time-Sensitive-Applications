///////////////////////////////////////////////////////////////////////////////
///  DijkstraShortestPathAlg.h
///  The implementation of Dijkstra algorithm to get the shortest path of 
///  a pair of vertices in a graph. 
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 5/30/2010
/// 
///  $Id: DijkstraShortestPathAlg.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///
///////////////////////////////////////////////////////////////////////////////

#pragma once

using namespace std;
#include <stack>
#include <algorithm>

class DijkstraShortestPathAlg
{
private: // members

    Graph* m_pDirectGraph;

    std::map<BaseVertex*, double> m_mpStartDistanceIndex; 
    std::map<BaseVertex*, BaseVertex*> m_mpPredecessorVertex; 
    
    std::map<int, std::set<int> > AG_Edges;
    std::vector<std::set<unsigned> > AGP_Edges;
    std::vector<map<int, std::set<int> > > AG_sg_edges;

    std::set<int> m_stDeterminedVertices;
    
    std::multiset<BaseVertex*, WeightLess<BaseVertex> > m_quCandidateVertices;
    BaseVertex* current_source;
    std::map<int,shared_ptr<BasePath> > m_mpPerimPaths;
    std::map<int,shared_ptr<BasePath> > m_mpFullPaths;
    std::set<int> forbidden_nodes;
    std::vector<BaseVertex*> AGpath;
    std::set<int> AGpathSet;
    std::set<int> AGseen;
    std::vector<pair<unsigned, unsigned> > prefixes;
    std::vector<set<unsigned> > paths_per_prefix;
    std::vector<std::vector<unsigned> > prefixes_per_path;
    std::map<pair<unsigned, unsigned>, set<unsigned> > dest_paths_per_prefix;
    map<unsigned, unsigned> DestToOrder;

public:
    BaseVertex* origin;
    set<int> Destinations;
    vector<int> Destinations_order;
    std::set<std::vector<BaseVertex*> > AGpaths;
    std::map<int, std::vector<BaseVertex*> > AGpaths2;
    std::set<pair<BaseVertex*,double> > AGnodes;
    std::set<pair<unsigned,double> > AGPnodes;
    std::set<pair<unsigned,double> > AGnodesID;
    std::vector<std::set<pair<BaseVertex*,double> > > AG_subgraph;
    std::vector<std::set<pair<unsigned, double> > > AGP_subgraph;
    std::set<pair<BaseVertex*,double> > AG_destinations;
    std::map<pair<int,double>,set<pair<int,double> > > AG_edges_bw;
    DijkstraShortestPathAlg(Graph* pGraph):m_pDirectGraph(pGraph){}
    ~DijkstraShortestPathAlg(void){clear();}
    DijkstraShortestPathAlg();

    void clear();
    
    void set_Destinations_Order(vector<int>* dests)
    {
        Destinations_order.assign(dests->begin(),dests->end());
    }

    std::shared_ptr<set<vector<BaseVertex* > > > getAGpaths(){
        auto ptr=make_shared<std::set<std::vector<BaseVertex*> > >(AGpaths);
        return ptr;
    }
    std::shared_ptr<map<int,vector<BaseVertex* > > > getAGpaths2(){
    
        auto ptr=make_shared<std::map<int,std::vector<BaseVertex*> > >(AGpaths2);
        return ptr;
    }
    //std::shared_ptr<set<pair<BaseVertex*,double> > > getAGNodes(){
    //	auto ptr=make_shared<std::set<pair<BaseVertex*,double> > >(AGnodes);
    //	return ptr;
    //}
    //std::shared_ptr<set<pair<int,double> > > getAGNodesID(){
    //	auto ptr=make_shared<std::set<pair<int,double> > >(AGnodesID);
    //	return ptr;
    //}
    std::shared_ptr<vector<set<pair<BaseVertex*,double> > > > getAGsubgraphs(){
        auto ptr=make_shared<vector<std::set<pair<BaseVertex*,double> > > >(AG_subgraph);
        return ptr;
    }
    std::shared_ptr<set<pair<BaseVertex*,double> > > getAGDestinationNodes(){
        auto ptr=make_shared<std::set<pair<BaseVertex*,double> > >(AG_destinations);
        return ptr;
    }
    std::shared_ptr<std::map<int,std::set<int> > > getAGEdges(){
        auto ptr=make_shared<std::map<int,std::set<int> > > (AG_Edges);
        return ptr;
    }
    std::shared_ptr<std::map<pair<int,double>,std::set<pair<int,double> > > > getAGEdgesBw(){
        auto ptr=make_shared<std::map<pair<int,double>,set<pair<int,double> > > > (AG_edges_bw);
        return ptr;
    }
    std::shared_ptr<vector<map<int,std::set<int> > > > getSGEdges(){
        auto ptr=make_shared<vector<map<int,set<int> > > > (AG_sg_edges);
        return ptr;
    }
    std::shared_ptr<vector<set<unsigned> > > getPrefixToPaths(){
        auto ptr = make_shared<std::vector<set<unsigned> > >(paths_per_prefix);
        return ptr;
    }
    std::shared_ptr<vector<vector<unsigned> > > getPathsToPrefix(){
        auto ptr = make_shared<std::vector<std::vector<unsigned> > >(prefixes_per_path);
        return ptr;
    }
    std::shared_ptr<vector<set<unsigned> > > getAGPEdges(){
        auto ptr = make_shared<std::vector<set<unsigned> > >(AGP_Edges);
        return ptr;
    }
    //std::shared_ptr<std::vector<unsigned> > linkPrefixToPaths();
    std::shared_ptr<vector<set<unsigned>>> linkPrefixToPaths();
    std::shared_ptr<map<pair<unsigned, unsigned>, set<unsigned> > > get_dest_paths_per_prefix(){
    //std::map<pair<unsigned, unsigned>, set<unsigned> > dest_paths_per_prefix;
        auto ptr = make_shared<std::map<pair<unsigned, unsigned>, set<unsigned>>>(dest_paths_per_prefix);
        return ptr;
    }

    BasePath* get_shortest_path(BaseVertex* source, BaseVertex* sink);
    std::shared_ptr<BasePath> get_shortest_path2(BaseVertex* source, BaseVertex* sink);

    void set_predecessor_vertex(BaseVertex* vt1, BaseVertex* vt2)
    {
        m_mpPredecessorVertex[vt1] = vt2;
    }

    double get_start_distance_at(BaseVertex* vertex)
    {
        return m_mpStartDistanceIndex.find(vertex)->second;
    }

    void set_start_distance_at(BaseVertex* vertex, double weight)
    {
        m_mpStartDistanceIndex[vertex] = weight;
    }

    void get_shortest_path_flower(BaseVertex* root)
    {
        determine_shortest_paths(NULL, root, false);
    }
    std::shared_ptr<std::vector<std::pair< unsigned, unsigned > > > get_prefixes(){
        auto ptr=make_shared<std::vector<pair<unsigned,unsigned> > >(prefixes);
        return ptr;
    }

    // The following two methods are prepared for the top-k shortest paths algorithm
    BasePath* update_cost_forward(BaseVertex* vertex);
    void correct_cost_backward(BaseVertex* vertex);
    void determine_all_shortest_paths( BaseVertex* source, float lambda_limit=INT_MAX);
        
    void get_perim_paths(std::map<int,shared_ptr<BasePath> > &current_perim_paths){current_perim_paths=m_mpPerimPaths;};
    void make_full_paths(std::map<int,shared_ptr<BasePath> > *FullOriginPaths);
    std::map<int,shared_ptr<BasePath> >* get_pt_paths(){return &m_mpPerimPaths;};
    void print_paths();
    void print_paths_costs();
    void get_perim_nodes(set<int> *perim_nodes,int lambda_limit);
    shared_ptr<BasePath> recover_shortest_full_path(int sink);
    shared_ptr<BasePath> recover_shortest_perim_path(int sink);
    void populate_final_paths(set<BasePath> *Final_paths);
    void dump_edges();
    void print_perim_paths();
    bool is_node_reachable(int sink);
    double get_cost(int sink);
    int get_perim_size();
    void improve_distances(unordered_set<int> *FO_nodes, vector<Link> *FO_edges,unordered_map<pair<int,int>,float,hash_pair> *best_d);
    void createAG( BaseVertex* source, BaseVertex* sink, float Budget=INT_MAX,bool acyclic=true,bool grandparent_check=false);
    double createAGYen( BaseVertex* source, BaseVertex* sink, float Budget=INT_MAX,bool writeToFile=true);
    void printAGFile();
    //void calculateMinPathStrategy();
    void populateAGnodes();
    void populateAGPnodes();
    double load_paths(vector<vector<unsigned>> input_paths);
    void print_prefix(unsigned pref){
        for (size_t i = 0; i <= prefixes[pref].second;i++)
            cout << AGpaths2[prefixes[pref].first][i]->getID()<<","<<flush;
    }
    unsigned get_prefix_ending(unsigned pref){
        size_t i = prefixes[pref].second;
        return AGpaths2[prefixes[pref].first][i]->getID();
    }
    void print_path(unsigned id){
        for (auto node : AGpaths2[id])
            cout<<node->getID()<<",";
    }
    size_t path_size(unsigned id){
        return AGpaths2[id].size();
    }
    float get_path_cost(std::vector<BaseVertex*> *path){
        float cost = 0;
        for (size_t i = 0; i < path->size() - 1; i++)
        {
            auto parent = (*path)[i];
            auto child = (*path)[i+1];
            //cout << "parent:" << parent->getID() << flush << endl;
            cost += m_pDirectGraph->get_edge_weight(parent, child);
        }
        //cost = cost * cost;
        return cost;
    }
    float get_path_edge_cost(unsigned path, int parent){
        auto parent_node = AGpaths2[path][parent];
        auto child_node = AGpaths2[path][parent+1];
        //cout << "parent:" << parent_node->getID() << ",child:," <<child_node->getID()<< endl;
        return m_pDirectGraph->get_edge_weight(parent_node, child_node);
    }
    int get_path_size(unsigned path){
        return AGpaths2[path].size();
    }
    float get_subpath_cost(unsigned path, int node)
    {
        float cost = 0;
        bool start_counting = AGpaths2[path][0]->getID() == node;
        for (size_t i = 0; i < AGpaths2[path].size() - 1; i++)
        {
            auto parent = AGpaths2[path][i];
            auto child = AGpaths2[path][i + 1];
            //cout << "parent:" << parent->getID() << flush << endl;
            if(!start_counting){
                if(child->getID()==node){
                    start_counting = true;
                    cost += m_pDirectGraph->get_edge_weight(parent, child);
                }
            }
            else{
                cost += m_pDirectGraph->get_edge_weight(parent, child);
            }

        }
        //cost = cost * cost;
        return cost;
    }

    float get_cost_to_dest(unsigned pref, unsigned path)
    {
        // cout << "hola cost_to_dest" << flush<<endl;
        float cost = 0;
        // size_t path = prefixes[pref].first;
        for (size_t i = prefixes[pref].second; i < AGpaths2[path].size() - 1; i++)
        {
            // cout << "i:" << i << flush << endl;
            auto parent = AGpaths2[path][i];
            // cout << "parent:" << parent->getID() << flush << endl;
            auto child = AGpaths2[path][i + 1];
            // cout << "child:" << child->getID() << flush << endl;
            cost += m_pDirectGraph->get_edge_weight(parent, child);
        }
        // cout << "\t\tcost from:"; print_prefix(pref); cout<<"for path:";
        // print_path(path);
        // cout << ", is:," << cost << endl;
        //cost = cost * cost;
        return cost;
    }
    float get_pref_cost(unsigned pref){
        //cout << "hola cost_to_dest" << flush<<endl;
        float cost = 0;
        size_t path = prefixes[pref].first;
        for (size_t i = 0; i<prefixes[pref].second; i++)
        {
            //cout << "i:" << i << flush << endl;
            auto parent = AGpaths2[path][i];
            //cout << "parent:" << parent->getID() << flush << endl;
            auto child = AGpaths2[path][i+1];
            //cout << "child:" << child->getID() << flush << endl;
            cost += m_pDirectGraph->get_edge_weight(parent, child);
        }
        cout << "\t\tcost to:"; print_prefix(pref); cout<<"for path:"; print_path(prefixes[pref].first);
        cout << ", is:," << cost << endl;
        //cost = cost * cost;
        return cost;
    }



template <typename T>
bool is_prefix(const std::vector<T>& prefix, const std::vector<T>& vec,size_t pref_size) {
    if (pref_size >= vec.size()) return false; // A prefix can't be longer than the vector itself.

    //for (size_t i = 0; i < pref_size; ++i) 
    for (size_t i = 0; i <= pref_size; i++) {
        //if (prefix[pref_size - 1 - i]->getID() != vec[vec.size() - 1 - i]->getID()) 
        //    return false;
        if(prefix[i]->getID()!=vec[i]->getID())
            return false;
    }
    return true;
    //return std::equal(prefix.begin(), prefix.end(), vec.begin());
}
    //void expand_AG(double Budget);
protected:
    bool stuck( BaseVertex* source, BaseVertex* sink);

    void determine_shortest_paths(BaseVertex* source, BaseVertex* sink, bool is_source2sink);

    void improve2vertex(BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit=INT_MAX);
    void improve2vertexKeepPaths( BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit );
    void add_full_paths(std::map<int,BasePath*> origin_optimal_paths);
};
void printMemoryUsage();
