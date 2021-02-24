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

class DijkstraShortestPathAlg
{
private: // members

	Graph* m_pDirectGraph;

	std::map<BaseVertex*, double> m_mpStartDistanceIndex; 
	std::map<BaseVertex*, BaseVertex*> m_mpPredecessorVertex; 

	std::set<int> m_stDeterminedVertices;
	
	std::multiset<BaseVertex*, WeightLess<BaseVertex> > m_quCandidateVertices;
	BaseVertex* current_source;
	std::map<int,shared_ptr<BasePath> > m_mpPerimPaths;
	std::map<int,shared_ptr<BasePath> > m_mpFullPaths;
	std::set<int> forbidden_nodes;
	
public:
	DijkstraShortestPathAlg(Graph* pGraph):m_pDirectGraph(pGraph){}
	~DijkstraShortestPathAlg(void){clear();}
	DijkstraShortestPathAlg();

	void clear();

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
protected:

	void determine_shortest_paths(BaseVertex* source, BaseVertex* sink, bool is_source2sink);

	void improve2vertex(BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit=INT_MAX);
	void improve2vertexKeepPaths( BaseVertex* cur_vertex_pt, bool is_source2sink, float lambda_limit );
	void add_full_paths(std::map<int,BasePath*> origin_optimal_paths);

};
