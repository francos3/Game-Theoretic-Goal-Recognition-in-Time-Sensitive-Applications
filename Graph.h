///////////////////////////////////////////////////////////////////////////////
///  Graph.h
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 8/18/2010
///  modified by Santiago Franco @date 17/06/2019
/// 
///  $Id: Graph.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////
#ifndef GRAPH
#define GRAPH

#include <iostream>
#include <memory>
#include <algorithm>
#include <iterator>
#include <bits/stdc++.h> 
#include "main.h"
#include "perfecthash_paths.hpp"



extern bool debug;
extern std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1l, 1000000000l> > > start_time;
extern std::chrono::time_point<std::chrono::_V2::system_clock, std::chrono::duration<long int, std::ratio<1l, 1000000000l> > > orig_start_time;
extern int N;
extern bool min_cal;
extern bool random_destinations;
extern int random_positions;
extern int random_seed;
extern float YenTime;
extern std::vector<int> OrigDests;
extern float avg_dest_intradist;
extern float avg_dest_orig_dest;
extern int iter_limit;
//extern int limit_hits;

struct Link
{
	int u;
	int v;
	double weight;
};

//typedef Link link;

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
  if ( !v.empty() ) {
    out << '[';
    std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ","));
    out << "]";
  }
  return out;
}




#pragma once

using namespace std;
//using namespace boost::timer;
//extern cpu_timer g_timer;

class Path : public BasePath
{
public: 

  Path(const vector<BaseVertex*>& vertex_list, double weight):BasePath(vertex_list,weight){}

	// display the content
	void PrintOut(ostream& out_stream) const
	{
		out_stream << "Cost: " << m_dWeight << " Length: " << m_vtVertexList.size() << endl;
		for(vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); pos!=m_vtVertexList.end();++pos)
		{
			out_stream << (*pos)->getID() << " ";
		}
		out_stream << endl <<  "*********************************************" << endl;	
	}

};
class PathMatrix
{
  vector<set<vector<int>,vect_comp_int > > Path_Matrix;
  vector<set<vector<int>,vect_comp_int > > Candidate_Path_Matrix;
  map<vector<int>,vector<set<vector<int>,vect_comp_int > >,vect_comp_int > cover_map;
  unordered_map<vector<int>,size_t,VectorHash> path_ids;
  unordered_map<int,unordered_map<vector<int>,size_t,VectorHash>::iterator > active_columns;
  unordered_map<int,unordered_map<vector<int>,size_t,VectorHash>::iterator > active_rows;
  vector<vector<int> > cover_matrix;
  vector<int> destinations;
  set<int> safe_nodes;
  map<int,int> safe_nodes2;
  vector<vector<vector<int> > > path_T;
  vector<vector<vector<int> > > path_L;//L for lambda
  vector<int> indiv_path_values_L;
  vector<int> removed_paths;
  double avg_lambda=0;
  int min_lambda=INT_MAX;
  int max_lambda=0;
  int old_max_lambda=INT_MAX;
  int lambda_points=0;
  set<int> W;//unobserved nodes 
  int C=-1;
  int iteration=0;
  int populated_lambda_values=0;
  bool min_cal=false;
  float MatrixTime=0;
  float InitialLambda=0;
  int final_paths=0;
  int InitialPaths=0;
  public:
  PathMatrix(){};
  
  void set_W(set<int>  *_W);
  void set_C(int _C);
  void add_paths_set(shared_ptr<set<vector<int>,vect_comp_int > > destination_paths,int dest);
  void eliminate_paths_C();
  void print_paths();
  //void get_cost(int node);
  void record_initial_paths_gperf();
  void record_initial_paths();
  void record_best_solution();
  void get_safe_nodes(int k);
  void get_safe_nodes2(int k,int min_dests);
  int calculate_lambda(int t, vector<int> path);
  bool erase_disclosing_paths(int k,int &disclosing_paths);
  void erase_false_undisclosing_paths(int k);
  void print_t_indexes(int i=0);
  void print_l_indexes(int i=0,bool full_print=false);
  void get_dependent_set_for_path(vector<int> cover_path);
  void populate_cover_matrix();
  void get_dependent_set_for_all_paths();
  void print_dependent_set();
  void elminate_max_lambda_paths_U();
  void eliminate_max_lambda_path_U(vector<int> elim_path);
  void set_path_values_from_cover_paths();
  void get_maximal_lambda_paths();
  void clear_lambda_values();
  int is_lambda_connected();
  int remove_paths_higher_lambda();
  void check_for_improvement_on_lambda_values();
  void print_Candidate_Path_Matrix();
  int get_number_paths_PM();
  int get_number_paths(vector<set<vector<int>,vect_comp_int > > *PM);
  int get_number_paths_lambda();
  void print_W();
  void print_final_results();
  bool vect_reverse_comparator(std::vector<int> v1,std::vector<int> v2) ;
  int get_max_lambda();
  int get_old_max_lambda();
  void set_max_lambda(int input);
  void set_old_max_lambda(int input);
  void set_InitialLambda(int input);
  void calculate_path_min_lambdas_from_cover_matrix();
  bool remove_maximal_lambda_paths();
  void eliminate_max_lambda_path_matrix(vector<int> next_path);
  void recalculate_Path_Matrix_from_lambdas();
  void set_min_cal(bool min_cal);
  void print_exit_statement();
  int get_remaining_dests();
};
class Graph
{
public: // members

	const static double DISCONNECT; 

	typedef set<BaseVertex*>::iterator VertexPtSetIterator;
	typedef map<BaseVertex*, set<BaseVertex*>*>::iterator BaseVertexPt2SetMapIterator;

protected: // members
	// Basic information
	map<BaseVertex*, set<BaseVertex*>*> m_mpFanoutVertices;
	map<BaseVertex*, set<BaseVertex*>*> m_mpFaninVertices;
	map<long, double> m_mpEdgeCodeWeight; 
	vector<BaseVertex*> m_vtVertices;
	int m_nEdgeNum;
	int m_nVertexNum;

	map<int, BaseVertex*> m_mpVertexIndex;

	// Members for graph modification
	set<int> m_stRemovedVertexIds;
	set<pair<int,int> > m_stRemovedEdge;

public:

	// Constructors and Destructor
	Graph(const string& file_name);
	Graph(const Graph& rGraph);
	Graph() { clear(); }
	~Graph(void);
	void dump_edges();

	void clear();

	void set_number_vertices( const int& nVertices );
	int get_number_vertices(){return m_nVertexNum;};
	void add_link( const int& s, const int& d, const double& weight );
	void setv();
	int getVertexNum();

	BaseVertex* get_vertex(int node_id);
	
	long get_edge_code(const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt) const;
	set<BaseVertex*>* get_vertex_set_pt(BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index);

	double get_original_edge_weight(const BaseVertex* source, const BaseVertex* sink);

	double get_edge_weight(const BaseVertex* source, const BaseVertex* sink);
	void get_adjacent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);
	void get_precedent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);

	/// Methods for changing graph
	void remove_edge(const pair<int,int> edge)
	{
		m_stRemovedEdge.insert(edge);
	}

	void remove_vertex(const int vertex_id)
	{
		m_stRemovedVertexIds.insert(vertex_id);
	}

	void recover_removed_edges()
	{
		m_stRemovedEdge.clear();
	}

	void recover_removed_vertices()
	{
		m_stRemovedVertexIds.clear();
	}

	void recover_removed_edge(const pair<int,int> edge)
	{
		m_stRemovedEdge.erase(m_stRemovedEdge.find(edge));
	}

	void recover_removed_vertex(int vertex_id)
	{
		m_stRemovedVertexIds.erase(m_stRemovedVertexIds.find(vertex_id));
	}
	
private:
	void _import_from_file(const string& file_name);

};
#endif
